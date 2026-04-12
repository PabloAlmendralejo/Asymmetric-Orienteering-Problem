"""
Visualization of orienteering results on real terrain maps.

Generates:
  1. Map overview with DEM contours
  2. Cost surface components (ACR, HCR, elevation, slope, valid mask)
  3. Hillshade + elevation
  4. Directional asymmetry heatmap
  5. Optimal route overlaid on map (from solver JSON output)

Usage:
    python visualize_results.py torremocha
    python visualize_results.py la_muela
    python visualize_results.py torremocha --output-dir figures/
"""
import sys
import os
import json
import numpy as np
import math
import heapq
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from PIL import Image, ImageDraw

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python'))

from core.cost_functions import (
    build_base_cost_from_omap, cost_to_veg_proxy, minetti_factor, WT
)
from core.terrain_analysis import build_base_cost_with_hcr
from core.pathfinding import compute_slope_dir
from core.control_placement import create_valid_mask
import argparse


# ── Helpers ──

def load_terrain(cache_path):
    with np.load(cache_path, allow_pickle=True) as t:
        data = {
            'elev': t['elevation'], 'slope': t['slope'],
            'cost_omap': t['cost_omap'], 'img': t['img'],
            'bounds': tuple(t['bounds']), 'resolution': float(t['resolution']),
        }
        for key in ['path_grid', 'wall_grid', 'water_grid']:
            if key in t:
                data[key] = t[key]
    return data


def make_extent(bounds):
    return [bounds[0], bounds[2], bounds[1], bounds[3]]


def pixel_to_utm(px, py, bounds, resolution):
    return bounds[0] + px * resolution, bounds[3] - py * resolution


def cost_str(cost_units, resolution=2.0, base_speed=2.5, ref_weight=0.3):
    s = cost_units * resolution / (base_speed * ref_weight)
    if s < 60:
        return f"{s:.0f}s"
    elif s < 3600:
        return f"{int(s // 60)}m{int(s % 60):02d}s"
    else:
        return f"{int(s // 3600)}h{int((s % 3600) // 60):02d}m"


# ── A* path tracer (for route visualization) ──

NEIGHBORS = [(-1, 0, 1.0), (1, 0, 1.0), (0, -1, 1.0), (0, 1, 1.0),
             (-1, -1, 1.414), (1, -1, 1.414), (-1, 1, 1.414), (1, 1, 1.414)]


def astar_path(cost_grid, start, goal):
    """A* pathfinding on the display cost grid."""
    h, w = cost_grid.shape
    sx, sy = start
    gx, gy = goal
    if not (0 <= sx < w and 0 <= sy < h and 0 <= gx < w and 0 <= gy < h):
        return math.inf, []
    mc = max(float(cost_grid[cost_grid < 5].min()), 0.01) if (cost_grid < 5).any() else 0.05

    def H(x, y):
        dx, dy = abs(x - gx), abs(y - gy)
        return mc * (1.414 * min(dx, dy) + abs(dx - dy))

    dist = np.full((h, w), np.inf)
    dist[sy, sx] = 0.0
    vis = np.zeros((h, w), dtype=bool)
    par = np.full((h, w, 2), -1, dtype=np.int32)
    pq = [(H(sx, sy), 0.0, sx, sy)]

    while pq:
        f, g, x, y = heapq.heappop(pq)
        if vis[y, x]:
            continue
        vis[y, x] = True
        if x == gx and y == gy:
            path = []
            cx, cy = gx, gy
            while not (cx == sx and cy == sy):
                path.append((cx, cy))
                px, py = par[cy, cx]
                cx, cy = int(px), int(py)
            path.append((sx, sy))
            path.reverse()
            return g, path
        for ndx, ndy, sl in NEIGHBORS:
            nx, ny = x + ndx, y + ndy
            if 0 <= nx < w and 0 <= ny < h and not vis[ny, nx]:
                c = cost_grid[ny, nx]
                if not np.isfinite(c):
                    continue
                ng = g + sl * 0.5 * (cost_grid[y, x] + c)
                if ng < dist[ny, nx]:
                    dist[ny, nx] = ng
                    par[ny, nx] = [x, y]
                    heapq.heappush(pq, (ng + H(nx, ny), ng, nx, ny))
    return math.inf, []


def trace_route_legs(disp_cost, nodes, route):
    """Trace A* paths for each leg of the route."""
    legs = []
    seq = [0] + list(route) + [0]
    for i in range(len(seq) - 1):
        _, path = astar_path(disp_cost, nodes[seq[i]], nodes[seq[i + 1]])
        legs.append(path)
    return legs


# ── Plot functions ──

ROUTE_COLOR = '#0066FF'


def plot_map_overview(terrain, extent, output_dir, area_name):
    """Original map + DEM contours."""
    fig, axes = plt.subplots(1, 2, figsize=(20, 10))
    axes[0].imshow(terrain['img'], extent=extent, aspect='equal')
    axes[0].set_title(f'{area_name} — Orienteering Map')
    axes[0].set_xlabel('UTM Easting (m)')
    axes[0].set_ylabel('UTM Northing (m)')

    axes[1].imshow(terrain['img'], extent=extent, aspect='equal')
    elev = terrain['elev'].copy()
    elev[elev <= 0] = np.nan
    e_min, e_max = np.nanmin(elev), np.nanmax(elev)
    if np.isfinite(e_min):
        h, w = elev.shape
        x = np.linspace(extent[0], extent[1], w)
        y = np.linspace(extent[3], extent[2], h)
        levels = np.arange(int(e_min) // 5 * 5, int(e_max) + 5, 5)
        axes[1].contour(x, y, elev, levels=levels, colors='blue',
                        linewidths=0.5, alpha=0.6)
        levels_major = np.arange(int(e_min) // 25 * 25, int(e_max) + 25, 25)
        cs2 = axes[1].contour(x, y, elev, levels=levels_major, colors='blue',
                              linewidths=1.5, alpha=0.8)
        axes[1].clabel(cs2, fmt='%d', fontsize=8)
    axes[1].set_title(f'{area_name} — Map + DEM Contours')
    axes[1].set_xlabel('UTM Easting (m)')
    plt.tight_layout()
    path = os.path.join(output_dir, f'{area_name.lower()}_map_overview.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✅ {path}")


def plot_cost_surfaces(terrain, extent, output_dir, area_name, cfg):
    """6-panel: ACR, HCR, combined, elevation, slope, valid mask."""
    valid = create_valid_mask(terrain['img'])
    acr = build_base_cost_from_omap(terrain['cost_omap'], valid)
    hcr_combined, components = build_base_cost_with_hcr(
        terrain['cost_omap'], terrain['elev'], cfg['resolution'], valid,
        cfg.get('hcr_tri_power', 1.0), cfg.get('hcr_pc_power', 0.5),
        cfg.get('hcr_slo_power', 1.0))

    fig, axes = plt.subplots(2, 3, figsize=(24, 16))
    panels = [
        (axes[0, 0], np.clip(acr, 0, 1), 'ACR (Symbol Cost)', 'RdYlGn_r', 0, 1),
        (axes[0, 1], components['hcr_raw'], 'HCR (Morphometric)', 'hot', None, None),
        (axes[0, 2], np.clip(hcr_combined, 0, 2), 'Combined (ACR + HCR)', 'RdYlGn_r', 0, 1.5),
        (axes[1, 0], terrain['elev'], 'Elevation (m)', 'terrain', None, None),
        (axes[1, 1], terrain['slope'], 'Slope (degrees)', 'Reds', 0, 45),
        (axes[1, 2], valid.astype(float), f'Valid Mask ({100 * valid.sum() / valid.size:.0f}%)',
         'RdYlGn', 0, 1),
    ]
    for ax, data, title, cmap, vmin, vmax in panels:
        show = data.copy()
        show[show >= 5] = np.nan
        kwargs = {'extent': extent, 'cmap': cmap, 'aspect': 'equal'}
        if vmin is not None:
            kwargs['vmin'] = vmin
            kwargs['vmax'] = vmax
        im = ax.imshow(show, **kwargs)
        ax.set_title(title)
        plt.colorbar(im, ax=ax, shrink=0.7)

    plt.suptitle(f'{area_name} — Cost Surface Components', fontsize=16, fontweight='bold')
    plt.tight_layout()
    path = os.path.join(output_dir, f'{area_name.lower()}_cost_surfaces.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✅ {path}")


def plot_hillshade(terrain, extent, output_dir, area_name):
    """Hillshade + elevation overlay."""
    ls = LightSource(azdeg=315, altdeg=45)
    elev = terrain['elev'].copy()
    elev[elev <= 0] = np.nan
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    hs = ls.hillshade(np.nan_to_num(elev), vert_exag=3,
                      dx=terrain['resolution'], dy=terrain['resolution'])
    ax.imshow(hs, extent=extent, cmap='gray', aspect='equal')
    im = ax.imshow(elev, extent=extent, cmap='terrain', alpha=0.5, aspect='equal')
    plt.colorbar(im, ax=ax, shrink=0.7, label='Elevation (m)')
    ax.set_title(f'{area_name} — Hillshade + Elevation')
    ax.set_xlabel('UTM Easting (m)')
    ax.set_ylabel('UTM Northing (m)')
    plt.tight_layout()
    path = os.path.join(output_dir, f'{area_name.lower()}_hillshade.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✅ {path}")


def plot_asymmetry(terrain, extent, output_dir, area_name, cfg):
    """Directional cost asymmetry heatmap."""
    sm = np.sqrt(terrain['slope'] ** 2) * np.pi / 180  # convert to radians for Minetti
    phi_up = minetti_factor(np.tan(sm))
    phi_down = minetti_factor(-np.tan(sm))
    asymmetry = np.abs(phi_up - phi_down) / np.maximum(phi_up + phi_down, 1e-9) * 200

    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    im = ax.imshow(np.clip(asymmetry, 0, 100), extent=extent,
                   cmap='RdYlBu_r', vmin=0, vmax=80, aspect='equal')
    plt.colorbar(im, ax=ax, shrink=0.7, label='Directional Asymmetry (%)')
    ax.set_title(f'{area_name} — Cost Asymmetry (Minetti up vs down)')
    ax.set_xlabel('UTM Easting (m)')
    ax.set_ylabel('UTM Northing (m)')
    plt.tight_layout()
    path = os.path.join(output_dir, f'{area_name.lower()}_asymmetry.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✅ {path}")


# ── Route plotting on map (inspired by existing pipeline) ──

ROUTE_COLOR = '#B400FF'
ROUTE_BG = (255, 100, 255)


def plot_route_on_map(terrain, extent, nodes, ctrls, route, legs,
                      output_dir, area_name, dist_name, results_info=None):
    """Plot optimal route overlaid on the orienteering map with A* traced paths."""
    bounds = terrain['bounds']
    resolution = terrain['resolution']

    fig, ax = plt.subplots(1, 1, figsize=(14, 12))
    ax.imshow(terrain['img'], extent=extent, origin='upper', aspect='equal')

    # Draw traced paths
    for path in legs:
        if len(path) < 2:
            continue
        utm_path = [pixel_to_utm(x, y, bounds, resolution) for x, y in path]
        step = max(1, len(utm_path) // 400)
        sub = utm_path[::step]
        if sub[-1] != utm_path[-1]:
            sub.append(utm_path[-1])
        if len(sub) >= 2:
            xs = [p[0] for p in sub]
            ys = [p[1] for p in sub]
            ax.plot(xs, ys, '-', color='white', linewidth=4, alpha=0.6)
            ax.plot(xs, ys, '-', color=ROUTE_COLOR, linewidth=2, alpha=0.9)

    # Draw all controls
    visited_set = set(route)
    for ci, (x, y, code, pv) in enumerate(ctrls):
        ux, uy = pixel_to_utm(x, y, bounds, resolution)
        is_visited = (ci + 1) in visited_set
        if is_visited:
            ax.plot(ux, uy, 'o', color=ROUTE_COLOR, markersize=11,
                    markeredgecolor='white', markeredgewidth=2, zorder=5)
            ax.text(ux + 20, uy, f"{code}({pv})", fontsize=6,
                    color=ROUTE_COLOR, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.15',
                              facecolor='white', alpha=0.7, edgecolor='none'))
        else:
            ax.plot(ux, uy, 'o', color='gray', markersize=5, alpha=0.4)

    # Draw HH
    hux, huy = pixel_to_utm(nodes[0][0], nodes[0][1], bounds, resolution)
    ax.plot(hux, huy, '^', color='red', markersize=16,
            markeredgecolor='white', markeredgewidth=2.5, zorder=6)
    ax.text(hux + 20, huy, 'HH', fontsize=9, color='red', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                      alpha=0.8, edgecolor='red'))

    # Draw visit order
    seq = [0] + list(route) + [0]
    for order, ni in enumerate(seq[1:-1], 1):
        nx, ny = nodes[ni]
        ux, uy = pixel_to_utm(nx, ny, bounds, resolution)
        ax.text(ux - 25, uy - 20, f"#{order}", fontsize=7, color='white',
                fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.15',
                          facecolor=ROUTE_COLOR, alpha=0.8, edgecolor='none'))

    # Title
    pts_total = sum(pv for _, _, _, pv in ctrls)
    pts_got = sum(ctrls[i - 1][3] for i in route if 1 <= i <= len(ctrls))
    title = f"{area_name} — {dist_name}\n{len(route)} controls, {pts_got}/{pts_total} pts"
    if results_info:
        title += f" | {results_info}"
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.set_xlabel('UTM Easting (m)')
    ax.set_ylabel('UTM Northing (m)')
    ax.set_aspect('equal')

    plt.tight_layout()
    fname = f'{area_name.lower()}_{dist_name.lower().replace(" ", "_")}_route.png'
    path_out = os.path.join(output_dir, fname)
    plt.savefig(path_out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✅ {path_out}")


def plot_map_overview(terrain, extent, output_dir, area_name):
    """Map image + DEM contours."""
    fig, axes = plt.subplots(1, 2, figsize=(20, 10))
    axes[0].imshow(terrain['img'], extent=extent, origin='upper', aspect='equal')
    axes[0].set_title(f'{area_name} — Orienteering Map')
    axes[0].set_xlabel('UTM Easting (m)')
    axes[0].set_ylabel('UTM Northing (m)')

    axes[1].imshow(terrain['img'], extent=extent, origin='upper', aspect='equal')
    elev = terrain['elev'].copy()
    elev[elev <= 0] = np.nan
    e_min, e_max = np.nanmin(elev), np.nanmax(elev)
    if np.isfinite(e_min):
        h, w = elev.shape
        x = np.linspace(extent[0], extent[1], w)
        y = np.linspace(extent[3], extent[2], h)
        levels = np.arange(int(e_min) // 5 * 5, int(e_max) + 5, 5)
        axes[1].contour(x, y, elev, levels=levels, colors='blue',
                        linewidths=0.5, alpha=0.6)
        levels_major = np.arange(int(e_min) // 25 * 25, int(e_max) + 25, 25)
        cs2 = axes[1].contour(x, y, elev, levels=levels_major, colors='blue',
                              linewidths=1.5, alpha=0.8)
        axes[1].clabel(cs2, fmt='%d', fontsize=8)
    axes[1].set_title(f'{area_name} — Map + DEM Contours')
    axes[1].set_xlabel('UTM Easting (m)')

    plt.tight_layout()
    path = os.path.join(output_dir, f'{area_name.lower()}_map_overview.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✅ {path}")


def plot_cost_surfaces(terrain, extent, output_dir, area_name, cfg):
    """6-panel: ACR, HCR, combined, elevation, slope, valid mask."""
    valid = create_valid_mask(terrain['img'])
    acr = build_base_cost_from_omap(terrain['cost_omap'], valid)
    hcr_combined, components = build_base_cost_with_hcr(
        terrain['cost_omap'], terrain['elev'], cfg['resolution'], valid,
        cfg.get('hcr_tri_power', 1.0), cfg.get('hcr_pc_power', 0.5),
        cfg.get('hcr_slo_power', 1.0))

    fig, axes = plt.subplots(2, 3, figsize=(24, 16))

    acr_show = acr.copy(); acr_show[acr_show >= 5] = np.nan
    im0 = axes[0, 0].imshow(np.clip(acr_show, 0, 1), extent=extent,
                             cmap='RdYlGn_r', vmin=0, vmax=1,
                             origin='upper', aspect='equal')
    axes[0, 0].set_title('ACR (Symbol Cost)')
    plt.colorbar(im0, ax=axes[0, 0], shrink=0.7)

    hcr_raw = components['hcr_raw']
    im1 = axes[0, 1].imshow(hcr_raw, extent=extent, cmap='hot',
                             origin='upper', aspect='equal')
    axes[0, 1].set_title('HCR (Morphometric)')
    plt.colorbar(im1, ax=axes[0, 1], shrink=0.7)

    comb = hcr_combined.copy(); comb[comb >= 5] = np.nan
    im2 = axes[0, 2].imshow(np.clip(comb, 0, 2), extent=extent,
                             cmap='RdYlGn_r', vmin=0, vmax=1.5,
                             origin='upper', aspect='equal')
    axes[0, 2].set_title('Combined (ACR + HCR)')
    plt.colorbar(im2, ax=axes[0, 2], shrink=0.7)

    im3 = axes[1, 0].imshow(terrain['elev'], extent=extent, cmap='terrain',
                             origin='upper', aspect='equal')
    axes[1, 0].set_title('Elevation (m)')
    plt.colorbar(im3, ax=axes[1, 0], shrink=0.7)

    im4 = axes[1, 1].imshow(terrain['slope'], extent=extent, cmap='Reds',
                             vmin=0, vmax=45, origin='upper', aspect='equal')
    axes[1, 1].set_title('Slope (degrees)')
    plt.colorbar(im4, ax=axes[1, 1], shrink=0.7)

    im5 = axes[1, 2].imshow(valid.astype(float), extent=extent,
                             cmap='RdYlGn', origin='upper', aspect='equal')
    axes[1, 2].set_title(f'Valid Mask ({100 * valid.sum() / valid.size:.0f}%)')
    plt.colorbar(im5, ax=axes[1, 2], shrink=0.7)

    plt.suptitle(f'{area_name} — Cost Surface Components',
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    path = os.path.join(output_dir, f'{area_name.lower()}_cost_surfaces.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✅ {path}")


def plot_hillshade(terrain, extent, output_dir, area_name):
    """Hillshade + elevation overlay."""
    ls = LightSource(azdeg=315, altdeg=45)
    elev = terrain['elev'].copy()
    elev[elev <= 0] = np.nan

    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    hs = ls.hillshade(np.nan_to_num(elev), vert_exag=3,
                      dx=terrain['resolution'], dy=terrain['resolution'])
    ax.imshow(hs, extent=extent, cmap='gray', origin='upper', aspect='equal')
    im = ax.imshow(elev, extent=extent, cmap='terrain', alpha=0.5,
                   origin='upper', aspect='equal')
    plt.colorbar(im, ax=ax, shrink=0.7, label='Elevation (m)')
    ax.set_title(f'{area_name} — Hillshade + Elevation')
    ax.set_xlabel('UTM Easting (m)')
    ax.set_ylabel('UTM Northing (m)')

    plt.tight_layout()
    path = os.path.join(output_dir, f'{area_name.lower()}_hillshade.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✅ {path}")


def plot_asymmetry(terrain, extent, output_dir, area_name, cfg):
    """Directional cost asymmetry heatmap."""
    sx, sy = compute_slope_dir(terrain['elev'], cfg['resolution'])
    sm = np.sqrt(sx ** 2 + sy ** 2)
    phi_up = minetti_factor(sm)
    phi_down = minetti_factor(-sm)
    asymmetry = np.abs(phi_up - phi_down) / np.maximum(phi_up + phi_down, 1e-9) * 200

    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    im = ax.imshow(np.clip(asymmetry, 0, 100), extent=extent,
                   cmap='RdYlBu_r', vmin=0, vmax=80,
                   origin='upper', aspect='equal')
    plt.colorbar(im, ax=ax, shrink=0.7, label='Directional Asymmetry (%)')
    ax.set_title(f'{area_name} — Cost Asymmetry (Minetti up vs down)')
    ax.set_xlabel('UTM Easting (m)')
    ax.set_ylabel('UTM Northing (m)')

    plt.tight_layout()
    path = os.path.join(output_dir, f'{area_name.lower()}_asymmetry.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✅ {path}")


# ── Main ──

def main():
    parser = argparse.ArgumentParser(description='Visualize orienteering results')
    parser.add_argument('area', choices=['torremocha', 'la_muela'],
                        help='Study area name')
    parser.add_argument('--output-dir', default='figures',
                        help='Output directory for figures')
    parser.add_argument('--terrain-cache', default=None,
                        help='Path to terrain .npz cache')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    if args.area == 'torremocha':
        from config.torremocha import (TARGET_RESOLUTION, TERRAIN_CACHE,
                                        HCR_TRI_POWER, HCR_PC_POWER, HCR_SLO_POWER)
        area_name = 'Torremocha'
    else:
        from config.la_muela import (TARGET_RESOLUTION, TERRAIN_CACHE,
                                      HCR_TRI_POWER, HCR_PC_POWER, HCR_SLO_POWER)
        area_name = 'La Muela'

    cfg = {
        'resolution': TARGET_RESOLUTION,
        'hcr_tri_power': HCR_TRI_POWER,
        'hcr_pc_power': HCR_PC_POWER,
        'hcr_slo_power': HCR_SLO_POWER,
    }

    cache_path = args.terrain_cache or TERRAIN_CACHE
    print(f"Loading terrain: {cache_path}")
    terrain = load_terrain(cache_path)
    extent = make_extent(terrain['bounds'])

    print(f"\nGenerating figures for {area_name}...")
    plot_map_overview(terrain, extent, args.output_dir, area_name)
    plot_cost_surfaces(terrain, extent, args.output_dir, area_name, cfg)
    plot_hillshade(terrain, extent, args.output_dir, area_name)
    plot_asymmetry(terrain, extent, args.output_dir, area_name, cfg)

    print(f"\n✅ All figures saved to {args.output_dir}/")


if __name__ == '__main__':
    main()
