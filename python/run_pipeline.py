"""
Main pipeline: load terrain data, build cost surface, optimize route, export JSONs.

Usage:
    python run_pipeline.py torremocha
    python run_pipeline.py la_muela
"""
import sys
import os
import json
import time
import numpy as np
import random

# Add parent to path for imports
sys.path.insert(0, os.path.dirname(__file__))

from core.cost_functions import (
    WT, minetti_factor, build_base_cost_from_omap, build_display_cost,
    cost_to_veg_proxy, cost_str, effective_budget
)
from core.terrain_analysis import build_base_cost_with_hcr
from core.pathfinding import (
    compute_slope_dir, compute_cost_matrix, rcost, rcost_fatigue,
    rcost_fatigue_detail, rpts
)
from core.route_optimizer import optimize_route
from core.control_placement import (
    create_valid_mask, find_valid_hh,
    place_controls_standard, place_controls_clustered,
    place_controls_ring, place_controls_sparse_far,
    place_controls_mixed_density, place_controls_path_biased,
    place_controls_elev_biased
)


def load_config(area_name):
    """Load location-specific configuration."""
    if area_name == 'torremocha':
        from config import torremocha as cfg_mod
    elif area_name == 'la_muela':
        from config import la_muela as cfg_mod
    elif area_name == 'pinofranqueado':
        from config import pinofranqueado as cfg_mod
    else:
        raise ValueError(f"Unknown area: {area_name}")

    return {
        'map_crs': cfg_mod.MAP_CRS, 'resolution': cfg_mod.TARGET_RESOLUTION,
        'seed': cfg_mod.SEED, 'num_controls': cfg_mod.NUM_CONTROLS,
        'base_speed': cfg_mod.BASE_SPEED, 'reference_weight': cfg_mod.REFERENCE_WEIGHT,
        'race_hours': cfg_mod.RACE_DURATION_HOURS, 'fatigue_rate': cfg_mod.FATIGUE_RATE,
        'route_budget': cfg_mod.ROUTE_BUDGET, 'cost_to_seconds': cfg_mod.COST_TO_SECONDS,
        'downsample': cfg_mod.DOWNSAMPLE,
        'hcr_tri_power': cfg_mod.HCR_TRI_POWER, 'hcr_pc_power': cfg_mod.HCR_PC_POWER,
        'hcr_slo_power': cfg_mod.HCR_SLO_POWER,
        'hcr_norm_low': cfg_mod.HCR_NORM_LOW, 'hcr_norm_high': cfg_mod.HCR_NORM_HIGH,
        'terrain_cache': cfg_mod.TERRAIN_CACHE, 'distributions': cfg_mod.DISTRIBUTIONS,
        # Preprocessing paths
        'tif_path': cfg_mod.TIF_PATH, 'mdt_path': cfg_mod.MDT_PATH,
        'omap_path': cfg_mod.OMAP_PATH,
        'world_A': cfg_mod.WORLD_A, 'world_B': cfg_mod.WORLD_B,
        'world_C': cfg_mod.WORLD_C, 'world_D': cfg_mod.WORLD_D,
        'world_E': cfg_mod.WORLD_E, 'world_F': cfg_mod.WORLD_F,
    }


def load_terrain(cache_path):
    """Load preprocessed terrain data from .npz cache."""
    print(f"Loading terrain: {cache_path}")
    with np.load(cache_path, allow_pickle=True) as terrain:
        data = {
            'elev': terrain['elevation'],
            'slope': terrain['slope'],
            'cost_omap': terrain['cost_omap'],
            'img': terrain['img'],
            'bounds': tuple(terrain['bounds']),
            'resolution': float(terrain['resolution']),
        }
        if 'path_grid' in terrain:
            data['path_grid'] = terrain['path_grid']
        if 'wall_grid' in terrain:
            data['wall_grid'] = terrain['wall_grid']
        if 'water_grid' in terrain:
            data['water_grid'] = terrain['water_grid']
    h, w = data['elev'].shape
    print(f"  Grid: {w}x{h} @ {data['resolution']:.1f}m/cell")
    print(f"  Elevation: {np.nanmin(data['elev']):.0f}–{np.nanmax(data['elev']):.0f} m")
    return data


PLACEMENT_FUNCS = {
    'standard': place_controls_standard,
    'clustered': place_controls_clustered,
    'ring': place_controls_ring,
    'path_biased': place_controls_path_biased,
    'elev_biased': place_controls_elev_biased,
    'sparse_far': place_controls_sparse_far,
    'mixed_density': place_controls_mixed_density,
}


def run_main_route(cfg, terrain):
    """Run the main rogaine route optimization."""
    print("=" * 65)
    print(f"  ROGAINE — {cfg['map_crs']} (ACR + HCR + Minetti + Fatigue)")
    print("=" * 65)

    # Valid mask
    print("\n[1/7] Valid map area...")
    valid = create_valid_mask(terrain['img'])

    # HCR + base cost
    print("\n[2/7] Computing HCR...")
    t0 = time.time()
    base_cost_acr = build_base_cost_from_omap(terrain['cost_omap'], valid)
    base_cost, hcr_components = build_base_cost_with_hcr(
        terrain['cost_omap'], terrain['elev'], cfg['resolution'], valid,
        cfg['hcr_tri_power'], cfg['hcr_pc_power'], cfg['hcr_slo_power'],
        cfg['hcr_norm_low'], cfg['hcr_norm_high'])
    t_hcr = time.time() - t0
    print(f"  HCR computation: {t_hcr:.1f}s")

    flat = base_cost[base_cost < 5]
    flat_acr = base_cost_acr[base_cost_acr < 5]
    print(f"  Combined cost: {flat.min():.3f} / {np.median(flat):.3f} / {flat.max():.3f}")
    print(f"  ACR-only cost: {flat_acr.min():.3f} / {np.median(flat_acr):.3f} / {flat_acr.max():.3f}")
    diff_pct = 100 * (np.median(flat) - np.median(flat_acr)) / np.median(flat_acr)
    print(f"  HCR median increase: {diff_pct:+.1f}%")

    # Slopes
    print("\n[3/7] Directional slopes...")
    sx, sy = compute_slope_dir(terrain['elev'], cfg['resolution'])
    sm = np.sqrt(sx ** 2 + sy ** 2)
    print(f"  Max slope: {np.degrees(np.arctan(sm.max())):.1f}°")

    # Controls
    veg_grid = cost_to_veg_proxy(base_cost)
    print("\n[4/7] Hash house & controls...")
    hh_x, hh_y = find_valid_hh(valid, veg_grid, terrain['slope'], cfg['resolution'])
    hh = (hh_x, hh_y)
    print(f"  HH at ({hh_x}, {hh_y}), elev={terrain['elev'][hh_y, hh_x]:.0f}m")

    ctrls = place_controls_standard(veg_grid, terrain['elev'], terrain['slope'],
                                     valid, cfg['num_controls'], hh, cfg['resolution'],
                                     seed=cfg['seed'])
    nodes = [hh]
    pa = [0.0]
    for cx, cy, code, pts in ctrls:
        nodes.append((cx, cy)); pa.append(float(pts))
    pa = np.array(pa)

    # Cost matrix
    print(f"\n[5/7] Cost matrix ({len(nodes)} nodes, DS={cfg['downsample']})...")
    t0 = time.time()
    cm = compute_cost_matrix(base_cost, sx, sy, nodes, cfg['downsample'])
    tm = time.time() - t0
    print(f"  Done in {tm:.1f}s")

    asym = np.abs(cm - cm.T)
    vld = np.isfinite(cm) & (cm > 0)
    apct = 100 * asym[vld].mean() / cm[vld].mean() if vld.any() else 0
    print(f"  Asymmetry: {apct:.1f}%")

    fm = np.isfinite(cm)
    ic = (~fm).sum() - len(nodes)
    if ic > 0:
        print(f"  ⚠️ {ic} unreachable pairs")
        cm[~fm] = cm[fm].max() * 3

    bud_raw = cfg['route_budget']
    bud_eff = effective_budget(bud_raw, cfg['fatigue_rate'])

    # Optimize
    print(f"\n[6/7] Optimizing route...")
    t0 = time.time()
    route = optimize_route(cm, pa, bud_eff, bud_raw, cfg['fatigue_rate'])
    to = time.time() - t0

    # Results
    tp = rpts(route, pa); cb = rcost(cm, route)
    cf = rcost_fatigue(cm, route, bud_raw, cfg['fatigue_rate'])
    mp = sum(pa)
    details = rcost_fatigue_detail(cm, route, bud_raw, cfg['fatigue_rate'])

    cs = lambda v: cost_str(v, cfg['resolution'], cfg['base_speed'], cfg['reference_weight'])
    print(f"\n{'=' * 65}")
    print(f"  RESULTS (with HCR)")
    print(f"{'=' * 65}")
    print(f"  Controls:  {len(route)}/{len(ctrls)}")
    print(f"  Points:    {int(tp)}/{int(mp)} ({100 * tp / mp:.1f}%)")
    print(f"  Base:      {cs(cb)}")
    print(f"  Fatigue:   {cs(cf)} (+{100 * (cf - cb) / max(cb, 1):.1f}%)")
    print(f"  Asymmetry: {apct:.1f}%")

    return {
        'valid': valid, 'base_cost': base_cost, 'sx': sx, 'sy': sy,
        'veg_grid': veg_grid, 'hh': hh, 'ctrls': ctrls,
        'route': route, 'cm': cm, 'pa': pa, 'nodes': nodes,
    }


def export_jsons(cfg, terrain, pipeline_result):
    """Export JSON files for the C++ B&C solver."""
    print("\n" + "=" * 65)
    print("  EXPORTING JSON FILES FOR C++ SOLVER")
    print("=" * 65)

    valid = pipeline_result['valid']
    base_cost = pipeline_result['base_cost']
    sx = pipeline_result['sx']
    sy = pipeline_result['sy']
    veg_grid = pipeline_result['veg_grid']
    hh = pipeline_result['hh']

    bud_raw = cfg['route_budget']
    bud_eff = effective_budget(bud_raw, cfg['fatigue_rate'])

    for dist_name, dist_cfg in cfg['distributions'].items():
        print(f"\n── {dist_name} ──")
        func_name = dist_cfg['func']
        seed = dist_cfg['seed']
        num = dist_cfg['num']

        if func_name not in PLACEMENT_FUNCS:
            print(f"  ⚠️ Unknown placement function: {func_name}, skipping")
            continue

        func = PLACEMENT_FUNCS[func_name]
        ctrls = func(veg_grid, terrain['elev'], terrain['slope'],
                     valid, num, hh, cfg['resolution'], seed=seed)
        print(f"  Placed {len(ctrls)} controls")

        nodes = [hh]; pa = [0.0]
        for cx, cy, code, pts in ctrls:
            nodes.append((cx, cy)); pa.append(float(pts))
        pa = np.array(pa)

        cm = compute_cost_matrix(base_cost, sx, sy, nodes, cfg['downsample'])
        fm = np.isfinite(cm)
        if (~fm).sum() - len(nodes) > 0:
            cm[~fm] = cm[fm].max() * 3

        filename = f"op_input_{dist_name}.json"
        op_data = {
            "cm": cm.tolist(),
            "pts": pa.tolist(),
            "bud_eff": float(bud_eff),
            "bud_raw": float(bud_raw),
            "fatigue_rate": float(cfg['fatigue_rate']),
        }
        with open(filename, 'w') as f:
            json.dump(op_data, f, indent=2)
        print(f"  ✅ Exported {filename} — {len(pa)} nodes, total pts={sum(pa):.0f}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python run_pipeline.py <torremocha|la_muela> [--preprocess]")
        sys.exit(1)

    area = sys.argv[1].lower()
    preprocess = '--preprocess' in sys.argv
    cfg = load_config(area)

    random.seed(cfg['seed'])
    np.random.seed(cfg['seed'])

    if preprocess or not os.path.exists(cfg['terrain_cache']):
        print("Running preprocessing (raw files → terrain cache)...")
        from core.preprocessing import run_preprocessing
        run_preprocessing(cfg)

    terrain = load_terrain(cfg['terrain_cache'])
    result = run_main_route(cfg, terrain)
    export_jsons(cfg, terrain, result)

    print("\n✅ Pipeline complete.")


if __name__ == "__main__":
    main()
