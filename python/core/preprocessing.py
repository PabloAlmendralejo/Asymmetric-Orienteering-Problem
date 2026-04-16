"""
Preprocessing: load map image, DEM, OMAP file → rasterize → build terrain cache (.npz).
This is the first stage of the pipeline, converting raw geospatial data into
the grid-based cost surface used by the pathfinding and optimization stages.
"""
import numpy as np
import math
import xml.etree.ElementTree as ET
from collections import Counter

import rasterio
from rasterio.windows import from_bounds
from scipy.ndimage import map_coordinates
from scipy.interpolate import griddata
from pyproj import Transformer
from PIL import Image, ImageDraw

from .omap_parser import (
    ISOM_CATEGORIES, CAT_COST, SYM2CAT, NAME_TO_CAT,
    get_category, parse_omap_coords_text
)
from .rasterization import (
    rast_line_omap, rast_line_cost_omap, rast_line_barrier_omap,
    rast_poly_omap, rast_poly_barrier_omap, rast_poly_bool_omap
)


def load_map_image(tif_path, world_A, world_B, world_C, world_D, world_E, world_F):
    """Load map TIF and compute UTM bounds from world file parameters."""
    img_pil = Image.open(tif_path)
    W, H = img_pil.size
    img_array = np.array(img_pil.convert('RGB'))

    def pixel_to_utm(col, row):
        x = world_A * col + world_B * row + world_C
        y = world_D * col + world_E * row + world_F
        return x, y

    def utm_to_pixel(x, y):
        det = world_A * world_E - world_B * world_D
        col = (world_E * (x - world_C) - world_B * (y - world_F)) / det
        row = (-world_D * (x - world_C) + world_A * (y - world_F)) / det
        return col, row

    corners_utm = [pixel_to_utm(c, r) for c, r in [(0, 0), (W, 0), (W, H), (0, H)]]
    all_x = [c[0] for c in corners_utm]
    all_y = [c[1] for c in corners_utm]
    bounds = (min(all_x), min(all_y), max(all_x), max(all_y))

    print(f"Map bounds (UTM): ({bounds[0]:.1f}, {bounds[1]:.1f}) → ({bounds[2]:.1f}, {bounds[3]:.1f})")
    return img_array, W, H, bounds, pixel_to_utm, utm_to_pixel


def load_dem(mdt_path, map_crs, bounds_utm):
    """Load DEM, crop to map extent, and resample to UTM grid.
    Handles DEMs in either EPSG:4326 or the same CRS as the map."""
    x_min, y_min, x_max, y_max = bounds_utm
    transformer = Transformer.from_crs(map_crs, "EPSG:4326", always_xy=True)

    corners_ll = [
        transformer.transform(x_min, y_min), transformer.transform(x_max, y_min),
        transformer.transform(x_max, y_max), transformer.transform(x_min, y_max),
    ]
    lons = [c[0] for c in corners_ll]
    lats = [c[1] for c in corners_ll]
    lon_min, lon_max = min(lons), max(lons)
    lat_min, lat_max = min(lats), max(lats)

    print(f"Map in lat/lon: ({lon_min:.6f}, {lat_min:.6f}) → ({lon_max:.6f}, {lat_max:.6f})")

    with rasterio.open(mdt_path) as src:
        dem_crs = src.crs
        # If DEM is in the same CRS as the map, crop in UTM coordinates
        if dem_crs is not None and dem_crs.to_epsg() == rasterio.crs.CRS.from_string(map_crs).to_epsg():
            pad = 100  # metres
            window = from_bounds(x_min - pad, y_min - pad,
                                 x_max + pad, y_max + pad, src.transform)
            print(f"DEM is in {map_crs} — cropping in UTM coordinates")
            dem_is_utm = True
        else:
            pad = 0.002  # degrees
            window = from_bounds(lon_min - pad, lat_min - pad,
                                 lon_max + pad, lat_max + pad, src.transform)
            print(f"DEM is in {dem_crs} — cropping in lat/lon coordinates")
            dem_is_utm = False

        elev_crop = src.read(1, window=window).astype(float)
        crop_transform = src.window_transform(window)
        crop_bounds = rasterio.windows.bounds(window, src.transform)
        nodata = src.nodata

    if nodata is not None:
        elev_crop[elev_crop == nodata] = np.nan

    print(f"MDT02 crop: {elev_crop.shape}")
    print(f"Elevation: {np.nanmin(elev_crop):.1f} – {np.nanmax(elev_crop):.1f} m")

    return elev_crop, crop_transform, crop_bounds, transformer, dem_is_utm


def build_utm_grid(bounds_utm, resolution, elev_crop, crop_transform, crop_bounds,
                   transformer_to_ll, img_array, utm_to_pixel_func, dem_is_utm=False):
    """Build UTM grid, resample DEM and image onto it."""
    x_min, y_min, x_max, y_max = bounds_utm

    grid_x = np.arange(x_min, x_max, resolution)
    grid_y = np.arange(y_max, y_min, -resolution)
    nx, ny = len(grid_x), len(grid_y)
    gx, gy = np.meshgrid(grid_x, grid_y)
    print(f"\nUTM grid: {nx} x {ny}")

    # Resample DEM
    crop_left, crop_bottom, crop_right, crop_top = crop_bounds
    if dem_is_utm:
        # DEM is already in UTM — map UTM coordinates directly to pixel indices
        pixel_size_x = crop_transform[0]
        pixel_size_y = abs(crop_transform[4])
        pixel_cols_f = (gx - crop_left) / pixel_size_x
        pixel_rows_f = (crop_top - gy) / pixel_size_y
    else:
        # DEM is in lat/lon — convert UTM grid to lat/lon first
        grid_lon, grid_lat = transformer_to_ll.transform(gx, gy)
        pixel_size_lon = crop_transform[0]
        pixel_size_lat = abs(crop_transform[4])
        pixel_cols_f = (grid_lon - crop_left) / pixel_size_lon
        pixel_rows_f = (crop_top - grid_lat) / pixel_size_lat

    elev_utm = map_coordinates(elev_crop, [pixel_rows_f, pixel_cols_f],
                                order=3, mode='nearest', cval=np.nan)
    print(f"Elevation grid: {elev_utm.shape}")
    print(f"NaN: {np.isnan(elev_utm).sum()}")
    print(f"Elevation: {np.nanmin(elev_utm):.1f} – {np.nanmax(elev_utm):.1f} m")

    if np.isnan(elev_utm).any():
        valid = ~np.isnan(elev_utm)
        coords_valid = np.column_stack([gx[valid], gy[valid]])
        coords_nan = np.column_stack([gx[~valid], gy[~valid]])
        elev_utm[~valid] = griddata(coords_valid, elev_utm[valid],
                                     coords_nan, method='nearest')

    # Remap image
    map_cols, map_rows = utm_to_pixel_func(gx, gy)
    img_utm = np.zeros((ny, nx, 3), dtype=np.uint8)
    for ch in range(3):
        img_utm[:, :, ch] = map_coordinates(
            img_array[:, :, ch].astype(float),
            [map_rows, map_cols], order=1, mode='constant', cval=255
        ).clip(0, 255).astype(np.uint8)
    print(f"Remapped image: {img_utm.shape}")

    return grid_x, grid_y, nx, ny, elev_utm, img_utm, gx, gy


def parse_omap_georef(omap_path, map_crs, world_D, world_A):
    """Parse OMAP XML file and extract georeferencing + symbol data."""
    NS = 'http://openorienteering.org/apps/mapper/xml/v2'
    NSP = f'{{{NS}}}'

    with open(omap_path, 'rb') as f:
        raw = f.read()
    text = raw.decode('utf-8', errors='replace').replace('\ufffd', '?')
    root = ET.fromstring(text)

    # Read georeferencing
    proj_x = proj_y = scale = None
    paper_ref_x = paper_ref_y = 0.0
    grivation = 0.0
    crs_spec = ""

    for georef in root.iter(f'{NSP}georeferencing'):
        scale = float(georef.get('scale', 15000))
        grivation = float(georef.get('grivation', '0'))
        rp = georef.find(f'{NSP}ref_point')
        if rp is not None:
            paper_ref_x = float(rp.get('x', 0))
            paper_ref_y = float(rp.get('y', 0))
        pcrs = georef.find(f'{NSP}projected_crs')
        if pcrs is not None:
            rp2 = pcrs.find(f'{NSP}ref_point')
            if rp2 is not None:
                proj_x = float(rp2.get('x', 0))
                proj_y = float(rp2.get('y', 0))
            spec_elem = pcrs.find(f'{NSP}spec')
            if spec_elem is not None and spec_elem.text:
                crs_spec = spec_elem.text.strip()

    if proj_x is None:
        raise ValueError("Could not find projection center in OMAP XML!")

    print(f"  Projection center: ({proj_x:.3f}, {proj_y:.3f})")
    print(f"  Scale: 1:{scale:.0f}")
    print(f"  Grivation: {grivation}°")

    # Detect OMAP CRS
    omap_crs = None
    spec_upper = crs_spec.upper()
    if 'WGS84' in spec_upper or 'wgs84' in crs_spec:
        # Try to detect zone
        import re
        zone_match = re.search(r'zone=(\d+)', crs_spec, re.IGNORECASE)
        if zone_match:
            zone = int(zone_match.group(1))
            omap_crs = f"EPSG:326{zone:02d}"
        else:
            omap_crs = "EPSG:32629"
    elif '25829' in crs_spec:
        omap_crs = "EPSG:25829"
    elif '25830' in crs_spec:
        omap_crs = "EPSG:25830"
    elif '23029' in crs_spec:
        omap_crs = "EPSG:23029"
    elif '23030' in crs_spec:
        omap_crs = "EPSG:23030"
    if omap_crs is None:
        omap_crs = "EPSG:32629"
        print(f"  ⚠️ Guessing OMAP CRS: {omap_crs}")
    else:
        print(f"  OMAP CRS: {omap_crs}")

    # CRS transformer
    omap_to_grid = None
    if omap_crs != map_crs:
        omap_to_grid = Transformer.from_crs(omap_crs, map_crs, always_xy=True)

    # World file rotation
    wf_rotation_deg = math.degrees(math.atan2(world_D, world_A))

    # Parse symbols
    omap_symbols = {}
    unparsed_codes = []
    for elem in root.iter():
        if elem.tag == f'{NSP}symbol':
            sid = elem.get('id', '')
            code = elem.get('code', '')
            name = elem.get('name', '')
            isom = None
            try:
                isom = float(code)
            except (ValueError, TypeError):
                clean = code.strip()
                for prefix in ['ISOM', 'ISSprOM', 'ISSOM', '#']:
                    clean = clean.replace(prefix, '').strip()
                filtered = ''
                dot_seen = False
                for ch in clean:
                    if ch.isdigit():
                        filtered += ch
                    elif ch == '.' and not dot_seen:
                        filtered += ch
                        dot_seen = True
                if filtered:
                    try:
                        isom = float(filtered)
                    except (ValueError, TypeError):
                        pass
            if isom is None and code.strip():
                unparsed_codes.append((sid, code, name))
            omap_symbols[sid] = {'isom': isom, 'name': name, 'code': code}

    if unparsed_codes:
        print(f"\n  ⚠️ Objects with unparseable symbol code:")
        type_counts = Counter()
        for sid, code, name in unparsed_codes:
            type_counts[code] += 1

    print(f"\n  Symbols: {len(omap_symbols)}")

    return {
        'root': root, 'NSP': NSP,
        'proj_x': proj_x, 'proj_y': proj_y, 'scale': scale,
        'paper_ref_x': paper_ref_x, 'paper_ref_y': paper_ref_y,
        'grivation': grivation, 'omap_crs': omap_crs,
        'omap_to_grid': omap_to_grid, 'wf_rotation_deg': wf_rotation_deg,
        'omap_symbols': omap_symbols,
    }


def rasterize_omap(omap_data, bounds_utm, resolution, nx, ny):
    """Rasterize OMAP objects onto the UTM grid to produce cost_omap and auxiliary grids."""
    root = omap_data['root']
    NSP = omap_data['NSP']
    omap_symbols = omap_data['omap_symbols']
    proj_x = omap_data['proj_x']
    proj_y = omap_data['proj_y']
    scale = omap_data['scale']
    paper_ref_x = omap_data['paper_ref_x']
    paper_ref_y = omap_data['paper_ref_y']
    wf_rotation_deg = omap_data['wf_rotation_deg']
    omap_to_grid_tf = omap_data['omap_to_grid']

    x_min, y_min, x_max, y_max = bounds_utm

    # Build coordinate transforms
    rot_rad = math.radians(wf_rotation_deg)
    cos_r = math.cos(rot_rad)
    sin_r = math.sin(rot_rad)

    def paper_to_utm(coords_mm):
        result = []
        for px, py in coords_mm:
            dx_mm = px - paper_ref_x
            dy_mm = py - paper_ref_y
            dx_m = dx_mm * scale / 1000.0
            dy_m = dy_mm * scale / 1000.0
            dx_rot = dx_m * cos_r + dy_m * sin_r
            dy_rot = -dx_m * sin_r + dy_m * cos_r
            utm_x = proj_x + dx_rot
            utm_y = proj_y - dy_rot
            result.append((utm_x, utm_y))
        if omap_to_grid_tf is not None:
            result = [omap_to_grid_tf.transform(x, y) for x, y in result]
        return result

    def utm_to_grid(utm_coords):
        return [((ux - x_min) / resolution,
                 (y_max - uy) / resolution)
                for ux, uy in utm_coords]

    # Initialize grids
    cost_omap = np.full((ny, nx), -1.0, dtype=np.float64)
    path_grid = np.zeros((ny, nx), dtype=np.uint8)
    water_grid = np.zeros((ny, nx), dtype=np.uint8)
    wall_grid = np.zeros((ny, nx), dtype=np.uint8)

    print("  Rasterizing...")

    # Process objects
    n_paths = n_areas = n_points = n_oob = 0
    n_total_obj = 0
    unknown_codes = Counter()

    all_objects = list(root.iter(f'{NSP}object'))
    print(f"  Total objects: {len(all_objects)}")

    for obj_idx, obj in enumerate(all_objects):
        if obj_idx % 500 == 0 and obj_idx > 0:
            print(f"  ... processed {obj_idx}/{len(all_objects)} objects")
        sym_id = obj.get('symbol', '')
        obj_type = int(obj.get('type', -1))
        sym_info = omap_symbols.get(sym_id, {})
        cat = get_category(sym_info)

        if cat is None or cat in ('contour', 'ignore'):
            continue

        cost_val = CAT_COST.get(cat)
        if cost_val is None:
            continue

        # Parse coordinates
        coords_elem = obj.find(f'{NSP}coords')
        if coords_elem is None or not coords_elem.text:
            continue
        paper_coords, flags = parse_omap_coords_text(coords_elem.text)
        if not paper_coords:
            continue

        utm_coords = paper_to_utm(paper_coords)
        grid_coords = utm_to_grid(utm_coords)

        is_barrier = cat in ('cliff', 'wall', 'fence', 'building',
                             'lake', 'water', 'river')
        is_path = cat in ('paved_road', 'road', 'track', 'footpath',
                          'small_path', 'faint_path')

        if obj_type == 1:  # Line
            buf = 2 if is_path else 1
            if is_barrier:
                rast_line_barrier_omap(cost_omap, grid_coords, cost_val, buf)
                if cat == 'wall':
                    rast_line_omap(wall_grid, grid_coords, 1, buf)
                elif cat in ('lake', 'water', 'river'):
                    rast_line_omap(water_grid, grid_coords, 1, buf)
            else:
                rast_line_cost_omap(cost_omap, grid_coords, cost_val, buf)
            if is_path:
                path_val = 2 if cat in ('paved_road', 'road', 'track') else 1
                rast_line_omap(path_grid, grid_coords, path_val, buf)
            n_paths += 1 if is_path else 0
            n_areas += 1

        elif obj_type == 0:  # Point
            gx, gy = grid_coords[0] if grid_coords else (0, 0)
            gi, gj = int(round(gy)), int(round(gx))
            if 0 <= gi < ny and 0 <= gj < nx:
                if is_barrier:
                    cost_omap[gi, gj] = max(cost_omap[gi, gj], cost_val)
                elif cost_omap[gi, gj] < 0 or cost_val < cost_omap[gi, gj]:
                    cost_omap[gi, gj] = cost_val
            n_points += 1

    # Fill unset cells with default forest cost
    cost_omap[cost_omap < 0] = CAT_COST['forest_good']

    non_default = (cost_omap != CAT_COST['forest_good']).sum()
    total = ny * nx
    print(f"\n  Non-default: {non_default:,} / {total:,} ({100 * non_default / total:.1f}%)")
    print(f"  Cost range: {cost_omap.min():.3f} – {cost_omap.max():.3f}")
    print(f"  Paths: {(path_grid > 0).sum():,}")
    print(f"  Water: {water_grid.sum():,}")
    print(f"  Wall pixels: {(wall_grid > 0).sum():,}")

    return cost_omap, path_grid, water_grid, wall_grid


def classify_vegetation(img_rgb):
    """Classify vegetation from RGB image (green channel proxy)."""
    r, g, b = img_rgb[:, :, 0].astype(float), img_rgb[:, :, 1].astype(float), img_rgb[:, :, 2].astype(float)
    green_ratio = g / np.maximum(r + g + b, 1.0)
    veg = np.clip((green_ratio - 0.25) / 0.25, 0.0, 1.0)
    return veg.astype(np.float32)


def save_terrain_cache(cache_path, nx, ny, resolution, bounds_utm,
                       elev_utm, img_utm, cost_omap, path_grid,
                       water_grid, wall_grid, veg_grid, slope_grid):
    """Save preprocessed terrain data to .npz cache."""
    np.savez_compressed(cache_path,
        vegetation=veg_grid, elevation=elev_utm, slope=slope_grid,
        cost_omap=cost_omap, path_grid=path_grid, wall_grid=wall_grid,
        water_grid=water_grid, img=img_utm,
        bounds=np.array(bounds_utm), resolution=resolution)

    print(f"\n✅ Saved {cache_path}")
    print(f"   Grid:       {nx} x {ny}")
    print(f"   Resolution: {resolution}m")
    print(f"   Elevation:  {np.nanmin(elev_utm):.0f} – {np.nanmax(elev_utm):.0f} m")
    print(f"   Slope:      {slope_grid.min():.1f}° – {slope_grid.max():.1f}°")
    print(f"   Cost:       {cost_omap.min():.3f} – {cost_omap.max():.3f}")


def run_preprocessing(cfg):
    """Full preprocessing pipeline: raw files → terrain cache .npz."""
    from .cost_functions import cost_to_veg_proxy

    print("=" * 65)
    print(f"  PREPROCESSING — {cfg['map_crs']}")
    print("=" * 65)

    # Load map image
    img_array, W, H, bounds, pixel_to_utm, utm_to_pixel = load_map_image(
        cfg['tif_path'], cfg['world_A'], cfg['world_B'], cfg['world_C'],
        cfg['world_D'], cfg['world_E'], cfg['world_F'])

    # Load DEM
    elev_crop, crop_tf, crop_bounds, tf_to_ll, dem_is_utm = load_dem(
        cfg['mdt_path'], cfg['map_crs'], bounds)

    # Build UTM grid
    grid_x, grid_y, nx, ny, elev_utm, img_utm, gx, gy = build_utm_grid(
        bounds, cfg['resolution'], elev_crop, crop_tf, crop_bounds,
        tf_to_ll, img_array, utm_to_pixel, dem_is_utm=dem_is_utm)

    # Parse and rasterize OMAP
    print(f"\n[3b] OMAP rasterization to UTM grid")
    omap_data = parse_omap_georef(
        cfg['omap_path'], cfg['map_crs'], cfg['world_D'], cfg['world_A'])
    cost_omap, path_grid, water_grid, wall_grid = rasterize_omap(
        omap_data, bounds, cfg['resolution'], nx, ny)

    # Vegetation and slope
    veg_grid = classify_vegetation(img_utm)
    dy_g, dx_g = np.gradient(elev_utm, cfg['resolution'])
    slope_grid = np.degrees(np.arctan(np.sqrt(dx_g ** 2 + dy_g ** 2)))

    # Clamp extreme slopes
    slope_max = 60.0
    n_extreme = (slope_grid > slope_max).sum()
    if n_extreme > 0:
        slope_grid = np.clip(slope_grid, 0, slope_max)
        print(f"\n  Slope clamped: {n_extreme:,} pixels > {slope_max}°")

    # DEM coverage mask
    dem_coverage = np.isfinite(elev_utm) & (elev_utm > 0)

    # Save cache
    save_terrain_cache(cfg['terrain_cache'], nx, ny, cfg['resolution'], bounds,
                       elev_utm, img_utm, cost_omap, path_grid,
                       water_grid, wall_grid, veg_grid, slope_grid)

    return {
        'elev': elev_utm, 'img': img_utm, 'cost_omap': cost_omap,
        'path_grid': path_grid, 'water_grid': water_grid,
        'wall_grid': wall_grid, 'veg_grid': veg_grid,
        'slope': slope_grid, 'bounds': bounds,
        'resolution': cfg['resolution'], 'nx': nx, 'ny': ny,
    }
