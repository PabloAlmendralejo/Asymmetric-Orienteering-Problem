"""
Merge and reproject Pinofranqueado elevation tiles.
Input: two DEM tiles in EPSG:4326 (WGS84)
Output: merged DEM in EPSG:25829 (UTM 29N) at 2m resolution
"""
import rasterio
from rasterio.merge import merge
from rasterio.warp import calculate_default_transform, reproject, Resampling
import numpy as np
import sys
import os

# Paths
east_tile = r"C:\Users\borrepa\Downloads\PinoFranqueado_este.tif"
west_tile = r"C:\Users\borrepa\Downloads\Pinofranquead_oeste.tif"
merged_4326 = r"C:\Users\borrepa\Downloads\MDT02_pinofranqueado_merged_4326.tif"
merged_utm = r"C:\Users\borrepa\Downloads\MDT02_pinofranqueado_merged.tif"

dst_crs = "EPSG:25829"
dst_res = 2.0  # metres

# Step 1: Merge tiles
print("Merging tiles...")
src1 = rasterio.open(east_tile)
src2 = rasterio.open(west_tile)
print(f"  East: {src1.bounds}")
print(f"  West: {src2.bounds}")

mosaic, mosaic_transform = merge([src1, src2])
print(f"  Merged shape: {mosaic.shape}")
elev = mosaic[mosaic > -9999]
print(f"  Elevation range: {elev.min():.1f} – {elev.max():.1f} m")

meta = src1.meta.copy()
meta.update({
    'height': mosaic.shape[1],
    'width': mosaic.shape[2],
    'transform': mosaic_transform,
})
with rasterio.open(merged_4326, 'w', **meta) as dst:
    dst.write(mosaic)
src1.close()
src2.close()
print(f"  Saved merged (4326): {merged_4326}")

# Step 2: Reproject to UTM 29N
print(f"\nReprojecting to {dst_crs} at {dst_res}m...")
with rasterio.open(merged_4326) as src:
    transform, width, height = calculate_default_transform(
        src.crs, dst_crs, src.width, src.height, *src.bounds,
        resolution=dst_res
    )
    meta = src.meta.copy()
    meta.update({
        'crs': dst_crs,
        'transform': transform,
        'width': width,
        'height': height,
    })
    with rasterio.open(merged_utm, 'w', **meta) as dst:
        for i in range(1, src.count + 1):
            reproject(
                source=rasterio.band(src, i),
                destination=rasterio.band(dst, i),
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=transform,
                dst_crs=dst_crs,
                resampling=Resampling.bilinear
            )

with rasterio.open(merged_utm) as src:
    print(f"  Output: {src.width}x{src.height}, res={src.res}")
    print(f"  Bounds: {src.bounds}")
    data = src.read(1)
    valid = data[data > -9999]
    print(f"  Elevation: {valid.min():.1f} – {valid.max():.1f} m")

print(f"\n✅ Done: {merged_utm}")
