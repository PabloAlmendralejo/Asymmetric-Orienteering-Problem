"""
Configuration for Pinofranqueado study area.
Rogaine map — Samper, November 2022.
"""
import os

# ── CRS and Projection ──
MAP_CRS = "EPSG:25829"  # ETRS89 / UTM zone 29N

# ── World file parameters (from Pinofranqueado_tif.tfw) ──
WORLD_A = 0.635
WORLD_B = 0.0
WORLD_C = 721261.3925
WORLD_D = 0.0
WORLD_E = -0.635
WORLD_F = 4467151.933

# ── Data files ──
DATA_DIR = os.environ.get('ORIENTEERING_DATA_DIR',
                          os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'pinofranqueado'))
TIF_PATH = os.path.join(DATA_DIR, 'Pinofranqueado_tif.tif')
MDT_PATH = os.path.join(DATA_DIR, 'MDT02_pinofranqueado_merged.tif')
OMAP_PATH = os.path.join(DATA_DIR, 'Pinofranqueado_omap.omap')
TERRAIN_CACHE = os.path.join(DATA_DIR, 'pinofranqueado_terrain.npz')

# ── Elevation tile merging ──
# The DEM sits across two tiles that must be merged and reprojected to UTM 29N.
# Run the merge script before the pipeline:
#   python merge_elevation.py pinofranqueado
MDT_TILE_EAST = os.path.join(DATA_DIR, 'PinoFranqueado_este.tif')
MDT_TILE_WEST = os.path.join(DATA_DIR, 'Pinofranquead_oeste.tif')

# ── Grid ──
TARGET_RESOLUTION = 2.0
DOWNSAMPLE = 8

# ── Race parameters ──
SEED = 42
NUM_CONTROLS = 40
BASE_SPEED = 2.5
REFERENCE_WEIGHT = 0.3
RACE_DURATION_HOURS = 1
FATIGUE_RATE = 0.20

# ── HCR Configuration ──
HCR_TRI_POWER = 1.0
HCR_PC_POWER = 0.5
HCR_SLO_POWER = 1.0
HCR_NORM_LOW = 0.1
HCR_NORM_HIGH = 1.0

# ── Derived constants ──
RACE_DURATION_SECONDS = RACE_DURATION_HOURS * 3600
COST_TO_SECONDS = TARGET_RESOLUTION / (BASE_SPEED * REFERENCE_WEIGHT)
ROUTE_BUDGET = RACE_DURATION_SECONDS / COST_TO_SECONDS

# ── JSON export distributions ──
DISTRIBUTIONS = {
    'standard':      {'func': 'standard',      'seed': 42,  'num': 40},
    'clustered':     {'func': 'clustered',     'seed': 123, 'num': 35},
    'ring':          {'func': 'ring',          'seed': 77,  'num': 30},
    'path_biased':   {'func': 'path_biased',   'seed': 99,  'num': 35},
    'elev_biased':   {'func': 'elev_biased',   'seed': 17,  'num': 35},
    'sparse_far':    {'func': 'sparse_far',    'seed': 31,  'num': 25},
    'mixed_density': {'func': 'mixed_density', 'seed': 55,  'num': 40},
}
