"""
OMAP file parsing: symbol classification, coordinate extraction, Bezier tessellation.
"""
import numpy as np
import math

# ── ISOM 2017-2 Symbol Categories ──
ISOM_CATEGORIES = {
    'paved_road':   [501, 501.1, 501.2, 529.1],
    'road':         [502, 502.1, 502.2, 503, 503.1],
    'track':        [504, 504.1, 831, 833, 835, 837],
    'footpath':     [505, 505.1, 832, 834],
    'small_path':   [506, 506.1, 836, 838],
    'faint_path':   [507, 507.1, 508, 508.1, 509, 509.1],
    'open':         [401, 401.1, 402, 402.1, 403, 403.1],
    'semi_open':    [404, 404.1],
    'forest_good':  [405, 405.1],
    'forest_slow':  [406, 406.1, 407, 407.1, 408, 408.1, 418, 419, 420],
    'forest_walk':  [409, 409.1, 410, 410.1],
    'forest_fight': [411, 411.1, 412, 412.1, 413],
    'cultivated':   [414, 414.1, 415, 415.1, 416, 416.1],
    'lake':         [301, 301.1, 301.2, 302, 302.1, 303, 303.1, 314],
    'river':        [304, 305, 305.1, 306, 306.1],
    'marsh':        [307, 308, 309, 309.1, 309.2, 310, 310.1,
                     311, 311.1, 312, 312.1, 313, 313.1],
    'cliff':        [201, 201.1, 201.2, 201.3, 201.4,
                     202, 202.1, 202.2, 203, 203.1],
    'rock':         [204, 205, 206, 206.1, 207, 208, 209],
    'stony':        [210, 210.1, 211, 212, 213, 213.1],
    'bare_rock':    [214, 214.1],
    'fence':        [522, 524, 524.1, 844],
    'wall':         [513, 513.1, 519, 519.1, 520, 520.1, 521, 521.1],
    'hedge':        [535, 536],
    'building':     [526, 526.1, 527, 529, 529.2],
    'railway':      [515, 516, 517, 518, 518.1, 534],
    'bridge':       [512, 512.1],
    'tree':         [417, 417.1],
    'veg_boundary': [416.1],
    'contour':      list(range(101, 116)),
    'ignore':       [525, 528, 530, 531, 532, 538, 539, 540,
                     601, 601.2, 601.3, 602, 603, 603.1,
                     701, 702, 703, 704, 705, 706, 707, 708, 709,
                     709.1, 711, 711.1, 712, 713, 799, 840, 843,
                     901, 950, 999, 2],
}

CAT_COST = {
    'paved_road': 0.08, 'road': 0.08, 'track': 0.10,
    'footpath': 0.12, 'small_path': 0.15, 'faint_path': 0.20,
    'open': 0.20, 'semi_open': 0.28, 'forest_good': 0.30,
    'forest_slow': 0.40, 'forest_walk': 0.55, 'forest_fight': 0.80,
    'cultivated': 0.25, 'lake': 10.0, 'river': 5.0,
    'marsh': 0.70, 'cliff': 10.0, 'rock': 0.60,
    'stony': 0.35, 'bare_rock': 0.30, 'fence': 2.0,
    'wall': 5.0, 'hedge': 1.5, 'building': 10.0,
    'railway': 1.0, 'bridge': 0.10, 'tree': 0.35,
    'veg_boundary': None, 'contour': None, 'ignore': None,
}

# Build lookup dictionaries
SYM2CAT = {}
for cat, syms in ISOM_CATEGORIES.items():
    for s in syms:
        SYM2CAT[s] = cat
        SYM2CAT[int(s)] = cat

NAME_TO_CAT = {
    'wall': 'wall', 'stone wall': 'wall', 'ruined wall': 'wall',
    'ruined stone wall': 'wall', 'fence': 'fence', 'cliff': 'cliff',
    'boulder': 'rock', 'path': 'small_path', 'track': 'track',
    'road': 'road', 'building': 'building', 'lake': 'lake',
    'marsh': 'marsh', 'prominent large tree': 'tree',
    'prominent tree': 'tree', 'spring': 'marsh', 'waterhole': 'lake',
    'bridge': 'bridge', 'hedge': 'hedge',
}


def get_category(sym_info):
    """Robust category lookup: ISOM code first, then name fallback."""
    isom = sym_info.get('isom')
    if isom is not None:
        cat = SYM2CAT.get(isom)
        if cat is not None:
            return cat
        cat = SYM2CAT.get(int(isom))
        if cat is not None:
            return cat
        cat = SYM2CAT.get(float(int(isom)))
        if cat is not None:
            return cat
    name = sym_info.get('name', '').strip().lower()
    if name:
        if name in NAME_TO_CAT:
            return NAME_TO_CAT[name]
        for key, cat in NAME_TO_CAT.items():
            if key in name:
                return cat
    return None


def parse_omap_coords_text(text):
    """Parse OMAP coordinate text with Bezier curve tessellation.

    OpenOrienteering Mapper's native XML format stores coordinates as
    integers in units of 1/1000 mm, not plain mm -- dividing here converts
    to the mm used everywhere else (paper_ref_x/y, the georeferencing
    scale math in rasterize_omap). Without this, every object landed
    ~1000x further from the map's reference point than intended, pushing
    all of them outside the UTM grid bounds: the whole ACR (map-symbol)
    cost layer silently rasterized to nothing, falling back to the
    uniform default cost everywhere (confirmed via a 0.0% non-default
    pixel count on a real run).
    """
    UNITS_PER_MM = 1000.0
    raw_tokens = text.strip().replace(';', ' ').split()
    coords = []
    flags = []
    i = 0
    while i < len(raw_tokens):
        tok = raw_tokens[i]
        flag = 0
        while not tok.replace('-', '', 1).replace('.', '', 1).isdigit():
            try:
                flag = int(tok)
            except ValueError:
                pass
            i += 1
            if i >= len(raw_tokens):
                return coords, flags
            tok = raw_tokens[i]
        x = float(tok) / UNITS_PER_MM
        i += 1
        if i >= len(raw_tokens):
            break
        y = float(raw_tokens[i]) / UNITS_PER_MM
        i += 1
        coords.append((x, y))
        flags.append(flag)

    # Tessellate Bezier curves
    out_coords = []
    idx = 0
    while idx < len(coords):
        if idx + 3 < len(coords) and flags[idx + 1] == 1 and flags[idx + 2] == 2:
            p0 = np.array(coords[idx])
            p1 = np.array(coords[idx + 1])
            p2 = np.array(coords[idx + 2])
            p3 = np.array(coords[idx + 3])
            chord = np.linalg.norm(p3 - p0)
            n_seg = max(4, int(chord / 50))
            for t_i in range(n_seg + 1):
                t = t_i / n_seg
                pt = ((1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * p1 +
                      3 * (1 - t) * t ** 2 * p2 + t ** 3 * p3)
                out_coords.append((float(pt[0]), float(pt[1])))
            idx += 3
        else:
            out_coords.append(coords[idx])
            idx += 1
    return out_coords, flags
