"""
Cost functions: Minetti metabolic model, fatigue, ISOM weights, base cost construction.
"""
import numpy as np

# ── ISOM 2017-2 Terrain Weights ──
WT = {
    'paved_road': 0.10, 'road': 0.10, 'footpath': 0.10,
    'small_footpath': 0.15, 'less_distinct_path': 0.20,
    'open_land': 0.20, 'rough_open': 0.25,
    'rough_open_bushes': 0.45, 'rough_open_thickets': 0.60,
    'forest': 0.30, 'slow_run': 0.35,
    'walk': 0.40, 'walk_good_vis': 0.45,
    'fight': 0.75, 'marsh': 0.30, 'hedge': 0.40, 'cliff': 0.70,
    'water': 10.0, 'building': 10.0,
    'impassable_cliff': 10.0, 'out_of_bounds': 10.0,
}

# ── Minetti et al. 2002 polynomial ──
_MINETTI_COEFFS = np.array([280.5, -58.7, -76.8, 51.9, 19.6, 2.5])
_MINETTI_FLAT = np.polyval(_MINETTI_COEFFS, 0.0)


def minetti_cost(slope_signed):
    """Metabolic cost (J/kg/m) as a function of slope gradient (rise/run)."""
    i = np.clip(slope_signed, -0.45, 0.45)
    return np.polyval(_MINETTI_COEFFS, i)


def minetti_factor(slope_signed):
    """Cost multiplier relative to flat terrain."""
    cost = minetti_cost(slope_signed)
    cost = np.maximum(cost, 0.5)
    return cost / _MINETTI_FLAT


def derive_mu():
    """Distance-to-elevation fatigue weight, derived from the Minetti curve.

    mu converts 1 metre of horizontal distance into "vertical-metre
    equivalent" fatigue units, for combining with gain/loss (phi_plus /
    phi_minus) in the asymmetric fatigue model's psi_ij = gain - rho*loss
    + mu*dist (Sec. 3.6 rework). Derived, not fitted: mu = C(0) / (C(i*)/i*),
    where C(i) is the Minetti metabolic cost (J/kg/m) at gradient i and i*
    is the gradient that minimises cost per VERTICAL metre climbed (the
    most energy-efficient climbing angle implied by the same curve already
    used for the base cost matrix). This keeps the distance term physically
    grounded in the same model as the rest of the cost function rather than
    a free/fitted parameter.
    """
    from scipy.optimize import minimize_scalar
    c0 = float(minetti_cost(0.0))

    def cost_per_vert_m(i):
        return float(minetti_cost(i)) / i

    res = minimize_scalar(cost_per_vert_m, bounds=(0.01, 0.45), method='bounded')
    cost_per_vert = cost_per_vert_m(res.x)
    return c0 / cost_per_vert


def fatigue_multiplier(elapsed_fraction, fatigue_rate, terrain_difficulty=1.0):
    """Linear fatigue model: cost increases with elapsed time."""
    eff = elapsed_fraction * terrain_difficulty
    return 1.0 + fatigue_rate * eff


def effective_budget(budget, fatigue_rate):
    """Approximate effective budget accounting for average fatigue."""
    avg_fatigue = 1.0 + fatigue_rate / 2.0
    return budget / avg_fatigue


def build_base_cost_from_omap(cost_omap, valid):
    """Build base cost grid from OMAP symbol costs only."""
    base = cost_omap.copy().astype(np.float32)
    base[~valid] = WT['out_of_bounds']
    base = np.maximum(base, 0.05)
    return base


def build_display_cost(base, slope_mag):
    """Display cost grid: base cost × Minetti slope factor."""
    sf = minetti_factor(slope_mag)
    cg = base * sf
    cg[slope_mag > 0.7] = np.maximum(cg[slope_mag > 0.7], WT['impassable_cliff'])
    return cg


def cost_to_veg_proxy(cost):
    """Convert cost grid to vegetation proxy [0, 1]."""
    low, high = 0.08, 0.75
    veg = np.clip(1.0 - (cost - low) / (high - low), 0.0, 1.0)
    veg[cost >= 5.0] = 0.0
    return veg.astype(np.float32)


def cost_str(cost_units, cell_meters, base_speed, reference_weight):
    """Format cost units as human-readable time string."""
    cost_to_seconds = cell_meters / (base_speed * reference_weight)
    s = cost_units * cost_to_seconds
    if s < 60:
        return f"{s:.0f}s"
    elif s < 3600:
        return f"{int(s // 60)}m{int(s % 60):02d}s"
    else:
        return f"{int(s // 3600)}h{int((s % 3600) // 60):02d}m"
