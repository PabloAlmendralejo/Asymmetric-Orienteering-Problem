"""
Generate parameterized asymmetric OP benchmark instances.

Produces JSON files compatible with the B&C solver (benchmark_solver.cpp).

Parameters varied:
  - n: number of controls (20, 30, 40, 50, 75, 100)
  - asymmetry: 0.0 (symmetric) to 0.5 (highly asymmetric)
  - fatigue_rate: 0.0, 0.1, 0.2, 0.3
  - budget_tightness: loose (~70%), medium (~50%), tight (~30%)

Usage:
    python generate_instances.py
    python generate_instances.py --output-dir instances/
"""
import numpy as np
import json
import os
import argparse
from itertools import product


def _nearest_neighbor_cost(cm, n):
    """Estimate cost to visit all nodes via nearest-neighbor heuristic."""
    total = n + 1
    visited = {0}
    cur = 0
    cost = 0.0
    for _ in range(n):
        best_j, best_c = -1, np.inf
        for j in range(total):
            if j not in visited and cm[cur][j] < best_c:
                best_j, best_c = j, cm[cur][j]
        if best_j < 0:
            break
        cost += best_c
        visited.add(best_j)
        cur = best_j
    cost += cm[cur][0]  # return to depot
    return cost


def generate_instance(n, asymmetry=0.25, fatigue_rate=0.2,
                      budget_tightness='medium', seed=42):
    """Generate a single asymmetric OP instance."""
    rng = np.random.RandomState(seed)
    total = n + 1

    # Place depot at center, controls randomly in [0, 1000]x[0, 1000]
    coords = np.zeros((total, 2))
    coords[0] = [500, 500]
    coords[1:, 0] = rng.uniform(50, 950, n)
    coords[1:, 1] = rng.uniform(50, 950, n)

    # Scores: correlated with distance from depot + noise
    dist_from_depot = np.sqrt(np.sum((coords[1:] - coords[0]) ** 2, axis=1))
    max_dist = max(dist_from_depot.max(), 1.0)
    raw_scores = 10 + 40 * (dist_from_depot / max_dist) + rng.uniform(-5, 5, n)
    pts = np.zeros(total)
    pts[1:] = (np.clip(raw_scores, 10, 50) // 10) * 10

    # Asymmetric cost matrix
    cm = np.zeros((total, total))
    for i in range(total):
        for j in range(total):
            if i == j:
                continue
            dx = coords[j, 0] - coords[i, 0]
            dy = coords[j, 1] - coords[i, 1]
            base_dist = np.sqrt(dx ** 2 + dy ** 2)

            # Directional asymmetry: "uphill" (positive dy) costs more
            slope_factor = dy / max(base_dist, 1e-9)
            asym_mult = 1.0 + asymmetry * slope_factor

            # Terrain noise (10%)
            terrain_noise = 1.0 + rng.uniform(-0.1, 0.1)
            cm[i, j] = base_dist * max(asym_mult * terrain_noise, 0.2)

    # Budget from nearest-neighbor estimate
    nn_cost = _nearest_neighbor_cost(cm, n)
    tightness_map = {'loose': 0.70, 'medium': 0.50, 'tight': 0.30}
    fraction = tightness_map[budget_tightness]
    bud_raw = nn_cost * fraction
    bud_eff = bud_raw / (1.0 + fatigue_rate / 2.0)

    # Measure actual asymmetry
    asym_vals = []
    for i in range(total):
        for j in range(i + 1, total):
            if cm[i, j] > 0 and cm[j, i] > 0:
                asym_vals.append(abs(cm[i, j] - cm[j, i]))
    mean_cost = np.mean(cm[cm > 0])
    actual_asym = 100 * np.mean(asym_vals) / mean_cost if mean_cost > 0 else 0

    return {
        'cm': cm.tolist(),
        'pts': pts.tolist(),
        'bud_eff': float(bud_eff),
        'bud_raw': float(bud_raw),
        'fatigue_rate': float(fatigue_rate),
    }, {
        'n': n,
        'asymmetry_param': asymmetry,
        'actual_asymmetry_pct': round(actual_asym, 1),
        'budget_tightness': budget_tightness,
        'seed': seed,
        'nn_cost': round(nn_cost, 1),
        'total_pts': int(sum(pts)),
    }


# ── Benchmark configurations ──

CONFIGS = {
    'nodes': [20, 30, 40, 50, 75, 100],
    'asymmetry': [0.0, 0.1, 0.25, 0.5],
    'fatigue_rate': [0.0, 0.1, 0.2, 0.3],
    'budget': ['loose', 'medium', 'tight'],
}

# Representative subset: vary one parameter at a time from a baseline
BASELINE = {'n': 40, 'asymmetry': 0.25, 'fatigue_rate': 0.2, 'budget': 'medium'}


def generate_all(output_dir='instances'):
    """Generate the full benchmark suite."""
    os.makedirs(output_dir, exist_ok=True)
    manifest = []
    instance_id = 0

    # 1. Vary node count (baseline asymmetry/fatigue/budget)
    for n in CONFIGS['nodes']:
        instance_id += 1
        name = f"bench_{instance_id:03d}_n{n}_a25_f20_med"
        data, meta = generate_instance(
            n, asymmetry=0.25, fatigue_rate=0.2,
            budget_tightness='medium', seed=1000 + instance_id)
        _save(output_dir, name, data, meta, manifest)

    # 2. Vary asymmetry (baseline n/fatigue/budget)
    for asym in CONFIGS['asymmetry']:
        instance_id += 1
        a_str = f"{int(asym * 100):02d}"
        name = f"bench_{instance_id:03d}_n40_a{a_str}_f20_med"
        data, meta = generate_instance(
            40, asymmetry=asym, fatigue_rate=0.2,
            budget_tightness='medium', seed=2000 + instance_id)
        _save(output_dir, name, data, meta, manifest)

    # 3. Vary fatigue rate (baseline n/asymmetry/budget)
    for fr in CONFIGS['fatigue_rate']:
        instance_id += 1
        f_str = f"{int(fr * 100):02d}"
        name = f"bench_{instance_id:03d}_n40_a25_f{f_str}_med"
        data, meta = generate_instance(
            40, asymmetry=0.25, fatigue_rate=fr,
            budget_tightness='medium', seed=3000 + instance_id)
        _save(output_dir, name, data, meta, manifest)

    # 4. Vary budget tightness (baseline n/asymmetry/fatigue)
    for bud in CONFIGS['budget']:
        instance_id += 1
        b_str = bud[0]  # l, m, t
        name = f"bench_{instance_id:03d}_n40_a25_f20_{b_str}"
        data, meta = generate_instance(
            40, asymmetry=0.25, fatigue_rate=0.2,
            budget_tightness=bud, seed=4000 + instance_id)
        _save(output_dir, name, data, meta, manifest)

    # 5. Extreme cases
    extremes = [
        (20, 0.0, 0.0, 'loose', 'symmetric_easy'),
        (100, 0.5, 0.3, 'tight', 'asymmetric_hard'),
        (50, 0.5, 0.0, 'medium', 'high_asym_no_fatigue'),
        (50, 0.0, 0.3, 'medium', 'symmetric_high_fatigue'),
    ]
    for n, asym, fr, bud, label in extremes:
        instance_id += 1
        name = f"bench_{instance_id:03d}_{label}"
        data, meta = generate_instance(
            n, asymmetry=asym, fatigue_rate=fr,
            budget_tightness=bud, seed=5000 + instance_id)
        _save(output_dir, name, data, meta, manifest)

    # Save manifest
    manifest_path = os.path.join(output_dir, 'manifest.json')
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    print(f"\n✅ Generated {len(manifest)} instances in {output_dir}/")
    print(f"   Manifest: {manifest_path}")
    return manifest


def _save(output_dir, name, data, meta, manifest):
    """Save instance JSON and append to manifest."""
    filename = f"op_input_{name}.json"
    filepath = os.path.join(output_dir, filename)
    with open(filepath, 'w') as f:
        json.dump(data, f)
    meta['filename'] = filename
    meta['name'] = name
    manifest.append(meta)
    n = meta['n']
    asym = meta['actual_asymmetry_pct']
    print(f"  {name}: n={n}, asym={asym:.1f}%, "
          f"pts={meta['total_pts']}, budget={meta['budget_tightness']}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate benchmark instances')
    parser.add_argument('--output-dir', default='instances',
                        help='Output directory for JSON files')
    args = parser.parse_args()
    generate_all(args.output_dir)
