"""Check how much the distance term (mu*dist) contributes to the fatigue
state F relative to the elevation term (gain - rho*loss), for real routes
found by the solver. This tests whether the paper's lambda anchor
calibration (which approximates F by cumulative elevation gain D+ alone,
ignoring the distance term) is a reasonable simplification.
"""
import glob
import json
import os

INSTANCES_DIR = os.path.join(os.path.dirname(__file__), 'instances')


def analyze_route(inp, route):
    gain = inp['gain']
    loss = inp['loss']
    dist = inp['dist']
    rho = inp.get('rho_default', 0.5)
    mu = inp.get('mu_default', 0.0)

    full = [0] + list(route) + [0]
    F = 0.0
    total_gain = 0.0     # sum of phi+ (raw, before clipping interacts)
    total_rho_loss = 0.0  # sum of rho*phi-
    total_mu_dist = 0.0   # sum of mu*delta
    clip_events = 0
    for k in range(len(full) - 1):
        i, j = full[k], full[k + 1]
        g = gain[i][j]
        l = loss[i][j]
        d = dist[i][j]
        psi = g - rho * l + mu * d
        total_gain += g
        total_rho_loss += rho * l
        total_mu_dist += mu * d
        new_F = F + psi
        if new_F < 0:
            clip_events += 1
        F = max(0.0, new_F)

    # Unclipped net elevation-only contribution vs distance-only contribution
    elev_net = total_gain - total_rho_loss
    return {
        'n_legs': len(full) - 1,
        'final_F': F,
        'total_gain': total_gain,
        'total_rho_loss': total_rho_loss,
        'elev_net': elev_net,
        'total_mu_dist': total_mu_dist,
        'clip_events': clip_events,
        'mu': mu,
        'rho': rho,
        'dist_share_of_unclipped_sum': total_mu_dist / (elev_net + total_mu_dist) if (elev_net + total_mu_dist) != 0 else float('nan'),
    }


if __name__ == '__main__':
    results = []
    for path in sorted(glob.glob(os.path.join(INSTANCES_DIR, 'op_output_*_lam0.000000.json'))):
        base = os.path.basename(path)
        in_name = base.replace('op_output_', 'op_input_').replace('_lam0.000000', '')
        in_path = os.path.join(INSTANCES_DIR, in_name)
        if not os.path.exists(in_path):
            print(f"  (no matching input for {base}, skipping)")
            continue
        inp = json.load(open(in_path))
        out = json.load(open(path))
        route = out.get('bnc_sa', {}).get('route') or out.get('sa', {}).get('route')
        if not route:
            continue
        r = analyze_route(inp, route)
        r['instance'] = in_name
        results.append(r)

    print(f"{'instance':<45} {'legs':>5} {'final_F':>9} {'elev_net':>10} {'mu*dist':>10} {'dist_share%':>11} {'clips':>6}")
    print("-" * 100)
    for r in results:
        print(f"{r['instance']:<45} {r['n_legs']:>5} {r['final_F']:>9.1f} {r['elev_net']:>10.1f} "
              f"{r['total_mu_dist']:>10.1f} {100*r['dist_share_of_unclipped_sum']:>10.1f}% {r['clip_events']:>6}")

    if results:
        import statistics as st
        shares = [r['dist_share_of_unclipped_sum'] for r in results]
        print(f"\nMean distance-term share of (elev_net + mu*dist): {100*st.mean(shares):.1f}%")
        print(f"Median: {100*st.median(shares):.1f}%  Min: {100*min(shares):.1f}%  Max: {100*max(shares):.1f}%")
