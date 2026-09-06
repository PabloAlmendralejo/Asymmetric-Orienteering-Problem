"""
Post-processing statistics for the heuristic comparison (compare_heuristics.cpp).
Reads heuristic_comparison_results.csv (per-instance mean scores over 10 seeds
for GA/ACO/SA, deterministic score for Greedy) and reports, for each of SA vs
{Greedy, GA, ACO}: exact/asymptotic Wilcoxon signed-rank p-value (scipy, not
the C++ program's own normal-approximation z-score), Holm-Bonferroni-corrected
p-value across the 3 comparisons, and the matched-pairs rank-biserial
correlation as an effect size.
"""
import csv
import sys
from scipy.stats import wilcoxon

def load(path):
    rows = []
    with open(path, newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            rows.append({k: (v if k == 'name' else float(v)) for k, v in r.items()})
    return rows

def rank_biserial(sa, other):
    diffs = [a - b for a, b in zip(sa, other) if abs(a - b) > 1e-9]
    ranks = sorted(range(len(diffs)), key=lambda i: abs(diffs[i]))
    rank_of = {}
    for pos, i in enumerate(ranks):
        rank_of[i] = pos + 1
    w_plus = sum(rank_of[i] for i in range(len(diffs)) if diffs[i] > 0)
    w_minus = sum(rank_of[i] for i in range(len(diffs)) if diffs[i] < 0)
    total = w_plus + w_minus
    return (w_plus - w_minus) / total if total > 0 else 0.0

def holm_bonferroni(pvals):
    order = sorted(range(len(pvals)), key=lambda i: pvals[i])
    m = len(pvals)
    adjusted = [0.0] * m
    running_max = 0.0
    for rank, i in enumerate(order):
        val = min(1.0, (m - rank) * pvals[i])
        running_max = max(running_max, val)
        adjusted[i] = running_max
    return adjusted

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'heuristic_comparison_results.csv'
    rows = load(path)
    n = len(rows)
    sa = [r['sa_mean'] for r in rows]
    methods = {
        'Greedy': [r['greedy'] for r in rows],
        'GA': [r['ga_mean'] for r in rows],
        'ACO': [r['aco_mean'] for r in rows],
    }

    raw_pvals = []
    stats = []
    for name, scores in methods.items():
        res = wilcoxon(sa, scores, zero_method='wilcox', method='auto')
        r = rank_biserial(sa, scores)
        raw_pvals.append(res.pvalue)
        stats.append((name, res.statistic, res.pvalue, r))

    adjusted = holm_bonferroni(raw_pvals)

    print(f"n = {n} instances")
    print(f"{'Method':<8} {'W':>10} {'p (raw)':>12} {'p (Holm)':>12} {'effect r':>10}")
    for (name, W, p, r), p_adj in zip(stats, adjusted):
        print(f"{name:<8} {W:>10.1f} {p:>12.3e} {p_adj:>12.3e} {r:>10.3f}")

if __name__ == '__main__':
    main()
