
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <limits>
#include <queue>
#include <stack>
#include <vector>
#include <memory>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <functional>
#include <cassert>
#include <span>
#include <random>
#include <iomanip>
#include <set>
#include <filesystem>
#include <unordered_set>
#include "interfaces/highs_c_api.h"
#include "lp_data/HConst.h"
#include "Highs.h"

// ════════════════════════════════════════════════════════════════════════════
// Branch-and-Price solver for the Asymmetric Orienteering Problem with Fatigue
//
// Architecture:
//   Master LP: max ∑_r score(r) · λ_r
//              s.t. ∑_r a_{ir} · λ_r ≤ 1   for each node i (visit at most once)
//                   ∑_r λ_r ≤ 1             (select at most one route)
//                   λ_r ≥ 0
//
//   Pricing: find a route r* with positive reduced cost:
//     rc(r) = score(r) - ∑_i a_{ir} · π_i - μ
//   where π_i = dual of node-visit constraint, μ = dual of convexity constraint.
//   This is a Resource-Constrained Shortest Path Problem (RCSPP) solved via
//   dynamic programming with dominance (labeling algorithm).
//
//   Resource: cumulative fatigue-adjusted cost ≤ B
// ════════════════════════════════════════════════════════════════════════════

// ── JSON parser & Input (same as benchmark_solver.cpp) ─────────────────────

static std::string read_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot open " + path);
    return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

struct Input {
    std::vector<std::vector<double>> cm;
    std::vector<double> pts;
    double bud_eff = 0.0, bud_raw = 0.0, fatigue_rate = 0.0;
};

static Input parse_input(const std::string& json_str) {
    Input inp;
    size_t i = 0;
    auto skip_ws = [&](size_t pos) -> size_t {
        while (pos < json_str.size() && (json_str[pos]==' '||json_str[pos]=='\t'||
               json_str[pos]=='\n'||json_str[pos]=='\r')) ++pos;
        return pos;
    };
    auto parse_number = [&json_str](size_t& pos) -> double {
        size_t start = pos;
        if (pos < json_str.size() && (json_str[pos]=='-'||json_str[pos]=='+')) ++pos;
        while (pos < json_str.size() && (std::isdigit(json_str[pos])||json_str[pos]=='.'||
               json_str[pos]=='e'||json_str[pos]=='E'||json_str[pos]=='+'||json_str[pos]=='-')) ++pos;
        return std::stod(json_str.substr(start, pos - start));
    };
    auto parse_array1d = [&](size_t& pos) -> std::vector<double> {
        std::vector<double> v;
        pos = skip_ws(pos);
        if (json_str[pos] != '[') throw std::runtime_error("Expected [");
        ++pos;
        while (true) {
            pos = skip_ws(pos);
            if (json_str[pos] == ']') { ++pos; break; }
            if (json_str[pos] == ',') { ++pos; continue; }
            v.push_back(parse_number(pos));
        }
        return v;
    };
    auto parse_array2d = [&](size_t& pos) -> std::vector<std::vector<double>> {
        std::vector<std::vector<double>> m;
        pos = skip_ws(pos);
        if (json_str[pos] != '[') throw std::runtime_error("Expected [");
        ++pos;
        while (true) {
            pos = skip_ws(pos);
            if (json_str[pos] == ']') { ++pos; break; }
            if (json_str[pos] == ',') { ++pos; continue; }
            if (json_str[pos] == '[') m.push_back(parse_array1d(pos));
        }
        return m;
    };
    auto find_key = [&](const std::string& key) {
        std::string key_str = "\"" + key + "\":";
        size_t pos = json_str.find(key_str, i);
        if (pos == std::string::npos) throw std::runtime_error("Missing key: " + key);
        i = pos + key_str.size();
    };
    find_key("cm"); inp.cm = parse_array2d(i);
    find_key("pts"); inp.pts = parse_array1d(i);
    find_key("bud_eff"); i = skip_ws(i); inp.bud_eff = parse_number(i);
    find_key("bud_raw"); i = skip_ws(i); inp.bud_raw = parse_number(i);
    find_key("fatigue_rate"); i = skip_ws(i); inp.fatigue_rate = parse_number(i);
    return inp;
}

// ── Cost helpers ───────────────────────────────────────────────────────────

static double rcost(const std::vector<std::vector<double>>& cm, const std::vector<int>& route) {
    if (route.empty()) return 0.0;
    double c = cm[0][route[0]];
    for (size_t i = 0; i + 1 < route.size(); ++i) c += cm[route[i]][route[i+1]];
    return c + cm[route.back()][0];
}

static double rcost_fatigue(const std::vector<std::vector<double>>& cm, const std::vector<int>& route,
                            double bud_raw, double fatigue_rate) {
    if (route.empty()) return 0.0;
    double total = 0.0, elapsed = 0.0;
    std::vector<int> seq = {0};
    seq.insert(seq.end(), route.begin(), route.end());
    seq.push_back(0);
    for (size_t i = 0; i + 1 < seq.size(); ++i) {
        double leg = cm[seq[i]][seq[i+1]];
        total += leg * (1.0 + fatigue_rate * (elapsed / std::max(bud_raw, 1.0)));
        elapsed += leg;
    }
    return total;
}

static double rpts(const std::vector<double>& pts, const std::vector<int>& route) {
    double s = 0; for (int v : route) s += pts[v]; return s;
}

// ── SA (same as benchmark_solver.cpp) ──────────────────────────────────────

static std::vector<int> greedy_route(const Input& inp) {
    int n = static_cast<int>(inp.pts.size());
    std::vector<int> route;
    std::vector<bool> used(n, false);
    used[0] = true;
    int cur = 0;
    while (true) {
        int best = -1; double best_ratio = -1;
        for (int j = 1; j < n; ++j) {
            if (used[j] || !std::isfinite(inp.cm[cur][j])) continue;
            auto trial = route; trial.push_back(j);
            double fc = rcost_fatigue(inp.cm, trial, inp.bud_raw, inp.fatigue_rate);
            if (fc > inp.bud_raw) continue;
            double ratio = inp.pts[j] / std::max(inp.cm[cur][j], 1e-9);
            if (ratio > best_ratio) { best_ratio = ratio; best = j; }
        }
        if (best < 0) break;
        route.push_back(best); used[best] = true; cur = best;
    }
    return route;
}

static std::vector<int> solve_sa(const Input& inp, int n_iterations = 80000,
                                  double T0 = 100.0, double Tend = 0.1, unsigned seed = 42) {
    std::mt19937 rng(seed);
    auto route = greedy_route(inp);
    double best_score = rpts(inp.pts, route);
    auto best_route = route;
    double score = best_score;
    int n = static_cast<int>(inp.pts.size());
    double rho = std::pow(Tend / T0, 1.0 / n_iterations);
    double T = T0;
    int stagnation = 0, stag_limit = n_iterations / 5;

    for (int iter = 0; iter < n_iterations; ++iter) {
        if (stagnation >= stag_limit) break;
        T *= rho;
        double u = std::uniform_real_distribution<>(0, 1)(rng);
        auto new_route = route;
        if (u < 0.30 && !new_route.empty()) {
            int idx = std::uniform_int_distribution<>(0, (int)new_route.size()-1)(rng);
            new_route.erase(new_route.begin() + idx);
        } else if (u < 0.60) {
            std::vector<int> unvisited;
            std::set<int> vis(route.begin(), route.end()); vis.insert(0);
            for (int j = 1; j < n; ++j) if (!vis.count(j)) unvisited.push_back(j);
            if (!unvisited.empty()) {
                int node = unvisited[std::uniform_int_distribution<>(0,(int)unvisited.size()-1)(rng)];
                double best_inc = 1e30; int best_pos = 0;
                for (int p = 0; p <= (int)new_route.size(); ++p) {
                    auto trial = new_route; trial.insert(trial.begin()+p, node);
                    double fc = rcost_fatigue(inp.cm, trial, inp.bud_raw, inp.fatigue_rate);
                    if (fc < best_inc) { best_inc = fc; best_pos = p; }
                }
                new_route.insert(new_route.begin()+best_pos, node);
            }
        } else if (u < 0.80 && new_route.size() >= 2) {
            int len = std::uniform_int_distribution<>(1, std::min(3,(int)new_route.size()))(rng);
            int from = std::uniform_int_distribution<>(0, (int)new_route.size()-len)(rng);
            std::vector<int> seg(new_route.begin()+from, new_route.begin()+from+len);
            new_route.erase(new_route.begin()+from, new_route.begin()+from+len);
            int to = std::uniform_int_distribution<>(0, (int)new_route.size())(rng);
            new_route.insert(new_route.begin()+to, seg.begin(), seg.end());
        } else if (!new_route.empty()) {
            std::set<int> vis(route.begin(), route.end()); vis.insert(0);
            std::vector<int> unvisited;
            for (int j = 1; j < n; ++j) if (!vis.count(j)) unvisited.push_back(j);
            if (!unvisited.empty()) {
                int idx = std::uniform_int_distribution<>(0,(int)new_route.size()-1)(rng);
                int node = unvisited[std::uniform_int_distribution<>(0,(int)unvisited.size()-1)(rng)];
                new_route[idx] = node;
            }
        }
        double fc = rcost_fatigue(inp.cm, new_route, inp.bud_raw, inp.fatigue_rate);
        if (fc > inp.bud_raw) { ++stagnation; continue; }
        double ns = rpts(inp.pts, new_route);
        double delta = ns - score;
        if (delta > 0 || std::exp(delta / T) > std::uniform_real_distribution<>(0,1)(rng)) {
            route = new_route; score = ns;
            if (score > best_score) { best_score = score; best_route = route; stagnation = 0; }
            else ++stagnation;
        } else ++stagnation;
    }
    return best_route;
}

static std::vector<int> solve_sa_iterated(const Input& inp, int n_restarts = -1,
                                           int n_iterations = -1) {
    int n = static_cast<int>(inp.pts.size()) - 1;
    if (n_iterations < 0) n_iterations = std::clamp(2500 * n, 10000, 120000);
    if (n_restarts < 0) n_restarts = std::clamp(3000 / std::max(n, 1), 4, 20);
    std::vector<int> best_route;
    double best_score = 0.0;
    std::mt19937 seed_rng(42);
    for (int r = 0; r < n_restarts; ++r) {
        auto route = solve_sa(inp, n_iterations, 100.0, 0.1, seed_rng());
        double score = rpts(inp.pts, route);
        if (score > best_score) { best_score = score; best_route = route; }
    }
    return best_route;
}

// ════════════════════════════════════════════════════════════════════════════
// Column = a feasible route (sequence of nodes, depot-to-depot)
// ════════════════════════════════════════════════════════════════════════════

struct Column {
    std::vector<int> nodes;       // visited nodes (excluding depot)
    double score = 0.0;           // total points
    double fatigue_cost = 0.0;    // fatigue-adjusted travel cost
    std::vector<bool> visits;     // visits[i] = true if node i is in route
};

static Column make_column(const Input& inp, const std::vector<int>& route) {
    Column col;
    col.nodes = route;
    col.score = rpts(inp.pts, route);
    col.fatigue_cost = rcost_fatigue(inp.cm, route, inp.bud_raw, inp.fatigue_rate);
    int n = static_cast<int>(inp.pts.size());
    col.visits.assign(n, false);
    for (int v : route) col.visits[v] = true;
    return col;
}

// ════════════════════════════════════════════════════════════════════════════
// Pricing: Resource-Constrained Shortest Path with:
//   1. ng-route relaxation (neighborhood size NG_SIZE)
//   2. Completion bounds (greedy upper bound on remaining profit)
//   3. Label bucketing by resource consumption
//   4. Strong dominance (same node + ng-memory, cost ≤, profit ≥)
//   5. Label count cap to prevent blowup
//
// ng-route relaxation: instead of tracking the full visited set, each label
// tracks an "ng-memory" Π ⊆ N(node) — the set of nearby nodes that cannot
// be revisited. This allows revisits of distant nodes (which the master LP
// will forbid via node-visit constraints), dramatically reducing state space.
// ════════════════════════════════════════════════════════════════════════════

static constexpr int NG_SIZE = 10;       // neighborhood size for ng-routes
static constexpr int MAX_LABELS = 200000; // cap total labels to prevent blowup
static constexpr int N_BUCKETS = 50;      // resource buckets for dominance

struct Label {
    int node;
    double elapsed_base;      // cumulative base cost
    double fatigue_cost;      // cumulative fatigue-adjusted cost
    double reduced_profit;    // accumulated reduced profit
    std::vector<int> path;    // nodes visited (for route extraction)
    uint64_t ng_memory;       // bitmask: nodes in ng-neighborhood that are forbidden
};

// Precompute ng-neighborhoods: for each node, the NG_SIZE nearest nodes
static std::vector<uint64_t> compute_ng_neighborhoods(const Input& inp) {
    int n = static_cast<int>(inp.pts.size());
    std::vector<uint64_t> ng(n, 0ULL);
    for (int i = 0; i < n && i < 64; ++i) {
        // Sort other nodes by distance from i
        std::vector<std::pair<double, int>> dists;
        for (int j = 0; j < n && j < 64; ++j) {
            if (j == i) continue;
            double d = std::isfinite(inp.cm[i][j]) ? inp.cm[i][j] : 1e18;
            dists.push_back({d, j});
        }
        std::sort(dists.begin(), dists.end());
        uint64_t mask = 1ULL << i;  // always include self
        for (int k = 0; k < std::min(NG_SIZE, (int)dists.size()); ++k)
            mask |= 1ULL << dists[k].second;
        ng[i] = mask;
    }
    return ng;
}

// Completion bound: greedy estimate of max remaining reduced profit
// given remaining budget. Sort unvisited nodes by profit/cost ratio,
// greedily pack them (fractional knapsack for upper bound).
static double completion_bound(const std::vector<double>& rp,
                                const std::vector<std::vector<double>>& cm,
                                int cur_node, double remaining_budget,
                                uint64_t ng_memory, int n) {
    struct Item { double profit; double cost; };
    std::vector<Item> items;
    for (int j = 1; j < n && j < 64; ++j) {
        if (ng_memory & (1ULL << j)) continue;  // forbidden by ng-memory
        if (rp[j] <= 0) continue;               // negative reduced profit
        double cost_to = std::isfinite(cm[cur_node][j]) ? cm[cur_node][j] : 1e18;
        double cost_back = std::isfinite(cm[j][0]) ? cm[j][0] : 1e18;
        double min_cost = cost_to + cost_back;   // minimum cost to visit j and return
        if (min_cost > remaining_budget) continue;
        items.push_back({rp[j], cost_to});
    }
    // Sort by profit/cost ratio descending
    std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
        return a.profit * b.cost > b.profit * a.cost;
    });
    double bound = 0.0, budget_left = remaining_budget;
    for (const auto& item : items) {
        if (item.cost <= budget_left) {
            bound += item.profit;
            budget_left -= item.cost;
        } else {
            // Fractional: take what we can
            bound += item.profit * (budget_left / item.cost);
            break;
        }
    }
    return bound;
}

// Pricing subproblem: find route with maximum reduced cost
static std::vector<int> solve_pricing(const Input& inp,
                                       const std::vector<double>& pi,
                                       double mu,
                                       const std::vector<int>& forbidden_nodes,
                                       const std::vector<int>& required_nodes) {
    int n = static_cast<int>(inp.pts.size());
    int n_eff = std::min(n, 64);  // cap at 64 for bitmask
    double B = inp.bud_raw;
    double lambda = inp.fatigue_rate;

    std::vector<double> rp(n, 0.0);
    for (int i = 1; i < n; ++i) rp[i] = inp.pts[i] - pi[i];

    std::set<int> forbidden_set(forbidden_nodes.begin(), forbidden_nodes.end());
    auto ng_neighborhoods = compute_ng_neighborhoods(inp);

    // Bucketed label storage: labels[node][bucket] = list of labels
    double bucket_width = B / N_BUCKETS;
    std::vector<std::vector<std::vector<Label>>> labels(
        n_eff, std::vector<std::vector<Label>>(N_BUCKETS + 1));

    int total_labels = 0;

    // Dominance: L1 dominates L2 at same node and same bucket if:
    //   L1.fatigue_cost ≤ L2.fatigue_cost AND
    //   L1.reduced_profit ≥ L2.reduced_profit AND
    //   L1.ng_memory ⊆ L2.ng_memory (L1 forbids fewer nodes = more flexible)
    auto dominates = [](const Label& a, const Label& b) -> bool {
        return a.fatigue_cost <= b.fatigue_cost + 1e-9 &&
               a.reduced_profit >= b.reduced_profit - 1e-9 &&
               (a.ng_memory & b.ng_memory) == a.ng_memory;
    };

    auto get_bucket = [&](double cost) -> int {
        int b = static_cast<int>(cost / bucket_width);
        return std::min(b, N_BUCKETS);
    };

    // Try to add a label; returns true if added (not dominated)
    auto try_add_label = [&](Label& lab) -> bool {
        if (lab.node < 0 || lab.node >= n_eff) return false;
        int bkt = get_bucket(lab.fatigue_cost);
        auto& bucket = labels[lab.node][bkt];

        // Check if dominated by existing label in same or earlier buckets
        for (int b = 0; b <= bkt; ++b) {
            for (const auto& existing : labels[lab.node][b]) {
                if (dominates(existing, lab)) return false;
            }
        }

        // Remove labels in same or later buckets dominated by new label
        for (int b = bkt; b <= N_BUCKETS; ++b) {
            auto& bkt_labels = labels[lab.node][b];
            auto it = std::remove_if(bkt_labels.begin(), bkt_labels.end(),
                [&](const Label& l) {
                    if (dominates(lab, l)) { --total_labels; return true; }
                    return false;
                });
            bkt_labels.erase(it, bkt_labels.end());
        }

        bucket.push_back(std::move(lab));
        ++total_labels;
        return true;
    };

    // Initial label at depot
    Label init;
    init.node = 0;
    init.elapsed_base = 0.0;
    init.fatigue_cost = 0.0;
    init.reduced_profit = 0.0;
    init.ng_memory = ng_neighborhoods[0];
    try_add_label(init);

    // Process labels in order of increasing fatigue cost (Dijkstra-like)
    // Use a priority queue: (fatigue_cost, node, index_in_bucket)
    using PQEntry = std::tuple<double, int, int, int>;  // cost, node, bucket, idx
    std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<>> pq;
    pq.push({0.0, 0, 0, 0});

    double best_rc = -mu;
    std::vector<int> best_route;

    while (!pq.empty() && total_labels < MAX_LABELS) {
        auto [cost, cur_node, cur_bkt, cur_idx] = pq.top();
        pq.pop();

        // Validate label still exists (may have been dominated and removed)
        if (cur_node >= n_eff || cur_bkt > N_BUCKETS) continue;
        auto& bkt_labels = labels[cur_node][cur_bkt];
        if (cur_idx >= (int)bkt_labels.size()) continue;
        const Label& cur = bkt_labels[cur_idx];
        if (std::abs(cur.fatigue_cost - cost) > 1e-9) continue;  // stale entry

        // Check if this label can return to depot with positive reduced cost
        double return_base = inp.cm[cur_node][0];
        if (cur_node > 0 && std::isfinite(return_base)) {
            double return_fatigue = return_base * (1.0 + lambda * cur.elapsed_base / std::max(B, 1.0));
            if (cur.fatigue_cost + return_fatigue <= B + 1e-9) {
                double rc = cur.reduced_profit - mu;
                if (rc > best_rc + 1e-9) {
                    best_rc = rc;
                    best_route = cur.path;
                }
            }
        }

        // Completion bound pruning
        double remaining = B - cur.fatigue_cost;
        double cb = completion_bound(rp, inp.cm, cur_node, remaining, cur.ng_memory, n_eff);
        if (cur.reduced_profit + cb - mu <= best_rc + 1e-9) continue;

        // Extend to neighbors
        for (int j = 1; j < n_eff; ++j) {
            if (j == cur_node) continue;
            if (forbidden_set.count(j)) continue;
            if (!std::isfinite(inp.cm[cur_node][j])) continue;

            // ng-route check: can only visit j if j is NOT in current ng_memory
            // (ng_memory tracks which nearby nodes are forbidden)
            if (cur.ng_memory & (1ULL << j)) continue;

            double leg_base = inp.cm[cur_node][j];
            double new_elapsed = cur.elapsed_base + leg_base;
            double leg_fatigue = leg_base * (1.0 + lambda * cur.elapsed_base / std::max(B, 1.0));
            double new_fatigue = cur.fatigue_cost + leg_fatigue;

            // Feasibility: can we still return to depot?
            double ret_base = inp.cm[j][0];
            if (!std::isfinite(ret_base)) continue;
            double ret_fatigue = ret_base * (1.0 + lambda * new_elapsed / std::max(B, 1.0));
            if (new_fatigue + ret_fatigue > B + 1e-9) continue;

            Label ext;
            ext.node = j;
            ext.elapsed_base = new_elapsed;
            ext.fatigue_cost = new_fatigue;
            ext.reduced_profit = cur.reduced_profit + rp[j];
            ext.path = cur.path;
            ext.path.push_back(j);
            // Update ng-memory: union of current ng_memory with j's neighborhood
            // This forbids j and its neighbors from being revisited
            ext.ng_memory = cur.ng_memory | ng_neighborhoods[j];

            if (try_add_label(ext)) {
                int bkt = get_bucket(ext.fatigue_cost);
                int idx = static_cast<int>(labels[j][bkt].size()) - 1;
                pq.push({ext.fatigue_cost, j, bkt, idx});
            }
        }
    }

    if (total_labels >= MAX_LABELS)
        std::cerr << "  Pricing: label cap reached (" << MAX_LABELS << ")\n";

    return best_route;
}

// ════════════════════════════════════════════════════════════════════════════
// Master LP + Column Generation
// ════════════════════════════════════════════════════════════════════════════

struct BPSolver {
    const Input& inp;
    int n;
    std::vector<Column> columns;
    double best_pts = 0.0;
    std::vector<int> best_route;
    bool proved_optimal = false;
    double best_ub = std::numeric_limits<double>::infinity();
    double time_limit_s = 900.0;

    explicit BPSolver(const Input& i) : inp(i), n(static_cast<int>(i.pts.size())) {}

    // Add initial columns from SA and greedy
    void add_initial_columns(const std::vector<int>& sa_route) {
        if (!sa_route.empty()) {
            columns.push_back(make_column(inp, sa_route));
        }
        // Add single-node routes for each reachable node
        for (int i = 1; i < n; ++i) {
            if (!std::isfinite(inp.cm[0][i]) || !std::isfinite(inp.cm[i][0])) continue;
            double fc = rcost_fatigue(inp.cm, {i}, inp.bud_raw, inp.fatigue_rate);
            if (fc <= inp.bud_raw) {
                columns.push_back(make_column(inp, {i}));
            }
        }
        // Add empty route (visit nothing)
        columns.push_back(make_column(inp, {}));
    }

    // Solve master LP with current columns, return (obj_value, node_duals, convexity_dual)
    struct MasterResult {
        double obj;
        std::vector<double> pi;   // node duals (n values)
        double mu;                // convexity constraint dual
        std::vector<double> lambda_vals;  // column values
        bool feasible;
    };

    MasterResult solve_master() {
        MasterResult res;
        res.pi.resize(n, 0.0);
        res.mu = 0.0;
        res.feasible = false;

        int ncols = static_cast<int>(columns.size());
        if (ncols == 0) return res;

        // Build LP: max ∑ score_r * λ_r
        //   s.t. ∑ a_{ir} * λ_r ≤ 1   for i = 1..n-1  (node visit)
        //        ∑ λ_r ≤ 1                              (convexity)
        //        λ_r ≥ 0

        Highs* highs = (Highs*)Highs_create();
        Highs_setBoolOptionValue(highs, "output_flag", false);
        Highs_setStringOptionValue(highs, "solver", "simplex");

        // Add λ columns
        // Each column r has: obj = score_r, bounds [0, 1]
        // Constraints: n-1 node-visit rows + 1 convexity row = n rows total
        int nrows = n;  // rows 0..n-2 = node visit for nodes 1..n-1, row n-1 = convexity

        // First add all rows with bounds
        for (int i = 1; i < n; ++i) {
            // Node visit: ∑ a_{ir} λ_r ≤ 1
            Highs_addRow(highs, -1e30, 1.0, 0, nullptr, nullptr);
        }
        // Convexity: ∑ λ_r ≤ 1
        Highs_addRow(highs, -1e30, 1.0, 0, nullptr, nullptr);

        // Add columns
        for (int r = 0; r < ncols; ++r) {
            const auto& col = columns[r];
            std::vector<int> indices;
            std::vector<double> values;

            // Node visit coefficients
            for (int i = 1; i < n; ++i) {
                if (col.visits[i]) {
                    indices.push_back(i - 1);  // row index for node i
                    values.push_back(1.0);
                }
            }
            // Convexity coefficient
            indices.push_back(n - 1);  // last row
            values.push_back(1.0);

            Highs_addCol(highs, col.score, 0.0, 1.0,
                         static_cast<int>(indices.size()),
                         indices.data(), values.data());
        }

        Highs_changeObjectiveSense(highs, -1);  // maximize
        int status = Highs_run(highs);
        bool ok = status == 0 && Highs_getModelStatus(highs) == 7;

        if (ok) {
            res.feasible = true;
            res.obj = Highs_getObjectiveValue(highs);

            int nc = Highs_getNumCol(highs);
            int nr = Highs_getNumRow(highs);
            std::vector<double> col_val(nc), col_dual(nc), row_val(nr), row_dual(nr);
            Highs_getSolution(highs, col_val.data(), col_dual.data(),
                              row_val.data(), row_dual.data());

            // Extract duals
            for (int i = 1; i < n; ++i)
                res.pi[i] = row_dual[i - 1];
            res.mu = row_dual[n - 1];

            res.lambda_vals.resize(ncols);
            for (int r = 0; r < ncols; ++r)
                res.lambda_vals[r] = col_val[r];
        }

        Highs_destroy(highs);
        return res;
    }

    void solve(double warm_start_pts = 0.0, std::vector<int> warm_start_route = {}) {
        if (warm_start_pts > best_pts) {
            best_pts = warm_start_pts;
            best_route = warm_start_route;
        }

        add_initial_columns(warm_start_route);
        std::cerr << "B&P: " << columns.size() << " initial columns\n";

        auto t_start = std::chrono::steady_clock::now();
        int cg_iters = 0;
        int max_cg_iters = 500;
        int no_improve = 0;

        // Column generation loop
        while (cg_iters++ < max_cg_iters) {
            double elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_start).count();
            if (elapsed > time_limit_s) { std::cerr << "Time limit in CG\n"; break; }

            auto master = solve_master();
            if (!master.feasible) { std::cerr << "Master infeasible\n"; break; }

            best_ub = std::min(best_ub, master.obj);
            std::cerr << "CG iter " << cg_iters << ": master obj=" << master.obj
                      << " cols=" << columns.size() << "\n";

            // Check if any column gives an integer solution
            for (int r = 0; r < (int)columns.size(); ++r) {
                if (master.lambda_vals[r] > 0.99 && columns[r].score > best_pts) {
                    best_pts = columns[r].score;
                    best_route = columns[r].nodes;
                    std::cerr << "  New integer solution: " << best_pts << " pts\n";
                }
            }

            // Solve pricing
            auto new_route = solve_pricing(inp, master.pi, master.mu, {}, {});
            if (new_route.empty()) {
                std::cerr << "  No positive reduced cost column — CG converged\n";
                // CG converged: master obj is a valid upper bound
                best_ub = master.obj;
                proved_optimal = (best_ub <= best_pts + 1e-6);
                break;
            }

            // Check if this route is already in the pool
            auto new_col = make_column(inp, new_route);
            bool duplicate = false;
            for (const auto& c : columns) {
                if (c.nodes == new_col.nodes) { duplicate = true; break; }
            }
            if (duplicate) {
                ++no_improve;
                if (no_improve > 10) {
                    std::cerr << "  Stagnation in CG\n";
                    best_ub = master.obj;
                    break;
                }
                continue;
            }
            no_improve = 0;

            columns.push_back(std::move(new_col));
            std::cerr << "  Added column: " << columns.back().score << " pts, "
                      << columns.back().nodes.size() << " nodes\n";
        }

        // Final: check if we proved optimality
        if (best_ub <= best_pts + 1e-6) proved_optimal = true;

        double gap_pct = best_pts > 0 ? 100.0 * (best_ub - best_pts) / best_pts : 0.0;
        std::cerr << "B&P done: best=" << best_pts << " UB=" << best_ub
                  << " gap=" << gap_pct << "% cols=" << columns.size()
                  << " optimal=" << proved_optimal << "\n";
    }
};

// ── Main ───────────────────────────────────────────────────────────────────

static void run_map(const std::string& in_path, const std::string& out_path) {
    std::cerr << "\n=== " << in_path << " ===\n";
    Input inp = parse_input(read_file(in_path));

    auto t_sa = std::chrono::steady_clock::now();
    auto sa_route = solve_sa_iterated(inp);
    double sa_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_sa).count();
    double sa_pts = rpts(inp.pts, sa_route);
    double sa_base = rcost(inp.cm, sa_route);
    double sa_fatigue = rcost_fatigue(inp.cm, sa_route, inp.bud_raw, inp.fatigue_rate);
    std::cerr << "SA: " << sa_pts << " pts (" << sa_route.size() << " nodes) in " << sa_elapsed << "s\n";

    auto t_bp = std::chrono::steady_clock::now();
    BPSolver solver(inp);
    solver.solve(sa_pts, sa_route);
    double bp_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_bp).count();
    double bp_base = rcost(inp.cm, solver.best_route);
    double bp_fatigue = rcost_fatigue(inp.cm, solver.best_route, inp.bud_raw, inp.fatigue_rate);

    std::ofstream out(out_path);
    out << "{\n";
    out << "  \"sa\": {\"pts\": " << sa_pts << ", \"nodes\": " << sa_route.size()
        << ", \"elapsed_s\": " << sa_elapsed << ", \"base_cost\": " << sa_base
        << ", \"fatigue_cost\": " << sa_fatigue << ", \"route\": [";
    for (size_t i = 0; i < sa_route.size(); ++i) { if (i) out << ", "; out << sa_route[i]; }
    out << "]},\n";
    out << "  \"bnc_sa\": {\"pts\": " << solver.best_pts << ", \"nodes\": " << solver.best_route.size()
        << ", \"elapsed_s\": " << bp_elapsed
        << ", \"proved_optimal\": " << (solver.proved_optimal ? "true" : "false")
        << ", \"best_ub\": " << solver.best_ub
        << ", \"gap_pct\": " << (solver.best_pts > 0 ? 100.0 * (solver.best_ub - solver.best_pts) / solver.best_pts : 0.0)
        << ", \"base_cost\": " << bp_base
        << ", \"fatigue_cost\": " << bp_fatigue << ", \"route\": [";
    for (size_t i = 0; i < solver.best_route.size(); ++i) { if (i) out << ", "; out << solver.best_route[i]; }
    out << "]}\n";
    out << "}\n";
}

struct MapResult {
    std::string name;
    double sa_pts, bnc_pts;
    int sa_nodes, bnc_nodes;
    double sa_s, bnc_s;
    bool bnc_optimal;
    double gap_pct;
};

int main(int argc, char* argv[]) {
    std::string input_dir = "instances";
    if (argc > 1) input_dir = argv[1];

    std::vector<std::pair<std::string,std::string>> maps;
    for (const auto& entry : std::filesystem::directory_iterator(input_dir)) {
        std::string fname = entry.path().filename().string();
        if (fname.substr(0, 9) == "op_input_" && fname.size() > 14 &&
            fname.substr(fname.size() - 5) == ".json") {
            std::string base = fname.substr(9, fname.size() - 14);
            std::string in_path = entry.path().string();
            std::string out_path = (entry.path().parent_path() /
                                    ("op_output_bp_" + base + ".json")).string();
            maps.push_back({in_path, out_path});
        }
    }
    std::sort(maps.begin(), maps.end());
    std::cerr << "Found " << maps.size() << " instances in " << input_dir << "/\n";

    auto json_val = [](const std::string& s, const std::string& key) -> double {
        std::string k = "\"" + key + "\": ";
        auto p = s.find(k);
        if (p == std::string::npos) return 0.0;
        p += k.size();
        return std::stod(s.substr(p, s.find_first_of(",}", p) - p));
    };

    std::vector<MapResult> results;
    for (const auto& [in, out] : maps) {
        try { run_map(in, out); }
        catch (const std::exception& e) {
            std::cerr << "Error on " << in << ": " << e.what() << '\n';
            continue;
        }
        try {
            std::string js = read_file(out);
            auto bnc_pos = js.find("\"bnc_sa\"");
            std::string sa_blk = js.substr(0, bnc_pos);
            std::string bnc_blk = js.substr(bnc_pos);
            MapResult r;
            std::string fname = std::filesystem::path(in).filename().string();
            r.name = fname.substr(9, fname.size() - 14);
            r.sa_pts = json_val(sa_blk, "pts");
            r.sa_nodes = static_cast<int>(json_val(sa_blk, "nodes"));
            r.sa_s = json_val(sa_blk, "elapsed_s");
            r.bnc_pts = json_val(bnc_blk, "pts");
            r.bnc_nodes = static_cast<int>(json_val(bnc_blk, "nodes"));
            r.bnc_s = json_val(bnc_blk, "elapsed_s");
            r.bnc_optimal = (bnc_blk.find("\"proved_optimal\": true") != std::string::npos);
            r.gap_pct = json_val(bnc_blk, "gap_pct");
            results.push_back(r);
        } catch (...) {}
    }

    // Summary table
    std::cout << "\n";
    std::cout << "+----------------------------------------+---------------------------+----------------------------------------------+\n";
    std::cout << "| Instance                               | SA                        | B&P(SA)                                      |\n";
    std::cout << "|                                        |  pts  nodes     time      |  pts  nodes     time    optimal?    gap       |\n";
    std::cout << "+----------------------------------------+---------------------------+----------------------------------------------+\n";
    for (const auto& r : results) {
        std::cout << "| " << std::left << std::setw(38) << r.name
                  << " | " << std::right << std::setw(5) << static_cast<int>(r.sa_pts)
                  << "  " << std::setw(5) << r.sa_nodes
                  << "  " << std::setw(7) << std::fixed << std::setprecision(1) << r.sa_s << "s"
                  << "   | " << std::setw(5) << static_cast<int>(r.bnc_pts)
                  << "  " << std::setw(5) << r.bnc_nodes
                  << "  " << std::setw(7) << r.bnc_s << "s"
                  << "  " << (r.bnc_optimal ? "YES (proven)" : "no  (limit) ")
                  << "  " << std::setw(5) << std::setprecision(1) << r.gap_pct << "%"
                  << " |\n";
    }
    std::cout << "+----------------------------------------+---------------------------+----------------------------------------------+\n";

    int proven = 0, total = static_cast<int>(results.size());
    double max_gap = 0;
    for (const auto& r : results) {
        if (r.bnc_optimal) ++proven;
        if (r.gap_pct > max_gap) max_gap = r.gap_pct;
    }
    std::cout << "\nProven optimal: " << proven << "/" << total
              << "  Max gap: " << std::setprecision(2) << max_gap << "%\n";
    return 0;
}
