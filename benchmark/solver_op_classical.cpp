// Classical exact Orienteering Problem solver (undirected, symmetric costs,
// no fatigue), following Fischetti, Salazar-Gonzalez & Toth (1998),
// "Solving the Orienteering Problem through Branch-and-Cut" -- the same
// paper this project's own directed/fatigue-extended cuts (SECs, cover
// cuts, path/cycle inequalities) are adapted from.
//
// Purpose: a baseline for the "no comparison against other exact OP
// methods" review feedback. Deliberately NOT the AOPF (no asymmetry, no
// fatigue, no state/flow variables) -- for an exact apples-to-apples check
// against this project's flow solver (both solving mathematically the same
// problem), only instances with BOTH symmetric costs AND fatigue_rate=0
// qualify: currently that's bench_061_symmetric_easy and
// bench_066_symmetric_easy_large. The a00 sweep points (bench_028-030) have
// asymmetry=0 but fatigue_rate=0.2 -- this solver still runs on them, but
// since it ignores fatigue entirely, its optimum is only a valid UPPER
// BOUND on the flow solver's fatigue-constrained optimum there, not an
// exact match.
//
// Model:
//   x_e in {0,1} for each undirected edge e={i,j}, i<j, with finite symmetric cost.
//   y_i in {0,1} for each node (y_0 fixed to 1).
//   max sum p_i y_i
//   s.t. sum_{e incident to i} x_e = 2 y_i           (degree)
//        sum_e c_e x_e <= B                          (budget)
//        sum_{e within S} x_e <= sum_{i in S} y_i - 1 (SEC, separated dynamically)
//
// Symmetric edge cost c_e = (cm[i][j] + cm[j][i]) / 2, since this project's
// underlying instances are directed/asymmetric by construction (terrain
// slope) -- averaging is the standard way to get a fair symmetric baseline
// from the same instance data rather than requiring separately-generated
// symmetric instances.

#include <iostream>
#include <fstream>
#include <chrono>
#include <limits>
#include <queue>
#include <stack>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <random>
#include <set>
#include <filesystem>
#include "interfaces/highs_c_api.h"
#include "lp_data/HConst.h"
#include "Highs.h"

// ── JSON parser (same format as the other solvers; only cm/pts/bud_raw used) ──

static std::string read_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot open " + path);
    return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

struct Input {
    std::vector<std::vector<double>> cm;
    std::vector<double> pts;
    double bud_raw = 0.0;
};

static Input parse_input(const std::string& json_str) {
    Input inp;
    size_t i = 0;
    auto skip_ws = [&](size_t pos) -> size_t {
        while (pos < json_str.size() && (json_str[pos] == ' ' || json_str[pos] == '\t' ||
               json_str[pos] == '\n' || json_str[pos] == '\r')) ++pos;
        return pos;
    };
    auto parse_number = [&json_str](size_t& pos) -> double {
        size_t start = pos;
        if (pos < json_str.size() && (json_str[pos] == '-' || json_str[pos] == '+')) ++pos;
        while (pos < json_str.size() && (std::isdigit(json_str[pos]) || json_str[pos] == '.' ||
               json_str[pos] == 'e' || json_str[pos] == 'E' || json_str[pos] == '+' || json_str[pos] == '-')) ++pos;
        return std::stod(json_str.substr(start, pos - start));
    };
    auto parse_array1d = [&](size_t& pos) -> std::vector<double> {
        std::vector<double> v;
        pos = skip_ws(pos);
        if (pos >= json_str.size() || json_str[pos] != '[') throw std::runtime_error("Expected [");
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
        if (pos >= json_str.size() || json_str[pos] != '[') throw std::runtime_error("Expected [");
        ++pos;
        while (true) {
            pos = skip_ws(pos);
            if (json_str[pos] == ']') { ++pos; break; }
            if (json_str[pos] == ',') { ++pos; continue; }
            if (json_str[pos] == '[') m.push_back(parse_array1d(pos));
        }
        return m;
    };
    auto find_key_opt = [&](const std::string& key) -> bool {
        std::string key_str = "\"" + key + "\":";
        size_t pos = json_str.find(key_str, 0);
        if (pos == std::string::npos) return false;
        i = pos + key_str.size();
        return true;
    };
    auto find_key = [&](const std::string& key) {
        if (!find_key_opt(key)) throw std::runtime_error("Missing key: " + key);
    };
    find_key("cm"); inp.cm = parse_array2d(i);
    find_key("pts"); inp.pts = parse_array1d(i);
    find_key("bud_raw"); i = skip_ws(i); inp.bud_raw = parse_number(i);
    return inp;
}

// Symmetric edge cost, averaged from the (possibly asymmetric) directed
// matrix -- see file header. inf if either direction is inf.
static std::vector<std::vector<double>> symmetrize(const std::vector<std::vector<double>>& cm) {
    int n = static_cast<int>(cm.size());
    std::vector<std::vector<double>> c(n, std::vector<double>(n, std::numeric_limits<double>::infinity()));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            if (std::isfinite(cm[i][j]) && std::isfinite(cm[j][i]))
                c[i][j] = c[j][i] = 0.5 * (cm[i][j] + cm[j][i]);
        }
    return c;
}

// ── LP model ────────────────────────────────────────────────────────────

struct Edge { int i, j; };  // i < j

struct LPModel {
    int n = 0;
    Highs* highs = nullptr;
    std::vector<Edge> edges;                  // edge index -> (i,j), i<j
    std::vector<std::vector<int>> edge_of;     // [i][j] (either order) -> edge index, or -1
    std::vector<int> x_col;                    // edge index -> LP column
    std::vector<int> y_col;                     // node -> LP column
    std::vector<double> col_ub_cache;
    std::vector<double> sol_cache;
    int n_rows_base = 0;

    LPModel() = default;
    ~LPModel() { if (highs) Highs_destroy(highs); }
    LPModel(const LPModel&) = delete;
    LPModel& operator=(const LPModel&) = delete;

    int add_col(double lb, double ub, double obj = 0.0) {
        Highs_addCol(highs, obj, lb, ub, 0, nullptr, nullptr);
        col_ub_cache.push_back(ub);
        return static_cast<int>(col_ub_cache.size()) - 1;
    }
    void add_row(double lhs, double rhs, const std::vector<int>& cols, const std::vector<double>& coeffs) {
        assert(cols.size() == coeffs.size());
        Highs_addRow(highs, lhs, rhs, static_cast<int>(cols.size()), cols.data(), coeffs.data());
    }
    void fix_col(int col, double val) {
        Highs_changeColBounds(highs, col, val, val);
        col_ub_cache[col] = val;
    }
    double get_col_ub(int col) const { return col_ub_cache[col]; }

    void build(const Input& inp, const std::vector<std::vector<double>>& c) {
        n = static_cast<int>(inp.pts.size());
        highs = (Highs*)Highs_create();
        Highs_setBoolOptionValue(highs, "output_flag", false);
        Highs_setStringOptionValue(highs, "presolve", "on");
        Highs_setStringOptionValue(highs, "solver", "simplex");

        edge_of.assign(n, std::vector<int>(n, -1));
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) {
                if (!std::isfinite(c[i][j])) continue;
                if (c[i][j] > inp.bud_raw) continue;
                // Arc elimination: even the cheapest possible route using
                // this edge must fit the budget. When neither endpoint is
                // the depot, that's the direct triangle 0-i-j-0 (c is
                // symmetric, so c[j][0]=c[0][j]). When i==0, edge (0,j) IS
                // the depot leg already -- there is no third leg, and using
                // c[0][i] here would read c[0][0], the diagonal, which is
                // left at +infinity (never set, since i==j is always
                // skipped) -- that silently excluded every depot edge,
                // making the whole LP infeasible from the root.
                if (i == 0) {
                    if (2.0 * c[0][j] > inp.bud_raw) continue;
                } else {
                    if (!std::isfinite(c[0][i]) || !std::isfinite(c[0][j])) continue;
                    if (c[0][i] + c[i][j] + c[0][j] > inp.bud_raw) continue;
                }
                int eidx = static_cast<int>(edges.size());
                edges.push_back({i, j});
                edge_of[i][j] = edge_of[j][i] = eidx;
            }

        x_col.assign(edges.size(), -1);
        for (size_t e = 0; e < edges.size(); ++e)
            x_col[e] = add_col(0.0, 1.0);

        y_col.assign(n, -1);
        for (int i = 0; i < n; ++i)
            y_col[i] = add_col(0.0, 1.0, inp.pts[i]);
        fix_col(y_col[0], 1.0);

        // Degree: sum of incident edges = 2*y_i
        for (int i = 0; i < n; ++i) {
            std::vector<int> cols; std::vector<double> coeffs;
            for (int j = 0; j < n; ++j) {
                if (edge_of[i][j] < 0) continue;
                cols.push_back(x_col[edge_of[i][j]]);
                coeffs.push_back(1.0);
            }
            cols.push_back(y_col[i]);
            coeffs.push_back(-2.0);
            add_row(0.0, 0.0, cols, coeffs);
        }

        // Budget: sum c_e x_e <= B
        {
            std::vector<int> cols; std::vector<double> coeffs;
            for (size_t e = 0; e < edges.size(); ++e) {
                cols.push_back(x_col[e]);
                coeffs.push_back(c[edges[e].i][edges[e].j]);
            }
            add_row(-1e30, inp.bud_raw, cols, coeffs);
        }

        n_rows_base = Highs_getNumRow(highs);
    }

    void clone_from(const LPModel& other, const std::vector<std::pair<int,double>>& fixings) {
        assert(highs == nullptr);
        n = other.n;
        edges = other.edges;
        edge_of = other.edge_of;
        x_col = other.x_col;
        y_col = other.y_col;
        col_ub_cache = other.col_ub_cache;
        n_rows_base = other.n_rows_base;

        highs = (Highs*)Highs_create();
        Highs_setBoolOptionValue(highs, "output_flag", false);
        Highs_setStringOptionValue(highs, "presolve", "on");
        Highs_setStringOptionValue(highs, "solver", "simplex");

        int nc = Highs_getNumCol(other.highs);
        int nr = Highs_getNumRow(other.highs);
        int nnz = Highs_getNumNz(other.highs);
        std::vector<double> costs(nc), lb(nc), ub(nc), rlb(nr), rub(nr);
        std::vector<int> astart(nc), aindex(nnz);
        std::vector<double> avalue(nnz);
        HighsInt sense, num_col, num_row, num_nz;
        double offset;
        std::vector<HighsInt> integrality(nc);
        Highs_getLp(other.highs, kHighsMatrixFormatColwise, &num_col, &num_row, &num_nz, &sense, &offset,
                    costs.data(), lb.data(), ub.data(), rlb.data(), rub.data(),
                    astart.data(), aindex.data(), avalue.data(), integrality.data());
        Highs_passLp(highs, nc, nr, nnz, kHighsMatrixFormatColwise, sense, offset,
                     costs.data(), lb.data(), ub.data(), rlb.data(), rub.data(),
                     astart.data(), aindex.data(), avalue.data());
        for (const auto& [col, val] : fixings) fix_col(col, val);
    }

    void add_sec(const std::vector<int>& S) {
        std::set<int> Sset(S.begin(), S.end());
        std::vector<int> cols; std::vector<double> coeffs;
        for (int i : S)
            for (int j : S) {
                if (i >= j) continue;
                if (edge_of[i][j] < 0) continue;
                cols.push_back(x_col[edge_of[i][j]]);
                coeffs.push_back(1.0);
            }
        for (int i : S) { cols.push_back(y_col[i]); coeffs.push_back(-1.0); }
        if (!cols.empty()) add_row(-1e30, -1.0, cols, coeffs);
    }

    bool solve() {
        Highs_changeObjectiveSense(highs, -1);  // maximize
        int status = Highs_run(highs);
        bool ok = status == 0 && Highs_getModelStatus(highs) == 7;
        if (ok) {
            int nc = Highs_getNumCol(highs);
            int nr = Highs_getNumRow(highs);
            sol_cache.resize(nc);
            std::vector<double> col_dual(nc), row_val(nr), row_dual(nr);
            Highs_getSolution(highs, sol_cache.data(), col_dual.data(), row_val.data(), row_dual.data());
        }
        return ok;
    }
    double obj() const { return Highs_getObjectiveValue(highs); }
    double prim(int col) const { return sol_cache[col]; }
};

// ── Connectivity separation (undirected: connected components via Union-Find) ──

struct DSU {
    std::vector<int> parent;
    explicit DSU(int n) : parent(n) { std::iota(parent.begin(), parent.end(), 0); }
    int find(int x) { return parent[x] == x ? x : parent[x] = find(parent[x]); }
    void unite(int a, int b) { a = find(a); b = find(b); if (a != b) parent[a] = b; }
};

// Returns node sets S (not containing the depot) that are violated
// subtours in the current (possibly fractional) support graph: connected
// components, formed from edges with x* > eps, that don't reach node 0.
static std::vector<std::vector<int>> find_subtours(const LPModel& model, double eps = 0.5) {
    DSU dsu(model.n);
    for (size_t e = 0; e < model.edges.size(); ++e)
        if (model.prim(model.x_col[e]) > eps)
            dsu.unite(model.edges[e].i, model.edges[e].j);

    std::vector<std::vector<int>> by_root(model.n);
    for (int i = 0; i < model.n; ++i)
        if (model.prim(model.y_col[i]) > eps)
            by_root[dsu.find(i)].push_back(i);

    std::vector<std::vector<int>> subtours;
    int depot_root = dsu.find(0);
    for (int r = 0; r < model.n; ++r)
        if (r != depot_root && !by_root[r].empty())
            subtours.push_back(by_root[r]);
    return subtours;
}

int find_and_add_sec_cuts(LPModel& lp, int max_cuts = 5) {
    auto subtours = find_subtours(lp);
    int added = 0;
    for (const auto& S : subtours) {
        if (added >= max_cuts) break;
        lp.add_sec(S);
        ++added;
    }
    return added;
}

// ── Route extraction / feasibility ─────────────────────────────────────

// Traces the depot's 2-cycle from x_e > eps edges. Assumes (post
// integrality-check) a genuine simple cycle through the depot.
static std::vector<int> extract_route(const LPModel& model, double eps = 0.5) {
    int n = model.n;
    std::vector<std::vector<int>> adj(n);
    for (size_t e = 0; e < model.edges.size(); ++e)
        if (model.prim(model.x_col[e]) > eps) {
            adj[model.edges[e].i].push_back(model.edges[e].j);
            adj[model.edges[e].j].push_back(model.edges[e].i);
        }
    std::vector<int> route;
    std::vector<bool> visited(n, false);
    int prev = -1, cur = 0;
    visited[0] = true;
    while (true) {
        int next = -1;
        for (int cand : adj[cur])
            if (cand != prev && (!visited[cand] || cand == 0)) { next = cand; break; }
        if (next == -1 || next == 0) break;
        route.push_back(next);
        visited[next] = true;
        prev = cur; cur = next;
        if (static_cast<int>(route.size()) >= n) break;
    }
    return route;
}

static double rcost(const std::vector<std::vector<double>>& c, const std::vector<int>& route) {
    if (route.empty()) return 0.0;
    double cost = c[0][route[0]];
    for (size_t i = 0; i + 1 < route.size(); ++i) cost += c[route[i]][route[i + 1]];
    return cost + c[route.back()][0];
}

static bool is_feasible_route(const std::vector<std::vector<double>>& c, const std::vector<int>& route, double bud_raw) {
    return rcost(c, route) <= bud_raw + 1e-6;
}

// ── Greedy warm start (simple nearest-fit; this is a reference baseline,
// not the sophisticated SA used by the project's own solvers) ──────────

static std::vector<int> greedy_route(const Input& inp, const std::vector<std::vector<double>>& c) {
    int n = static_cast<int>(inp.pts.size());
    std::vector<bool> visited(n, false);
    visited[0] = true;
    std::vector<int> route;
    double elapsed = 0.0;
    int cur = 0;
    while (true) {
        int best = -1; double best_score = -1e30;
        for (int j = 1; j < n; ++j) {
            if (visited[j] || !std::isfinite(c[cur][j])) continue;
            double added = c[cur][j] + c[j][0] - c[cur][0];
            if (elapsed + c[cur][j] + c[j][0] > inp.bud_raw) continue;
            double score = inp.pts[j] / std::max(added, 1e-6);
            if (score > best_score) { best_score = score; best = j; }
        }
        if (best == -1) break;
        elapsed += c[cur][best];
        cur = best; visited[cur] = true;
        route.push_back(cur);
    }
    return route;
}

// ── Branch-and-Cut ─────────────────────────────────────────────────────

struct BNCNode {
    std::vector<std::pair<int,double>> fixings;
    double ub = std::numeric_limits<double>::infinity();
    bool operator<(const BNCNode& o) const { return ub < o.ub; }  // max-heap by ub
};

struct Solver {
    const Input& inp;
    const std::vector<std::vector<double>>& c;
    LPModel root;
    double best_pts = 0.0;
    std::vector<int> best_route;
    int max_cuts = 20;
    int max_depth = 200;
    double time_limit_s = 900.0;
    bool depth_limit_hit = false;
    bool proved_optimal = false;
    double best_ub = std::numeric_limits<double>::infinity();
    int nodes_explored = 0;

    Solver(const Input& i, const std::vector<std::vector<double>>& cc) : inp(i), c(cc) {
        root.build(inp, c);
    }

    void solve(double warm_start_pts, std::vector<int> warm_start_route) {
        if (warm_start_pts > best_pts) { best_pts = warm_start_pts; best_route = std::move(warm_start_route); }
        auto gr = greedy_route(inp, c);
        double gr_pts = std::accumulate(gr.begin(), gr.end(), 0.0, [&](double s, int v){ return s + inp.pts[v]; });
        if (gr_pts > best_pts && is_feasible_route(c, gr, inp.bud_raw)) { best_pts = gr_pts; best_route = std::move(gr); }
        std::cerr << "B&C warm start: " << best_pts << " pts\n";

        std::priority_queue<BNCNode> node_pq;
        BNCNode root_node;
        bool root_ok = root.solve();
        root_node.ub = root_ok ? root.obj() : -std::numeric_limits<double>::infinity();
        if (root_node.ub > best_pts) node_pq.push(std::move(root_node));

        int nodes = 0;
        auto t_start = std::chrono::steady_clock::now();
        while (!node_pq.empty() && nodes++ < 10000) {
            double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
            if (elapsed > time_limit_s) { std::cerr << "Time limit reached\n"; break; }
            BNCNode node = std::move(const_cast<BNCNode&>(node_pq.top()));
            node_pq.pop();
            if (node.ub <= best_pts + 1e-6) continue;
            process_node(std::move(node), node_pq);
        }
        proved_optimal = node_pq.empty() && nodes < 10000 && !depth_limit_hit;
        nodes_explored = nodes;
        if (proved_optimal) {
            best_ub = best_pts;
        } else {
            best_ub = best_pts;
            auto tmp = node_pq;
            while (!tmp.empty()) { if (tmp.top().ub > best_ub) best_ub = tmp.top().ub; tmp.pop(); }
        }
        std::cerr << "Processed " << nodes << " nodes, best: " << best_pts << " pts, UB: " << best_ub << "\n";
    }

    void process_node(BNCNode node, std::priority_queue<BNCNode>& node_pq) {
        if (node.ub <= best_pts + 1e-6) return;
        if (node.fixings.size() > static_cast<size_t>(max_depth)) { depth_limit_hit = true; return; }

        LPModel lp;
        lp.clone_from(root, node.fixings);

        for (int iter = 0; iter < max_cuts; ++iter) {
            if (!lp.solve()) return;
            double lp_ub = lp.obj();
            if (lp_ub <= best_pts + 1e-6) return;
            int added = find_and_add_sec_cuts(lp);
            if (added == 0) break;
        }
        if (!lp.solve()) return;
        double lp_ub = lp.obj();
        if (lp_ub <= best_pts + 1e-6) return;

        bool integer_sol = true;
        int branch_col = -1;
        double max_frac = 0.0;
        for (int i = 1; i < root.n; ++i) {
            double v = lp.prim(root.y_col[i]);
            double frac = std::min(v, 1.0 - v);
            if (frac > 1e-5) { integer_sol = false; if (frac > max_frac) { max_frac = frac; branch_col = root.y_col[i]; } }
        }
        // y integral alone doesn't certify an integer route -- same lesson
        // as the project's other solvers: check x_e integrality too before
        // accepting a leaf.
        if (integer_sol) {
            for (size_t e = 0; e < root.edges.size() && integer_sol; ++e) {
                double v = lp.prim(root.x_col[e]);
                double frac = std::min(v, 1.0 - v);
                if (frac > 1e-5) { integer_sol = false; if (frac > max_frac) { max_frac = frac; branch_col = root.x_col[e]; } }
            }
        }

        if (integer_sol) {
            auto route = extract_route(lp);
            double pts = std::accumulate(route.begin(), route.end(), 0.0, [&](double s, int v){ return s + inp.pts[v]; });
            if (is_feasible_route(c, route, inp.bud_raw) && pts > best_pts + 1e-6) {
                best_pts = pts; best_route = std::move(route);
                std::cerr << "New best: " << best_pts << " pts\n";
            }
        } else if (branch_col >= 0) {
            BNCNode node0 = node, node1 = node;
            node0.fixings.emplace_back(branch_col, 0.0); node0.ub = lp_ub;
            node1.fixings.emplace_back(branch_col, 1.0); node1.ub = lp_ub;
            node_pq.push(std::move(node0));
            node_pq.push(std::move(node1));
        }
    }
};

// ── Main ────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {
    std::string input_dir = "instances";
    for (int a = 1; a < argc; ++a) input_dir = argv[a];

    std::vector<std::pair<std::string,std::string>> maps;
    for (const auto& entry : std::filesystem::directory_iterator(input_dir)) {
        std::string fname = entry.path().filename().string();
        if (fname.substr(0, 9) == "op_input_" && fname.size() > 14 && fname.substr(fname.size() - 5) == ".json") {
            std::string base = fname.substr(9, fname.size() - 14);
            std::string in_path = entry.path().string();
            std::string out_path = (entry.path().parent_path() / ("op_output_classical_" + base + ".json")).string();
            maps.push_back({in_path, out_path});
        }
    }
    std::sort(maps.begin(), maps.end());
    std::cerr << "Found " << maps.size() << " instances in " << input_dir << "/\n";

    for (const auto& [in_path, out_path] : maps) {
        std::cerr << "\n=== " << in_path << " ===\n";
        Input inp = parse_input(read_file(in_path));
        auto c = symmetrize(inp.cm);

        // Cheap symmetric-cost greedy for the warm start pts/route
        auto t0 = std::chrono::steady_clock::now();
        auto gr_route = greedy_route(inp, c);
        double gr_pts = std::accumulate(gr_route.begin(), gr_route.end(), 0.0,
            [&](double s, int v){ return s + inp.pts[v]; });
        double gr_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::cerr << "Greedy warm start: " << gr_pts << " pts (" << gr_route.size() << " nodes) in " << gr_elapsed << "s\n";

        auto t1 = std::chrono::steady_clock::now();
        Solver solver(inp, c);
        solver.solve(gr_pts, gr_route);
        double bnc_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
        std::cerr << "B&C: " << solver.best_pts << " pts (" << solver.best_route.size() << " nodes) in " << bnc_elapsed << "s\n";

        std::ofstream out(out_path);
        out << "{\n";
        out << "  \"model\": \"classical_op_fischetti1998\",\n";
        out << "  \"greedy\": {\"pts\": " << gr_pts << ", \"nodes\": " << gr_route.size()
            << ", \"elapsed_s\": " << gr_elapsed << ", \"base_cost\": " << rcost(c, gr_route) << ", \"route\": [";
        for (size_t i = 0; i < gr_route.size(); ++i) { if (i) out << ", "; out << gr_route[i]; }
        out << "]},\n";
        out << "  \"bnc\": {\"pts\": " << solver.best_pts << ", \"nodes\": " << solver.best_route.size()
            << ", \"elapsed_s\": " << bnc_elapsed
            << ", \"proved_optimal\": " << (solver.proved_optimal ? "true" : "false")
            << ", \"best_ub\": " << solver.best_ub
            << ", \"gap_pct\": " << (solver.best_pts > 0 ? 100.0 * (solver.best_ub - solver.best_pts) / solver.best_pts : 0.0)
            << ", \"nodes_explored\": " << solver.nodes_explored
            << ", \"base_cost\": " << rcost(c, solver.best_route) << ", \"route\": [";
        for (size_t i = 0; i < solver.best_route.size(); ++i) { if (i) out << ", "; out << solver.best_route[i]; }
        out << "]}\n";
        out << "}\n";
    }
    return 0;
}
