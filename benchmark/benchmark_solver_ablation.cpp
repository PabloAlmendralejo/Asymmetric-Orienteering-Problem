
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
#include "interfaces/highs_c_api.h"
#include "lp_data/HConst.h"
//  #include "glpk.h"

// ── Global ablation flags ──────────────────────────────────────────────────
static bool g_use_tightened_coupling = true;
static bool g_use_fatigue_covers = true;

// ── JSON parser (improved from original) ───────────────────────────────────

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
    // Simple but robust JSON parser for this exact format
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
    for (size_t i = 0; i + 1 < route.size(); ++i) c += cm[route[i]][route[i + 1]];
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
        double leg = cm[seq[i]][seq[i + 1]];
        total += leg * (1.0 + fatigue_rate * (elapsed / std::max(bud_raw, 1.0)));
        elapsed += leg;
    }
    return total;
}

#include "Highs.h"  // replace #include "glpk.h"

struct LPModel {
    int n = 0;
    Highs* highs = nullptr;
    std::vector<std::vector<int>> x_col;   // arc variables
    std::vector<int> y_col;                // node visit variables
    std::vector<std::vector<int>> f_col;   // flow variables (cumulative time)
    std::vector<double> col_ub_cache;
    std::vector<double> sol_cache;
    int n_cols_base = 0;
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

    void add_row(double lhs, double rhs,
                 const std::vector<int>& cols,
                 const std::vector<double>& coeffs) {
        assert(cols.size() == coeffs.size());
        Highs_addRow(highs, lhs, rhs,
                     static_cast<int>(cols.size()),
                     cols.data(), coeffs.data());
    }

    void fix_col(int col, double val) {
        Highs_changeColBounds(highs, col, val, val);
        col_ub_cache[col] = val;
    }

    double get_col_ub(int col) const {
        return col_ub_cache[col];
    }

    // ── build (flow-based formulation) ─────────────────────────────────────
    //
    // Variables:
    //   x_ij ∈ {0,1}  — arc traversal
    //   y_i  ∈ {0,1}  — node visit
    //   f_ij ≥ 0      — cumulative elapsed time flowing along arc (i,j)
    //
    // The flow f_ij carries the arrival time at node i when arc (i,j) is used:
    //   f_ij = t_i · x_ij  (exactly, not an approximation)
    //
    // Flow conservation gives arrival times:
    //   t_j = ∑_i f_ij + C_ij · x_ij  (for visited j)
    //   equivalently: ∑_i (f_ij + C_ij · x_ij) = ∑_k f_jk  (time in = time out)
    //
    // Fatigue budget (LINEAR — no McCormick needed):
    //   ∑ C_ij · x_ij + (λ/B) · ∑ C_ij · f_ij ≤ B
    //
    // This replaces MTZ + McCormick with flow constraints that give a
    // tighter LP relaxation.

    void build(const Input& inp) {
        n = static_cast<int>(inp.pts.size());
        highs = (Highs*)Highs_create();
        Highs_setBoolOptionValue(highs, "output_flag", false);
        Highs_setStringOptionValue(highs, "presolve", "on");
        Highs_setStringOptionValue(highs, "solver", "simplex");
        Highs_setStringOptionValue(highs, "simplex_strategy", "1");

        x_col.assign(n, std::vector<int>(n, -1));
        y_col.resize(n, -1);
        f_col.assign(n, std::vector<int>(n, -1));

        // x[i][j] — arc variables with structural elimination
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (i == j || !std::isfinite(inp.cm[i][j])) continue;
                bool infeasible = !std::isfinite(inp.cm[j][0]) ||
                                  inp.cm[i][j] + inp.cm[j][0] > inp.bud_raw ||
                                  inp.cm[0][i] + inp.cm[i][j] + inp.cm[j][0] > inp.bud_raw;
                if (!infeasible && inp.fatigue_rate > 0) {
                    double t_i = inp.cm[0][i];
                    double t_j = t_i + inp.cm[i][j];
                    double fat_cost = inp.cm[0][i]
                                    + inp.cm[i][j] * (1.0 + inp.fatigue_rate * t_i / inp.bud_raw)
                                    + inp.cm[j][0] * (1.0 + inp.fatigue_rate * t_j / inp.bud_raw);
                    if (fat_cost > inp.bud_raw) infeasible = true;
                }
                x_col[i][j] = add_col(0.0, infeasible ? 0.0 : 1.0);
            }

        // y[i] — node visit variables (objective coefficients = scores)
        for (int i = 0; i < n; ++i)
            y_col[i] = add_col(0.0, 1.0, inp.pts[i]);
        fix_col(y_col[0], 1.0);

        // f[i][j] — flow variables: f_ij carries arrival time at i along arc (i,j)
        //   0 ≤ f_ij ≤ B · x_ij (enforced via coupling constraints below)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (x_col[i][j] < 0) continue;
                f_col[i][j] = add_col(0.0, inp.bud_raw);
            }

        add_flow_conservation();
        add_time_flow_coupling(inp.cm, inp.bud_raw);
        add_time_flow_propagation(inp.cm, inp.bud_raw);
        add_fatigue_budget_flow(inp.cm, inp.bud_raw, inp.fatigue_rate);

        n_rows_base = Highs_getNumRow(highs);
    }

    void clone_from(const LPModel& other,
                    const std::vector<std::pair<int,double>>& fixings) {
        assert(highs == nullptr);
        n               = other.n;
        x_col           = other.x_col;
        y_col           = other.y_col;
        f_col           = other.f_col;
        col_ub_cache    = other.col_ub_cache;
        n_rows_base     = other.n_rows_base;
        // added_secs starts empty — each cloned node tracks its own cuts

        // Deep-copy the HiGHS model
        highs = (Highs*)Highs_create();
        Highs_setBoolOptionValue(highs, "output_flag", false);
        Highs_setStringOptionValue(highs, "presolve", "on");
        Highs_setStringOptionValue(highs, "solver", "simplex");

        // Deep-copy via Highs_passLp — no temp file, faster, thread-safe
        {
            int nc = Highs_getNumCol(other.highs);
            int nr = Highs_getNumRow(other.highs);
            int nnz = Highs_getNumNz(other.highs);

            std::vector<double> costs(nc), lb(nc), ub(nc), rlb(nr), rub(nr);
            std::vector<int> astart(nc), aindex(nnz);
            std::vector<double> avalue(nnz);
            HighsInt sense, num_col, num_row, num_nz;
            double offset;
            std::vector<HighsInt> integrality(nc);

            Highs_getLp(other.highs, kHighsMatrixFormatColwise,
                        &num_col, &num_row, &num_nz, &sense, &offset,
                        costs.data(), lb.data(), ub.data(),
                        rlb.data(), rub.data(),
                        astart.data(), aindex.data(), avalue.data(),
                        integrality.data());

            Highs_passLp(highs, nc, nr, nnz,
                         kHighsMatrixFormatColwise, sense, offset,
                         costs.data(), lb.data(), ub.data(),
                         rlb.data(), rub.data(),
                         astart.data(), aindex.data(), avalue.data());
        }

        for (const auto& [col, val] : fixings)
            fix_col(col, val);
    }

    // ── constraints (flow-based) ──────────────────────────────────────────

    void add_flow_conservation() {
        // Arc in-flow: sum_j x[j][i] = y[i]
        for (int i = 0; i < n; ++i) {
            std::vector<int> cols; std::vector<double> coeffs;
            for (int j = 0; j < n; ++j)
                if (x_col[j][i] >= 0) { cols.push_back(x_col[j][i]); coeffs.push_back(1.0); }
            cols.push_back(y_col[i]); coeffs.push_back(-1.0);
            add_row(0.0, 0.0, cols, coeffs);
        }
        // Arc out-flow: sum_j x[i][j] = y[i]
        for (int i = 0; i < n; ++i) {
            std::vector<int> cols; std::vector<double> coeffs;
            for (int j = 0; j < n; ++j)
                if (x_col[i][j] >= 0) { cols.push_back(x_col[i][j]); coeffs.push_back(1.0); }
            cols.push_back(y_col[i]); coeffs.push_back(-1.0);
            add_row(0.0, 0.0, cols, coeffs);
        }
    }

    void add_time_flow_coupling(const std::vector<std::vector<double>>& cm, double bud_raw) {
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (f_col[i][j] < 0) continue;

                double ub_coeff = bud_raw;
                if (g_use_tightened_coupling) {
                    // Tightened: f[i][j] <= (B - C_ij - C_j0) * x[i][j]
                    if (std::isfinite(cm[i][j]) && j < (int)cm.size()
                        && std::isfinite(cm[j][0])) {
                        ub_coeff = std::max(0.0, bud_raw - cm[i][j] - cm[j][0]);
                    }
                }
                // else: original loose coupling f[i][j] <= B * x[i][j]
                add_row(-1e30, 0.0,
                        {f_col[i][j], x_col[i][j]},
                        {1.0,         -ub_coeff});

                if (g_use_tightened_coupling) {
                    // Lower bound: f[i][j] >= C_0i * x[i][j]
                    if (i > 0 && std::isfinite(cm[0][i]) && cm[0][i] > 1e-9) {
                        add_row(0.0, 1e30,
                                {f_col[i][j], x_col[i][j]},
                                {1.0,         -cm[0][i]});
                    }
                }
            }
    }

    void add_time_flow_propagation(const std::vector<std::vector<double>>& cm, double bud_raw) {
        // At each non-depot node j:
        //   sum_i (f[i][j] + C[i][j] * x[i][j]) = sum_k f[j][k]
        for (int j = 1; j < n; ++j) {
            std::vector<int> cols; std::vector<double> coeffs;
            for (int i = 0; i < n; ++i) {
                if (f_col[i][j] >= 0) { cols.push_back(f_col[i][j]); coeffs.push_back(1.0); }
                if (x_col[i][j] >= 0) { cols.push_back(x_col[i][j]); coeffs.push_back(cm[i][j]); }
            }
            for (int k = 0; k < n; ++k) {
                if (f_col[j][k] >= 0) { cols.push_back(f_col[j][k]); coeffs.push_back(-1.0); }
            }
            if (!cols.empty())
                add_row(0.0, 0.0, cols, coeffs);
        }
        // Depot outgoing flow = 0
        for (int j = 0; j < n; ++j) {
            if (f_col[0][j] >= 0)
                fix_col(f_col[0][j], 0.0);
        }
    }

    void add_fatigue_budget_flow(const std::vector<std::vector<double>>& cm,
                                 double bud_raw, double fatigue_rate) {
        // sum C_ij * x_ij + (lambda/B) * sum C_ij * f_ij <= B
        std::vector<int> cols; std::vector<double> coeffs;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (x_col[i][j] >= 0) {
                    cols.push_back(x_col[i][j]);
                    coeffs.push_back(cm[i][j]);
                }
                if (f_col[i][j] >= 0) {
                    cols.push_back(f_col[i][j]);
                    coeffs.push_back((fatigue_rate / bud_raw) * cm[i][j]);
                }
            }
        add_row(-1e30, bud_raw, cols, coeffs);
    }

    std::set<std::vector<int>> added_secs;

    void add_sec(const std::vector<int>& S) {
        std::vector<int> key(S.begin(), S.end());
        std::sort(key.begin(), key.end());
        if (!added_secs.insert(key).second) return;  // already added
        std::vector<int> Sset(S.begin(), S.end());
        std::sort(Sset.begin(), Sset.end());
        auto in_S = [&](int v) {
            return std::binary_search(Sset.begin(), Sset.end(), v);
        };

        // 1. Directed SEC: sum_{i,j in S} x[i][j] <= |S|-1
        {
            std::vector<int> cols; std::vector<double> coeffs;
            for (int ii : S)
                for (int jj : S)
                    if (ii != jj && x_col[ii][jj] >= 0) {
                        cols.push_back(x_col[ii][jj]); coeffs.push_back(1.0);
                    }
            if (!cols.empty())
                add_row(-1e30, static_cast<double>(S.size()-1), cols, coeffs);
        }

        // Precompute outgoing/incoming arc lists
        std::vector<int> out_c, in_c;
        std::vector<double> out_v, in_v;
        for (int ii : S)
            for (int j = 0; j < n; ++j)
                if (!in_S(j) && x_col[ii][j] >= 0) {
                    out_c.push_back(x_col[ii][j]); out_v.push_back(1.0);
                }
        for (int i = 0; i < n; ++i)
            if (!in_S(i))
                for (int jj : S)
                    if (x_col[i][jj] >= 0) {
                        in_c.push_back(x_col[i][jj]); in_v.push_back(1.0);
                    }

        // 2+3. Outgoing + incoming cuts for each k in S
        for (int k : S) {
            if (!out_c.empty()) {
                auto cols = out_c; auto coeffs = out_v;
                cols.push_back(y_col[k]); coeffs.push_back(-1.0);
                add_row(0.0, 1e30, cols, coeffs);
            }
            if (!in_c.empty()) {
                auto cols = in_c; auto coeffs = in_v;
                cols.push_back(y_col[k]); coeffs.push_back(-1.0);
                add_row(0.0, 1e30, cols, coeffs);
            }
        }

        // 4. Combined: sum_out + sum_in >= 2*y[S[0]]
        if (!out_c.empty() || !in_c.empty()) {
            std::vector<int> cols; std::vector<double> coeffs;
            for (size_t t = 0; t < out_c.size(); ++t) {
                cols.push_back(out_c[t]); coeffs.push_back(1.0);
            }
            for (size_t t = 0; t < in_c.size(); ++t) {
                cols.push_back(in_c[t]); coeffs.push_back(1.0);
            }
            cols.push_back(y_col[S[0]]); coeffs.push_back(-2.0);
            add_row(0.0, 1e30, cols, coeffs);
        }
    }

    // ── solve ─────────────────────────────────────────────────────────────

    bool solve() {
        Highs_changeObjectiveSense(highs, -1);  // -1 = kHighsObjSenseMaximize
        int status = Highs_run(highs);
        bool ok = status == 0 && Highs_getModelStatus(highs) == 7;
        if (ok) {
            // Cache solution once so prim() is O(1)
            int nc = Highs_getNumCol(highs);
            int nr = Highs_getNumRow(highs);
            sol_cache.resize(nc);
            std::vector<double> col_dual(nc), row_val(nr), row_dual(nr);
            Highs_getSolution(highs, sol_cache.data(), col_dual.data(),
                              row_val.data(), row_dual.data());
        }
        return ok;
    }

    double obj() const {
        return Highs_getObjectiveValue(highs);
    }

    double prim(int col) const {
        return sol_cache[col];
    }

    void delete_extra_rows(int base_rows) {
        int cur = Highs_getNumRow(highs);
        if (cur <= base_rows) return;
        std::vector<int> to_del;
        for (int r = base_rows; r < cur; ++r) to_del.push_back(r);
        Highs_deleteRowsByRange(highs, base_rows, cur - 1);
    }
};


// ── Max-flow separation (push-relabel) ─────────────────────────────────────
// For each visited node k, compute max-flow from depot to k on the fractional
// graph. If max-flow < y_k, the min-cut gives a violated connectivity cut.

struct MaxFlow {
    struct Edge { int to, rev; double cap; };
    int n;
    std::vector<std::vector<Edge>> graph;
    std::vector<int> level, iter;

    MaxFlow(int n) : n(n), graph(n), level(n), iter(n) {}

    void add_edge(int from, int to, double cap) {
        graph[from].push_back({to, (int)graph[to].size(), cap});
        graph[to].push_back({from, (int)graph[from].size() - 1, 0.0});
    }

    bool bfs(int s, int t) {
        std::fill(level.begin(), level.end(), -1);
        std::queue<int> q;
        level[s] = 0;
        q.push(s);
        while (!q.empty()) {
            int v = q.front(); q.pop();
            for (auto& e : graph[v])
                if (e.cap > 1e-9 && level[e.to] < 0) {
                    level[e.to] = level[v] + 1;
                    q.push(e.to);
                }
        }
        return level[t] >= 0;
    }

    double dfs(int v, int t, double f) {
        if (v == t) return f;
        for (int& i = iter[v]; i < (int)graph[v].size(); ++i) {
            Edge& e = graph[v][i];
            if (e.cap > 1e-9 && level[v] < level[e.to]) {
                double d = dfs(e.to, t, std::min(f, e.cap));
                if (d > 1e-9) {
                    e.cap -= d;
                    graph[e.to][e.rev].cap += d;
                    return d;
                }
            }
        }
        return 0.0;
    }

    double max_flow(int s, int t) {
        double flow = 0.0;
        while (bfs(s, t)) {
            std::fill(iter.begin(), iter.end(), 0);
            double d;
            while ((d = dfs(s, t, 1e18)) > 1e-9)
                flow += d;
        }
        return flow;
    }

    // After max_flow, nodes reachable from s in residual graph form the S side of min-cut
    std::vector<int> min_cut_S(int s) {
        std::vector<int> S;
        for (int i = 0; i < n; ++i)
            if (level[i] >= 0) S.push_back(i);
        return S;
    }
};

std::vector<std::vector<int>> find_maxflow_cuts(const LPModel& model) {
    const int n = model.n;
    std::vector<std::vector<int>> cuts;

    // Build fractional capacity graph
    // For each visited node k, check if max-flow(0, k) < y_k
    for (int k = 1; k < n; ++k) {
        double yk = model.prim(model.y_col[k]);
        if (yk < 0.1) continue;  // skip unvisited nodes

        MaxFlow mf(n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (i == j || model.x_col[i][j] < 0) continue;
                double cap = model.prim(model.x_col[i][j]);
                if (cap > 1e-9)
                    mf.add_edge(i, j, cap);
            }

        double flow = mf.max_flow(0, k);
        if (flow < yk - 1e-6) {
            // Violated: min-cut separates depot from k
            auto S_side = mf.min_cut_S(0);
            // The cut set is V \ S_side (nodes not reachable from depot in residual)
            std::set<int> S_set(S_side.begin(), S_side.end());
            std::vector<int> cut_nodes;
            for (int i = 1; i < n; ++i)
                if (!S_set.count(i)) cut_nodes.push_back(i);
            if (!cut_nodes.empty())
                cuts.push_back(std::move(cut_nodes));
        }
    }
    return cuts;
}

// ── Subtour detection (Kosaraju SCC for asymmetric) ───────────────────────

// Returns sets of nodes that are either:
//   (a) a strongly-connected component not containing depot, or
//   (b) reachable from depot on the fractional graph but with no path back
// Both cases need a directed outgoing cut.
std::vector<std::vector<int>> find_subtours(const LPModel& model, double eps = 0.5) {
    const int n = model.n;

    // Build directed adjacency from fractional x values
    std::vector<std::vector<int>> adj(n), radj(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j && model.x_col[i][j] >= 0 && model.prim(model.x_col[i][j]) > eps) {
                adj[i].push_back(j);
                radj[j].push_back(i);
            }

    // Kosaraju pass 1: finish-order DFS on forward graph (iterative)
    std::vector<bool> visited(n, false);
    std::vector<int> order;
    order.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (visited[i]) continue;
        std::stack<std::pair<int,int>> stk;  // (node, adj_index)
        stk.push({i, 0});
        visited[i] = true;
        while (!stk.empty()) {
            auto& [u, idx] = stk.top();
            if (idx < static_cast<int>(adj[u].size())) {
                int v = adj[u][idx++];
                if (!visited[v]) { visited[v] = true; stk.push({v, 0}); }
            } else {
                order.push_back(u);
                stk.pop();
            }
        }
    }

    // Kosaraju pass 2: assign SCCs on reverse graph (iterative)
    std::vector<int> scc(n, -1);
    int ns = 0;
    for (int i = n - 1; i >= 0; --i) {
        int u = order[i];
        if (scc[u] >= 0) continue;
        std::stack<int> stk;
        stk.push(u);
        scc[u] = ns;
        while (!stk.empty()) {
            int cur = stk.top(); stk.pop();
            for (int v : radj[cur])
                if (scc[v] < 0) { scc[v] = ns; stk.push(v); }
        }
        ++ns;
    }

    // Collect SCCs not containing depot (node 0)
    int depot_scc = scc[0];
    std::vector<std::vector<int>> subtours;
    for (int c = 0; c < ns; ++c) {
        if (c == depot_scc) continue;
        std::vector<int> S;
        for (int i = 0; i < n; ++i)
            if (scc[i] == c) S.push_back(i);
        if (S.size() >= 2) subtours.push_back(std::move(S));
    }
    return subtours;
}

// ── Depot-unreachable detection (asymmetric only) ─────────────────────────
// Finds visited nodes from which depot (0) is not reachable following directed arcs.
// These are not SCCs but still need an outgoing cut toward the depot side.
std::vector<std::vector<int>> find_depot_unreachable(const LPModel& model, double eps = 0.5) {
    const int n = model.n;

    // Reverse graph: can we reach depot going backwards from depot?
    // i.e. which nodes can reach depot in the forward graph?
    std::vector<std::vector<int>> radj(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j && model.x_col[i][j] >= 0 && model.prim(model.x_col[i][j]) > eps)
                radj[j].push_back(i);

    // BFS backward from depot on reverse graph = all nodes that can reach depot
    std::vector<bool> can_reach_depot(n, false);
    std::queue<int> q;
    can_reach_depot[0] = true;
    q.push(0);
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : radj[u])
            if (!can_reach_depot[v]) { can_reach_depot[v] = true; q.push(v); }
    }

    // Any visited node that cannot reach depot is a violation
    std::vector<int> bad;
    for (int i = 1; i < n; ++i)
        if (!can_reach_depot[i] && model.prim(model.y_col[i]) > eps)
            bad.push_back(i);

    if (bad.empty()) return {};
    return {bad};  // Treat as one subset to cut
}

// ── Lifted cover cuts on the budget knapsack ──────────────────────────────

int find_and_add_cover_cuts(LPModel& lp, const Input& inp, int max_covers = 3) {
    const int n = lp.n;
    double B = inp.bud_raw;
    int cuts_added = 0;

    struct ArcInfo {
        int col;
        double weight;
        double lp_val;
    };
    std::vector<ArcInfo> arcs;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            if (lp.x_col[i][j] < 0) continue;
            double val = lp.prim(lp.x_col[i][j]);
            if (val < 1e-6) continue;
            double w = inp.cm[i][j];
            if (!std::isfinite(w) || w <= 0) continue;
            // Fatigue-aware weight: minimum possible fatigue cost of this arc
            // (arrival at i is at least C_0i on any feasible route)
            if (g_use_fatigue_covers && inp.fatigue_rate > 0 && i > 0
                && std::isfinite(inp.cm[0][i]) && inp.cm[0][i] > 0) {
                w *= (1.0 + inp.fatigue_rate * inp.cm[0][i] / inp.bud_raw);
            }
            arcs.push_back({lp.x_col[i][j], w, val});
        }

    if (arcs.empty()) return 0;

    std::sort(arcs.begin(), arcs.end(),
              [](const ArcInfo& a, const ArcInfo& b) { return a.lp_val > b.lp_val; });

    for (int attempt = 0; attempt < 5 && cuts_added < max_covers; ++attempt) {
        std::vector<int> cover_cols;
        double cover_weight = 0.0;
        double cover_lp_sum = 0.0;

        size_t start = attempt * (arcs.size() / 5);
        for (size_t idx = start; idx < arcs.size(); ++idx) {
            cover_cols.push_back(arcs[idx].col);
            cover_weight += arcs[idx].weight;
            cover_lp_sum += arcs[idx].lp_val;
            if (cover_weight > B) break;
        }

        if (cover_weight <= B) continue;

        int C_size = static_cast<int>(cover_cols.size());
        if (cover_lp_sum <= C_size - 1 + 1e-6) continue;

        std::vector<double> coeffs(C_size, 1.0);
        lp.add_row(-1e30, static_cast<double>(C_size - 1), cover_cols, coeffs);
        ++cuts_added;

        std::set<int> cover_set(cover_cols.begin(), cover_cols.end());
        double rhs = C_size - 1;

        for (const auto& arc : arcs) {
            if (cover_set.count(arc.col)) continue;
            if (arc.lp_val < 0.1) continue;
            double slack = cover_weight - B;
            int alpha = std::min(1, static_cast<int>(slack / std::max(arc.weight, 1e-9)));
            if (alpha <= 0) continue;
            auto lifted_cols = cover_cols;
            auto lifted_coeffs = std::vector<double>(C_size, 1.0);
            lifted_cols.push_back(arc.col);
            lifted_coeffs.push_back(static_cast<double>(alpha));
            lp.add_row(-1e30, rhs, lifted_cols, lifted_coeffs);
            ++cuts_added;
            if (cuts_added >= max_covers) break;
        }
    }

    if (cuts_added > 0)
        std::cerr << "  Added " << cuts_added << " cover cuts\n";
    return cuts_added;
}

// ── Hungarian algorithm (assignment relaxation) ────────────────────────────
// Returns the minimum cost assignment on an n×n cost matrix.
// Used to compute a lower bound on the TSP cost through a node subset.

double hungarian(const std::vector<std::vector<double>>& cost) {
    const int n = static_cast<int>(cost.size());
    if (n == 0) return 0.0;
    // Pad to 1-indexed
    const double INF = 1e18;
    std::vector<double> u(n + 1, 0), v(n + 1, 0);
    std::vector<int> p(n + 1, 0), way(n + 1, 0);

    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        std::vector<double> minv(n + 1, INF);
        std::vector<bool> used(n + 1, false);
        do {
            used[j0] = true;
            int i0 = p[j0], j1 = 0;
            double delta = INF;
            for (int j = 1; j <= n; ++j) {
                if (used[j]) continue;
                double cur = cost[i0 - 1][j - 1] - u[i0] - v[j];
                if (cur < minv[j]) { minv[j] = cur; way[j] = j0; }
                if (minv[j] < delta) { delta = minv[j]; j1 = j; }
            }
            for (int j = 0; j <= n; ++j) {
                if (used[j]) { u[p[j]] += delta; v[j] -= delta; }
                else { minv[j] -= delta; }
            }
            j0 = j1;
        } while (p[j0] != 0);
        do {
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0);
    }
    double total = 0.0;
    for (int j = 1; j <= n; ++j)
        total += cost[p[j] - 1][j - 1];
    return total;
}

// ── Routing Infeasibility Cuts ─────────────────────────────────────────────
// For a subset S of nodes, if the assignment relaxation lower bound on the
// cost of any cycle through S ∪ {depot} exceeds the budget B, then
// sum_{i in S} y_i <= |S| - 1 is a valid inequality.

int find_and_add_routing_cuts(LPModel& lp, const Input& inp, int max_cuts = 3) {
    const int n = lp.n;
    int cuts_added = 0;

    // Collect nodes the LP wants to visit, sorted by y* descending
    struct NodeInfo { int id; double y_val; };
    std::vector<NodeInfo> candidates;
    for (int i = 1; i < n; ++i) {
        double yv = lp.prim(lp.y_col[i]);
        if (yv > 0.5) candidates.push_back({i, yv});
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const NodeInfo& a, const NodeInfo& b) { return a.y_val > b.y_val; });

    if (candidates.size() < 3) return 0;

    // Greedily build infeasible sets
    std::vector<int> S_nodes;
    S_nodes.push_back(0); // depot always included in the assignment

    for (const auto& c : candidates) {
        S_nodes.push_back(c.id);

        if (S_nodes.size() < 4) continue; // need at least depot + 3 nodes

        // Build cost submatrix for S_nodes
        int sz = static_cast<int>(S_nodes.size());
        std::vector<std::vector<double>> sub_cost(sz, std::vector<double>(sz, 1e18));
        for (int i = 0; i < sz; ++i)
            for (int j = 0; j < sz; ++j) {
                if (i == j) { sub_cost[i][j] = 1e18; continue; }
                double cij = inp.cm[S_nodes[i]][S_nodes[j]];
                if (std::isfinite(cij)) sub_cost[i][j] = cij;
            }

        double ap_cost = hungarian(sub_cost);

        // Conservative fatigue adjustment: average fatigue multiplier is (1 + lambda/2)
        double fatigue_lb = ap_cost * (1.0 + inp.fatigue_rate / 2.0);

        if (fatigue_lb > inp.bud_raw) {
            // Check if the cut is violated: sum y_i > |S| - 1
            // S is S_nodes without depot
            std::vector<int> S_only; // nodes without depot
            double y_sum = 0.0;
            for (int k = 1; k < sz; ++k) { // skip depot at index 0
                S_only.push_back(S_nodes[k]);
                y_sum += lp.prim(lp.y_col[S_nodes[k]]);
            }
            int S_size = static_cast<int>(S_only.size());

            if (y_sum > S_size - 1 + 1e-6) {
                // Try to minimize S: remove nodes with smallest y* while still infeasible
                for (int attempt = 0; attempt < S_size && S_only.size() > 2; ++attempt) {
                    // Find node with smallest y* in S_only
                    int min_idx = 0;
                    double min_y = lp.prim(lp.y_col[S_only[0]]);
                    for (int k = 1; k < (int)S_only.size(); ++k) {
                        double yk = lp.prim(lp.y_col[S_only[k]]);
                        if (yk < min_y) { min_y = yk; min_idx = k; }
                    }
                    // Try removing it
                    std::vector<int> trial = {0}; // depot
                    for (int k = 0; k < (int)S_only.size(); ++k)
                        if (k != min_idx) trial.push_back(S_only[k]);

                    int tsz = static_cast<int>(trial.size());
                    std::vector<std::vector<double>> tc(tsz, std::vector<double>(tsz, 1e18));
                    for (int i = 0; i < tsz; ++i)
                        for (int j = 0; j < tsz; ++j) {
                            if (i == j) continue;
                            double cij = inp.cm[trial[i]][trial[j]];
                            if (std::isfinite(cij)) tc[i][j] = cij;
                        }
                    double trial_cost = hungarian(tc) * (1.0 + inp.fatigue_rate / 2.0);
                    if (trial_cost > inp.bud_raw) {
                        // Still infeasible, keep the removal
                        S_only.erase(S_only.begin() + min_idx);
                    } else {
                        break; // Can't shrink further
                    }
                }

                // Recheck violation after shrinking
                double y_sum2 = 0.0;
                for (int node : S_only) y_sum2 += lp.prim(lp.y_col[node]);
                int S2 = static_cast<int>(S_only.size());

                if (y_sum2 > S2 - 1 + 1e-6) {
                    // Add cut: sum y_i <= |S| - 1
                    std::vector<int> cols;
                    std::vector<double> coeffs;
                    for (int node : S_only) {
                        cols.push_back(lp.y_col[node]);
                        coeffs.push_back(1.0);
                    }
                    lp.add_row(-1e30, static_cast<double>(S2 - 1), cols, coeffs);
                    ++cuts_added;
                    std::cerr << "  Routing infeasibility cut: " << S2
                              << " nodes, violation=" << (y_sum2 - (S2 - 1)) << "\n";
                    if (cuts_added >= max_cuts) break;
                }
            }
            // Reset and try a different starting point
            S_nodes.clear();
            S_nodes.push_back(0);
        }
    }

    if (cuts_added > 0)
        std::cerr << "  Added " << cuts_added << " routing infeasibility cuts\n";
    return cuts_added;
}

// ── Cycle Cover Cuts (directed adaptation of Fischetti 1998) ───────────────
// If a set of directed arcs F forms a cycle with total cost > B, then:
//   sum_{(i,j) in F} x_ij <= sum_{v in V(F)} y_v - 1

int find_and_add_cycle_cover_cuts(LPModel& lp, const Input& inp, int max_cuts = 3) {
    const int n = lp.n;
    double B = inp.bud_raw;
    int cuts_added = 0;

    // Enumerate short directed cycles (length 3-5) in the support graph
    // For each cycle, check if cost > B and inequality is violated

    // Build adjacency from fractional solution
    std::vector<std::vector<int>> adj(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j && lp.x_col[i][j] >= 0 && lp.prim(lp.x_col[i][j]) > 0.1)
                adj[i].push_back(j);

    // Find triangles (3-cycles)
    for (int a = 0; a < n && cuts_added < max_cuts; ++a) {
        for (int b : adj[a]) {
            if (b <= a) continue; // avoid duplicates
            for (int c : adj[b]) {
                if (c <= a) continue;
                // Check if c -> a exists
                bool ca_exists = false;
                for (int nb : adj[c]) if (nb == a) { ca_exists = true; break; }
                if (!ca_exists) continue;

                // Cycle: a -> b -> c -> a
                double cost = inp.cm[a][b] + inp.cm[b][c] + inp.cm[c][a];
                // Fatigue-adjusted: conservative lower bound
                double fat_cost = cost * (1.0 + inp.fatigue_rate / 2.0);
                if (fat_cost <= B) continue;

                // Check violation: sum x >= sum y (would need sum x <= sum y - 1)
                double x_sum = lp.prim(lp.x_col[a][b]) + lp.prim(lp.x_col[b][c])
                             + lp.prim(lp.x_col[c][a]);
                double y_sum = lp.prim(lp.y_col[a]) + lp.prim(lp.y_col[b])
                             + lp.prim(lp.y_col[c]);
                if (x_sum <= y_sum - 1 + 1e-6) continue;

                // Violated — add cut
                std::vector<int> cols = {lp.x_col[a][b], lp.x_col[b][c], lp.x_col[c][a]};
                std::vector<double> coeffs = {1.0, 1.0, 1.0};
                // RHS: sum y_v - 1, so: sum x - sum y <= -1
                cols.push_back(lp.y_col[a]); coeffs.push_back(-1.0);
                cols.push_back(lp.y_col[b]); coeffs.push_back(-1.0);
                cols.push_back(lp.y_col[c]); coeffs.push_back(-1.0);
                lp.add_row(-1e30, -1.0, cols, coeffs);
                ++cuts_added;
                std::cerr << "  Cycle cover cut (3-cycle): nodes " << a << "," << b << "," << c
                          << " cost=" << fat_cost << "\n";
                if (cuts_added >= max_cuts) break;
            }
            if (cuts_added >= max_cuts) break;
        }
    }

    // Find 4-cycles
    for (int a = 0; a < n && cuts_added < max_cuts; ++a) {
        for (int b : adj[a]) {
            if (b == a) continue;
            for (int c : adj[b]) {
                if (c == a || c == b) continue;
                for (int d : adj[c]) {
                    if (d == a || d == b || d == c) continue;
                    // Check if d -> a exists
                    bool da_exists = false;
                    for (int nb : adj[d]) if (nb == a) { da_exists = true; break; }
                    if (!da_exists) continue;

                    double cost = inp.cm[a][b] + inp.cm[b][c] + inp.cm[c][d] + inp.cm[d][a];
                    double fat_cost = cost * (1.0 + inp.fatigue_rate / 2.0);
                    if (fat_cost <= B) continue;

                    double x_sum = lp.prim(lp.x_col[a][b]) + lp.prim(lp.x_col[b][c])
                                 + lp.prim(lp.x_col[c][d]) + lp.prim(lp.x_col[d][a]);
                    double y_sum = lp.prim(lp.y_col[a]) + lp.prim(lp.y_col[b])
                                 + lp.prim(lp.y_col[c]) + lp.prim(lp.y_col[d]);
                    if (x_sum <= y_sum - 1 + 1e-6) continue;

                    std::vector<int> cols = {lp.x_col[a][b], lp.x_col[b][c],
                                             lp.x_col[c][d], lp.x_col[d][a]};
                    std::vector<double> coeffs = {1.0, 1.0, 1.0, 1.0};
                    cols.push_back(lp.y_col[a]); coeffs.push_back(-1.0);
                    cols.push_back(lp.y_col[b]); coeffs.push_back(-1.0);
                    cols.push_back(lp.y_col[c]); coeffs.push_back(-1.0);
                    cols.push_back(lp.y_col[d]); coeffs.push_back(-1.0);
                    lp.add_row(-1e30, -1.0, cols, coeffs);
                    ++cuts_added;
                    std::cerr << "  Cycle cover cut (4-cycle): nodes "
                              << a << "," << b << "," << c << "," << d
                              << " cost=" << fat_cost << "\n";
                    if (cuts_added >= max_cuts) goto done_cycles;
                }
            }
        }
    }
    done_cycles:

    if (cuts_added > 0)
        std::cerr << "  Added " << cuts_added << " cycle cover cuts\n";
    return cuts_added;
}

// ── Path Inequality Cuts (directed adaptation of Fischetti 1998) ───────────
// For a directed path P = (i1 -> i2 -> ... -> ik) through non-depot nodes,
// define W(P) = {v : depot -> i1 -> ... -> ik -> v -> depot fits in budget}.
// Cut: sum_{arcs in P} x - sum_{v in V(P)} y + y_i1 + y_ik - sum_{v in W(P)} x_{ik,v} <= 0

int find_and_add_path_cuts(LPModel& lp, const Input& inp, int max_cuts = 3) {
    const int n = lp.n;
    int cuts_added = 0;

    // Build adjacency from fractional solution
    struct ArcInfo { int to; double x_val; };
    std::vector<std::vector<ArcInfo>> adj(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j && lp.x_col[i][j] >= 0) {
                double xv = lp.prim(lp.x_col[i][j]);
                if (xv > 0.1) adj[i].push_back({j, xv});
            }

    // Enumerate paths of length 2-4 starting from non-depot nodes
    for (int i1 = 1; i1 < n && cuts_added < max_cuts; ++i1) {
        if (lp.prim(lp.y_col[i1]) < 0.3) continue;

        for (const auto& [i2, x12] : adj[i1]) {
            if (i2 == 0) continue;

            // Path of length 2: i1 -> i2
            {
                // Cost from depot to i1 to i2
                double base_cost = inp.cm[0][i1] + inp.cm[i1][i2];
                // Fatigue-adjusted path cost
                double t0 = 0.0;
                double fat_path = inp.cm[0][i1] * (1.0 + inp.fatigue_rate * t0 / inp.bud_raw);
                double t1 = inp.cm[0][i1];
                fat_path += inp.cm[i1][i2] * (1.0 + inp.fatigue_rate * t1 / inp.bud_raw);
                double t2 = t1 + inp.cm[i1][i2];

                // Find W(P): nodes v such that i2 -> v -> 0 fits remaining budget
                double remaining = inp.bud_raw - fat_path;
                std::vector<int> W;
                for (int v = 1; v < n; ++v) {
                    if (v == i1 || v == i2) continue;
                    if (lp.x_col[i2][v] < 0) continue;
                    double leg_iv = inp.cm[i2][v] * (1.0 + inp.fatigue_rate * t2 / inp.bud_raw);
                    double t3 = t2 + inp.cm[i2][v];
                    double leg_v0 = inp.cm[v][0] * (1.0 + inp.fatigue_rate * t3 / inp.bud_raw);
                    if (leg_iv + leg_v0 <= remaining + 1e-9) W.push_back(v);
                }

                // Also check direct return: i2 -> 0
                double direct_return = inp.cm[i2][0] * (1.0 + inp.fatigue_rate * t2 / inp.bud_raw);
                bool can_return_direct = (direct_return <= remaining + 1e-9);

                // If W is empty and can't return directly, the path is infeasible
                // (handled by routing infeasibility cuts). Skip.
                // Path cut is useful when W is small but non-empty.
                if (W.empty() || can_return_direct) continue;

                // Check violation:
                // x_{i1,i2} - y_{i1} - y_{i2} + y_{i1} + y_{i2} - sum_{v in W} x_{i2,v} <= 0
                // Simplifies to: x_{i1,i2} - sum_{v in W} x_{i2,v} <= 0
                double lhs = lp.prim(lp.x_col[i1][i2]);
                for (int v : W) lhs -= lp.prim(lp.x_col[i2][v]);
                if (lhs <= 1e-6) continue;

                // Violated — add cut
                std::vector<int> cols = {lp.x_col[i1][i2]};
                std::vector<double> coeffs = {1.0};
                for (int v : W) {
                    cols.push_back(lp.x_col[i2][v]);
                    coeffs.push_back(-1.0);
                }
                lp.add_row(-1e30, 0.0, cols, coeffs);
                ++cuts_added;
                std::cerr << "  Path cut (len 2): " << i1 << "->" << i2
                          << " W=" << W.size() << " viol=" << lhs << "\n";
                if (cuts_added >= max_cuts) goto done_paths;
            }

            // Path of length 3: i1 -> i2 -> i3
            for (const auto& [i3, x23] : adj[i2]) {
                if (i3 == 0 || i3 == i1) continue;
                if (cuts_added >= max_cuts) goto done_paths;

                double t0 = 0.0;
                double fat_path = inp.cm[0][i1] * (1.0 + inp.fatigue_rate * t0 / inp.bud_raw);
                double t1 = inp.cm[0][i1];
                fat_path += inp.cm[i1][i2] * (1.0 + inp.fatigue_rate * t1 / inp.bud_raw);
                double t2 = t1 + inp.cm[i1][i2];
                fat_path += inp.cm[i2][i3] * (1.0 + inp.fatigue_rate * t2 / inp.bud_raw);
                double t3 = t2 + inp.cm[i2][i3];

                double remaining = inp.bud_raw - fat_path;
                if (remaining < 0) continue; // path itself exceeds budget

                std::vector<int> W;
                for (int v = 1; v < n; ++v) {
                    if (v == i1 || v == i2 || v == i3) continue;
                    if (lp.x_col[i3][v] < 0) continue;
                    double leg_iv = inp.cm[i3][v] * (1.0 + inp.fatigue_rate * t3 / inp.bud_raw);
                    double t4 = t3 + inp.cm[i3][v];
                    double leg_v0 = inp.cm[v][0] * (1.0 + inp.fatigue_rate * t4 / inp.bud_raw);
                    if (leg_iv + leg_v0 <= remaining + 1e-9) W.push_back(v);
                }

                double direct_return = inp.cm[i3][0] * (1.0 + inp.fatigue_rate * t3 / inp.bud_raw);
                bool can_return_direct = (direct_return <= remaining + 1e-9);
                if (W.empty() || can_return_direct) continue;

                // Path inequality for P = (i1, i2, i3):
                // x_{i1,i2} + x_{i2,i3} - y_{i1} - y_{i2} - y_{i3} + y_{i1} + y_{i3}
                //   - sum_{v in W} x_{i3,v} <= 0
                // Simplifies to: x_{i1,i2} + x_{i2,i3} - y_{i2} - sum_{v in W} x_{i3,v} <= 0
                double lhs = lp.prim(lp.x_col[i1][i2]) + lp.prim(lp.x_col[i2][i3])
                           - lp.prim(lp.y_col[i2]);
                for (int v : W) lhs -= lp.prim(lp.x_col[i3][v]);
                if (lhs <= 1e-6) continue;

                std::vector<int> cols = {lp.x_col[i1][i2], lp.x_col[i2][i3]};
                std::vector<double> coeffs = {1.0, 1.0};
                cols.push_back(lp.y_col[i2]); coeffs.push_back(-1.0);
                for (int v : W) {
                    cols.push_back(lp.x_col[i3][v]);
                    coeffs.push_back(-1.0);
                }
                lp.add_row(-1e30, 0.0, cols, coeffs);
                ++cuts_added;
                std::cerr << "  Path cut (len 3): " << i1 << "->" << i2 << "->" << i3
                          << " W=" << W.size() << " viol=" << lhs << "\n";
            }
        }
    }
    done_paths:

    if (cuts_added > 0)
        std::cerr << "  Added " << cuts_added << " path cuts\n";
    return cuts_added;
}

// ── Route extraction & validation ──────────────────────────────────────────

std::vector<int> extract_route(const LPModel& model, double eps = 0.5) {
    const int n = model.n;
    std::vector<int> succ(n, -1);
    for (int i = 0; i < n; ++i) {
        int found = 0;
        for (int j = 0; j < n; ++j) {
            if (i != j && model.x_col[i][j] >= 0 && model.prim(model.x_col[i][j]) > eps) {
                if (succ[i] == -1) succ[i] = j;
                ++found;
            }
        }
        if (found > 1)
            std::cerr << "Warning: node " << i << " has " << found << " outgoing arcs > eps\n";
    }
    std::vector<int> route;
    int cur = succ[0];
    while (cur > 0 && cur != -1 && static_cast<int>(route.size()) < n) {
        route.push_back(cur);
        cur = succ[cur];
    }
    return route;
}

bool is_feasible_route(const Input& inp, const std::vector<int>& route) {
    return rcost_fatigue(inp.cm, route, inp.bud_raw, inp.fatigue_rate) <= inp.bud_raw;
}

// ── Greedy heuristic ───────────────────────────────────────────────────────

std::vector<int> greedy_route(const Input& inp) {
    int n = static_cast<int>(inp.pts.size());
    std::vector<bool> visited(n, false);
    visited[0] = true;
    std::vector<int> route;
    double cost = 0.0, elapsed = 0.0;
    int cur = 0;
    // Use worst-case fatigue budget for consistency with LP
    double bud_lp = inp.bud_raw / (1.0 + inp.fatigue_rate);

    while (true) {
        int best_j = -1;
        double best_ratio = -1.0;
        for (int j = 1; j < n; ++j) {
            if (visited[j]) continue;
            double go = inp.cm[cur][j], back = inp.cm[j][0];
            if (!std::isfinite(go) || !std::isfinite(back)) continue;
            if (cost + go + back > bud_lp) continue;

            double fm_go = 1.0 + inp.fatigue_rate * (elapsed / std::max(inp.bud_raw, 1.0));
            double fat_go = go * fm_go;
            double fm_back = 1.0 + inp.fatigue_rate * ((elapsed + go) / std::max(inp.bud_raw, 1.0));
            double fat_back = fat_go + back * fm_back;
            if (fat_back > inp.bud_raw) continue;

            double ratio = inp.pts[j] / std::max(go, 1e-9);
            if (ratio > best_ratio) {
                best_ratio = ratio;
                best_j = j;
            }
        }
        if (best_j < 0) break;

        double go = inp.cm[cur][best_j];
        elapsed += go;
        cost += go;
        visited[best_j] = true;
        route.push_back(best_j);
        cur = best_j;
    }
    return route;
}

// ── Simulated Annealing ────────────────────────────────────────────────────

static double rpts(const std::vector<double>& pts, const std::vector<int>& route) {
    double s = 0.0;
    for (int v : route) s += pts[v];
    return s;
}

std::vector<int> solve_sa(const Input& inp, int n_iterations = 80000,
                          double temp_start = 100.0, double temp_end = 0.1,
                          unsigned seed = 42) {
    const int n = static_cast<int>(inp.pts.size());
    std::mt19937 rng(seed);
    auto randu = [&]() { return std::uniform_real_distribution<double>(0,1)(rng); };
    auto randi = [&](int lo, int hi) { return std::uniform_int_distribution<int>(lo,hi)(rng); };

    std::vector<int> all_ctrls;
    for (int i = 1; i < n; ++i) all_ctrls.push_back(i);

    auto route = greedy_route(inp);
    std::vector<bool> visited(n, false);
    visited[0] = true;
    for (int v : route) visited[v] = true;

    std::vector<int> best_route = route;
    double best_score = rpts(inp.pts, route);
    double cur_score = best_score;
    double cur_cost  = rcost(inp.cm, route);
    double temp = temp_start;
    double decay = std::pow(temp_end / temp_start, 1.0 / std::max(n_iterations, 1));

    // move probabilities matching Python tuned: (0.30, 0.30, 0.20, 0.20)
    const double t1 = 0.30, t2 = 0.60, t3 = 0.80;

    for (int it = 0; it < n_iterations; ++it) {
        temp *= decay;
        auto new_route   = route;
        auto new_visited = visited;
        double mv = randu();

        if (mv < t1 && !new_route.empty()) {
            // remove random node
            int idx = randi(0, static_cast<int>(new_route.size()) - 1);
            new_visited[new_route[idx]] = false;
            new_route.erase(new_route.begin() + idx);

        } else if (mv < t2) {
            // insert unvisited node at best position using O(n) delta computation
            std::vector<int> unv;
            for (int j : all_ctrls) if (!new_visited[j]) unv.push_back(j);
            if (unv.empty()) continue;
            int j = unv[randi(0, static_cast<int>(unv.size()) - 1)];
            int bp = 0; double bc = std::numeric_limits<double>::infinity();
            for (int pos = 0; pos <= static_cast<int>(new_route.size()); ++pos) {
                int prev = (pos == 0) ? 0 : new_route[pos - 1];
                int next = (pos == static_cast<int>(new_route.size())) ? 0 : new_route[pos];
                double delta = inp.cm[prev][j] + inp.cm[j][next] - inp.cm[prev][next];
                if (delta < bc) { bc = delta; bp = pos; }
            }
            new_route.insert(new_route.begin() + bp, j);
            new_visited[j] = true;

        } else if (mv < t3 && new_route.size() >= 2) {
            // Or-opt: remove a segment of length 1, 2, or 3 and reinsert elsewhere
            int seg_len = std::min(randi(1, 3), static_cast<int>(new_route.size()));
            if (static_cast<int>(new_route.size()) > seg_len) {
                int seg_start = randi(0, static_cast<int>(new_route.size()) - seg_len);
                // extract segment
                std::vector<int> seg(new_route.begin() + seg_start,
                                     new_route.begin() + seg_start + seg_len);
                new_route.erase(new_route.begin() + seg_start,
                                new_route.begin() + seg_start + seg_len);
                // reinsert at a different position
                int insert_pos = randi(0, static_cast<int>(new_route.size()));
                new_route.insert(new_route.begin() + insert_pos, seg.begin(), seg.end());
            }

        } else {
            // swap visited for unvisited
            std::vector<int> unv;
            for (int j : all_ctrls) if (!new_visited[j]) unv.push_back(j);
            if (unv.empty() || new_route.empty()) continue;
            int oi  = randi(0, static_cast<int>(new_route.size()) - 1);
            int old = new_route[oi];
            int nw  = unv[randi(0, static_cast<int>(unv.size()) - 1)];
            new_route[oi] = nw;
            new_visited[old] = false;
            new_visited[nw]  = true;
        }

        double nc = rcost(inp.cm, new_route);
        if (nc > inp.bud_eff) continue;
        if (rcost_fatigue(inp.cm, new_route, inp.bud_raw, inp.fatigue_rate) > inp.bud_raw) continue;

        double ns    = rpts(inp.pts, new_route);
        double delta = ns - cur_score;
        if (delta > 0 || (delta == 0 && nc < cur_cost) ||
            (delta < 0 && temp > 1e-6 && randu() < std::exp(delta / temp))) {
            route = std::move(new_route);
            visited = std::move(new_visited);
            cur_score = ns; cur_cost = nc;
            if (ns > best_score) { best_score = ns; best_route = route; }
        }
    }

    // repair: drop nodes that violate fatigue budget
    while (!best_route.empty() &&
           rcost_fatigue(inp.cm, best_route, inp.bud_raw, inp.fatigue_rate) > inp.bud_raw) {
        int worst = static_cast<int>(std::min_element(best_route.begin(), best_route.end(),
            [&](int a, int b){ return inp.pts[a] < inp.pts[b]; }) - best_route.begin());
        best_route.erase(best_route.begin() + worst);
    }
    return best_route;
}

std::vector<int> solve_sa_iterated(const Input& inp, int n_restarts = -1,
                                    int n_iterations = -1) {
    const int n = static_cast<int>(inp.pts.size());
    if (n_iterations < 0) n_iterations = std::max(10000, std::min(120000, n * 2500));
    if (n_restarts   < 0) n_restarts   = std::max(40,    std::min(200,    3000 / n));
    std::cerr << "SA: n=" << n << " iters=" << n_iterations
              << " restarts=" << n_restarts
              << " total=" << (long long)n_iterations * n_restarts << "\n";
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

// ── Branch-and-Cut Node ────────────────────────────────────────────────────

struct BNCNode {
    std::vector<std::pair<int, double>> fixings;  // col -> value (0 or 1)
    double ub = 0.0;
};

struct Solver {
    const Input& inp;
    LPModel root;
    double best_pts = 0.0;
    std::vector<int> best_route;
    int max_cuts = 20;
    int max_depth = 15;
    double time_limit_s = 900.0;

    bool proved_optimal = false;
    double best_ub = std::numeric_limits<double>::infinity();

    // Ablation flags — toggle cut families on/off
    bool use_tightened_coupling = true;
    bool use_fatigue_covers = true;
    bool use_routing_infeasibility = true;
    bool use_cycle_covers = true;
    bool use_path_ineq = true;

    explicit Solver(const Input& i) : inp(i) {
        root.build(inp);
    }

    void solve(double warm_start_pts = 0.0, std::vector<int> warm_start_route = {}) {
        // Warm start from provided solution (e.g. SA result)
        if (warm_start_pts > best_pts) {
            best_pts = warm_start_pts;
            best_route = std::move(warm_start_route);
        }
        // Also try greedy
        auto gr = greedy_route(inp);
        double gr_pts = std::accumulate(gr.begin(), gr.end(), 0.0,
            [&](double s, int v) { return s + inp.pts[v]; });
        if (gr_pts > best_pts) { best_pts = gr_pts; best_route = std::move(gr); }
        std::cerr << "B&C warm start: " << best_pts << " pts\n";

        // B&C search
        std::stack<BNCNode> node_stack;
        BNCNode root_node;
        std::cerr << "Solving root LP...\n";
        bool root_ok = root.solve();
        std::cerr << "Root LP solve returned: " << root_ok
                  << "  model_status=" << Highs_getModelStatus(root.highs)
                  << "  obj=" << (root_ok ? root.obj() : -1.0) << "\n";
        root_node.ub = root_ok ? root.obj() : -std::numeric_limits<double>::infinity();
        std::cerr << "Root UB=" << root_node.ub << "  best_pts=" << best_pts << "\n";
        if (root_node.ub > best_pts) node_stack.push(std::move(root_node));

        int nodes = 0;
        auto t_start = std::chrono::steady_clock::now();
        while (!node_stack.empty() && nodes++ < 10000) {
            double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
            if (elapsed > time_limit_s) { std::cerr << "Time limit reached\n"; break; }
            BNCNode node = std::move(node_stack.top());
            node_stack.pop();
            process_node(std::move(node), node_stack);
        }
        proved_optimal = node_stack.empty() && nodes < 10000;
        if (proved_optimal) {
            best_ub = best_pts;
        } else {
            best_ub = best_pts;
            std::stack<BNCNode> tmp = node_stack;
            while (!tmp.empty()) {
                if (tmp.top().ub > best_ub) best_ub = tmp.top().ub;
                tmp.pop();
            }
        }
        double gap_pct = best_pts > 0 ? 100.0 * (best_ub - best_pts) / best_pts : 0.0;
        std::cerr << "Processed " << nodes << " nodes, best: " << best_pts
                  << " pts, UB: " << best_ub << ", gap: " << gap_pct << "%\n";
    }

    void process_node(BNCNode node, std::stack<BNCNode>& node_stack) {
        if (node.ub <= best_pts + 1e-6 || node.fixings.size() > max_depth) return;

        // Clone LP and apply fixings
        LPModel lp;
        lp.clone_from(root, node.fixings);

        // Cut loop — accumulate cuts, do NOT delete them between iterations
        for (int cut_iter = 0; cut_iter < max_cuts; ++cut_iter) {
            if (!lp.solve()) return;
            double lp_ub = lp.obj();
            std::cerr << "LP UB: " << lp_ub << " (gap: " << (lp_ub-best_pts) << ")\n";  
            if (lp_ub <= best_pts + 1e-6) return;

            auto subtours = find_subtours(lp);
            auto unreachable = find_depot_unreachable(lp);
            subtours.insert(subtours.end(), unreachable.begin(), unreachable.end());

            // Max-flow separation: find violated connectivity cuts missed by rounding
            auto mf_cuts = find_maxflow_cuts(lp);
            subtours.insert(subtours.end(), mf_cuts.begin(), mf_cuts.end());

            int cover_cuts = find_and_add_cover_cuts(lp, inp);

            // Routing infeasibility cuts (on y variables)
            int routing_cuts = use_routing_infeasibility
                ? find_and_add_routing_cuts(lp, inp) : 0;

            // Cycle cover cuts (directed Fischetti 1998)
            int cycle_cuts = use_cycle_covers
                ? find_and_add_cycle_cover_cuts(lp, inp) : 0;

            // Path inequality cuts (directed Fischetti 1998)
            int path_cuts = use_path_ineq
                ? find_and_add_path_cuts(lp, inp) : 0;

            if (subtours.empty() && cover_cuts == 0 && routing_cuts == 0
                && cycle_cuts == 0 && path_cuts == 0) break;

            for (const auto& S : subtours)
                lp.add_sec(S);
        }

        if (!lp.solve()) return;
        double lp_ub = lp.obj();
        if (lp_ub <= best_pts + 1e-6) return;

        // Check integrality and select branching variable
        bool integer_sol = true;
        std::vector<std::pair<int,int>> candidates;  // (node_index, col)
        for (int i = 1; i < root.n; ++i) {
            double v = lp.prim(root.y_col[i]);
            double frac = std::min(v, 1.0 - v);
            if (frac > 1e-5) {
                integer_sol = false;
                candidates.push_back({i, root.y_col[i]});
            }
        }

        if (integer_sol) {
            // Extract and validate route
            auto route = extract_route(lp);
            double pts = std::accumulate(route.begin(), route.end(), 0.0,
                [&](double sum, int v) { return sum + inp.pts[v]; });
            if (is_feasible_route(inp, route) && pts > best_pts + 1e-6) {
                best_pts = pts;
                best_route = std::move(route);
                std::cerr << "New best: " << best_pts << " pts (" << best_route.size() << " nodes)\n";
            }
        } else if (!candidates.empty()) {
            int branch_col = -1;

            // Reliability branching: use strong branching for first few decisions,
            // then fall back to pseudocost branching
            if (node.fixings.size() < 8 && candidates.size() <= 20) {
                // Strong branching: tentatively branch on each candidate, pick best
                double best_score = -1.0;
                for (const auto& [ni, col] : candidates) {
                    // Try fixing to 0
                    auto fixings0 = node.fixings;
                    fixings0.emplace_back(col, 0.0);
                    LPModel lp0;
                    lp0.clone_from(root, fixings0);
                    double ub0 = lp0.solve() ? lp0.obj() : -1e30;

                    // Try fixing to 1
                    auto fixings1 = node.fixings;
                    fixings1.emplace_back(col, 1.0);
                    LPModel lp1;
                    lp1.clone_from(root, fixings1);
                    double ub1 = lp1.solve() ? lp1.obj() : -1e30;

                    // Score: product of bound improvements (standard strong branching score)
                    double down = std::max(lp_ub - ub0, 1e-6);
                    double up   = std::max(lp_ub - ub1, 1e-6);
                    double score = (1e-6 + down) * (1e-6 + up);
                    if (score > best_score) {
                        best_score = score;
                        branch_col = col;
                    }
                }
            } else {
                // Most-fractional fallback
                double max_frac = 0.0;
                for (const auto& [ni, col] : candidates) {
                    double v = lp.prim(col);
                    double frac = std::min(v, 1.0 - v);
                    if (frac > max_frac) { max_frac = frac; branch_col = col; }
                }
            }

            if (branch_col > 0) {
                BNCNode node0 = node, node1 = node;
                node0.fixings.emplace_back(branch_col, 0.0);
                node0.ub = lp_ub;
                node1.fixings.emplace_back(branch_col, 1.0);
                node1.ub = lp_ub;
                node_stack.push(std::move(node0));
                node_stack.push(std::move(node1));
            }
        }
    }
};

// ── Main ───────────────────────────────────────────────────────────────────

static void run_map(const std::string& in_path, const std::string& out_path) {
    std::cerr << "\n=== " << in_path << " ===\n";
    Input inp = parse_input(read_file(in_path));

    auto t_sa = std::chrono::steady_clock::now();
    auto sa_route = solve_sa_iterated(inp);
    double sa_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_sa).count();
    double sa_pts     = rpts(inp.pts, sa_route);
    double sa_base    = rcost(inp.cm, sa_route);
    double sa_fatigue = rcost_fatigue(inp.cm, sa_route, inp.bud_raw, inp.fatigue_rate);
    std::cerr << "SA: " << sa_pts << " pts (" << sa_route.size() << " nodes) in " << sa_elapsed << "s\n";

    auto t_bnc = std::chrono::steady_clock::now();
    Solver solver(inp);
    solver.solve(sa_pts, sa_route);
    double bnc_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_bnc).count();
    double bnc_base    = rcost(inp.cm, solver.best_route);
    double bnc_fatigue = rcost_fatigue(inp.cm, solver.best_route, inp.bud_raw, inp.fatigue_rate);
    std::cerr << "B&C(SA): " << solver.best_pts << " pts (" << solver.best_route.size() << " nodes) in " << bnc_elapsed << "s\n";

    std::ofstream out(out_path);
    out << "{\n";
    out << "  \"sa\": {\"pts\": " << sa_pts << ", \"nodes\": " << sa_route.size()
        << ", \"elapsed_s\": " << sa_elapsed << ", \"base_cost\": " << sa_base
        << ", \"fatigue_cost\": " << sa_fatigue << ", \"route\": [";
    for (size_t i = 0; i < sa_route.size(); ++i) { if (i) out << ", "; out << sa_route[i]; }
    out << "]},\n";
    out << "  \"bnc_sa\": {\"pts\": " << solver.best_pts << ", \"nodes\": " << solver.best_route.size()
        << ", \"elapsed_s\": " << bnc_elapsed
        << ", \"proved_optimal\": " << (solver.proved_optimal ? "true" : "false")
        << ", \"best_ub\": " << solver.best_ub
        << ", \"gap_pct\": " << (solver.best_pts > 0 ? 100.0 * (solver.best_ub - solver.best_pts) / solver.best_pts : 0.0)
        << ", \"base_cost\": " << bnc_base
        << ", \"fatigue_cost\": " << bnc_fatigue << ", \"route\": [";
    for (size_t i = 0; i < solver.best_route.size(); ++i) { if (i) out << ", "; out << solver.best_route[i]; }
    out << "]}\n";
    out << "}\n";
}

struct MapResult {
    std::string name;
    double sa_pts, bnc_pts;
    int    sa_nodes, bnc_nodes;
    double sa_s, bnc_s;
    bool   bnc_optimal;
    double gap_pct;
};

int main(int argc, char* argv[]) {
    std::string input_dir = "instances";
    std::string config_name = "all_cuts";
    std::string csv_file;

    // Parse command-line flags
    bool flag_routing = true, flag_cycle = true, flag_path = true;
    for (int a = 1; a < argc; ++a) {
        std::string arg = argv[a];
        if (arg == "--coupling=off")        g_use_tightened_coupling = false;
        else if (arg == "--coupling=on")    g_use_tightened_coupling = true;
        else if (arg == "--fatigue-covers=off") g_use_fatigue_covers = false;
        else if (arg == "--fatigue-covers=on")  g_use_fatigue_covers = true;
        else if (arg == "--routing=off")    flag_routing = false;
        else if (arg == "--routing=on")     flag_routing = true;
        else if (arg == "--cycle=off")      flag_cycle = false;
        else if (arg == "--cycle=on")       flag_cycle = true;
        else if (arg == "--path=off")       flag_path = false;
        else if (arg == "--path=on")        flag_path = true;
        else if (arg.rfind("--config=", 0) == 0) config_name = arg.substr(9);
        else if (arg.rfind("--csv=", 0) == 0)    csv_file = arg.substr(6);
        else input_dir = arg;
    }

    std::cerr << "Config: " << config_name
              << " coupling=" << (g_use_tightened_coupling ? "on" : "off")
              << " fatigue-covers=" << (g_use_fatigue_covers ? "on" : "off")
              << " routing=" << (flag_routing ? "on" : "off")
              << " cycle=" << (flag_cycle ? "on" : "off")
              << " path=" << (flag_path ? "on" : "off") << "\n";

    // Scan input directory for op_input_*.json files
    std::vector<std::pair<std::string,std::string>> maps;
    for (const auto& entry : std::filesystem::directory_iterator(input_dir)) {
        std::string fname = entry.path().filename().string();
        if (fname.substr(0, 9) == "op_input_" && fname.size() > 14 &&
            fname.substr(fname.size() - 5) == ".json") {
            std::string base = fname.substr(9, fname.size() - 14);
            std::string in_path = entry.path().string();
            std::string out_path = (entry.path().parent_path() /
                                    ("op_output_" + base + ".json")).string();
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
        try {
            std::cerr << "\n=== " << in << " ===\n";
            Input inp = parse_input(read_file(in));

            // SA
            auto t0 = std::chrono::steady_clock::now();
            auto sa_route = solve_sa_iterated(inp);
            double sa_s = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t0).count();
            double sa_pts = 0;
            for (int v : sa_route) sa_pts += inp.pts[v];

            // B&C with ablation flags
            auto t1 = std::chrono::steady_clock::now();
            Solver solver(inp);
            solver.use_routing_infeasibility = flag_routing;
            solver.use_cycle_covers = flag_cycle;
            solver.use_path_ineq = flag_path;
            solver.solve(sa_pts, sa_route);
            double bnc_s = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t1).count();

            MapResult r;
            std::string fname = std::filesystem::path(in).filename().string();
            r.name = fname.substr(9, fname.size() - 14);
            r.sa_pts = sa_pts;
            r.sa_nodes = static_cast<int>(sa_route.size());
            r.sa_s = sa_s;
            r.bnc_pts = solver.best_pts;
            r.bnc_nodes = static_cast<int>(solver.best_route.size());
            r.bnc_s = bnc_s;
            r.bnc_optimal = solver.proved_optimal;
            r.gap_pct = solver.best_pts > 0
                ? 100.0 * (solver.best_ub - solver.best_pts) / solver.best_pts : 0.0;
            results.push_back(r);
        }
        catch (const std::exception& e) {
            std::cerr << "Error on " << in << ": " << e.what() << '\n';
        }
    }

    // ── Summary table ──
    std::cout << "\n";
    std::cout << "+----------------------------------------+---------------------------+----------------------------------------------+\n";
    std::cout << "| Instance                               | SA                        | B&C(SA)                                      |\n";
    std::cout << "|                                        |  pts  nodes     time      |  pts  nodes     time    optimal?    gap       |\n";
    std::cout << "+----------------------------------------+---------------------------+----------------------------------------------+\n";
    for (const auto& r : results) {
        std::cout << "| " << std::left  << std::setw(38) << r.name
                  << " | " << std::right << std::setw(5)  << static_cast<int>(r.sa_pts)
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

    // Summary stats
    int proven = 0, total = static_cast<int>(results.size());
    double sum_gap = 0, max_gap = 0;
    for (const auto& r : results) {
        if (r.bnc_optimal) ++proven;
        sum_gap += r.gap_pct;
        if (r.gap_pct > max_gap) max_gap = r.gap_pct;
    }
    double mean_gap = total > 0 ? sum_gap / total : 0;
    std::cout << "\nProven optimal: " << proven << "/" << total
              << "  Mean gap: " << std::setprecision(2) << mean_gap << "%"
              << "  Max gap: " << max_gap << "%\n";

    // CSV output for batch collection
    if (!csv_file.empty()) {
        bool file_exists = std::filesystem::exists(csv_file);
        std::ofstream csv(csv_file, std::ios::app);
        if (!file_exists) {
            csv << "config,instance,sa_pts,bnc_pts,gap_pct,optimal,time_s\n";
        }
        for (const auto& r : results) {
            csv << config_name << "," << r.name << ","
                << static_cast<int>(r.sa_pts) << ","
                << static_cast<int>(r.bnc_pts) << ","
                << std::fixed << std::setprecision(2) << r.gap_pct << ","
                << (r.bnc_optimal ? "YES" : "no") << ","
                << std::setprecision(1) << r.bnc_s << "\n";
        }
    }

    return 0;
}