
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
#include <array>
#include "interfaces/highs_c_api.h"
#include "lp_data/HConst.h"
//  #include "glpk.h"

// ΓöÇΓöÇ JSON parser (improved from original) ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ

static std::string read_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot open " + path);
    return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

struct Input {
    std::vector<std::vector<double>> cm;
    std::vector<std::vector<double>> gain;   // phi_plus: cumulative uphill metres per arc (i->j)
    std::vector<std::vector<double>> loss;   // phi_minus: cumulative downhill metres per arc (i->j)
    std::vector<std::vector<double>> dist;   // cumulative horizontal path length (metres) per arc (i->j)
    std::vector<double> pts;
    double bud_eff = 0.0, bud_raw = 0.0, fatigue_rate = 0.0;
    double rho = 0.5;  // recovery fraction on downhill; sensitivity parameter, see main()
    double mu = 0.0;   // distance-to-elevation fatigue weight, derived not swept (see solver_flow_torremocha.cpp)
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
    find_key("bud_eff"); i = skip_ws(i); inp.bud_eff = parse_number(i);
    find_key("bud_raw"); i = skip_ws(i); inp.bud_raw = parse_number(i);
    find_key("fatigue_rate"); i = skip_ws(i); inp.fatigue_rate = parse_number(i);

    const int n = static_cast<int>(inp.pts.size());
    if (find_key_opt("gain")) {
        inp.gain = parse_array2d(i);
    } else {
        std::cerr << "  [warn] input has no 'gain' matrix (old format) -- "
                     "asymmetric fatigue model disabled, phi_plus=0\n";
        inp.gain.assign(n, std::vector<double>(n, 0.0));
    }
    if (find_key_opt("loss")) {
        inp.loss = parse_array2d(i);
    } else {
        std::cerr << "  [warn] input has no 'loss' matrix (old format) -- phi_minus=0\n";
        inp.loss.assign(n, std::vector<double>(n, 0.0));
    }
    if (find_key_opt("dist")) {
        inp.dist = parse_array2d(i);
    } else {
        std::cerr << "  [warn] input has no 'dist' matrix (old format) -- distance fatigue term=0\n";
        inp.dist.assign(n, std::vector<double>(n, 0.0));
    }
    if (find_key_opt("rho_default")) {
        i = skip_ws(i); inp.rho = parse_number(i);
    }
    if (find_key_opt("mu_default")) {
        i = skip_ws(i); inp.mu = parse_number(i);
    }

    return inp;
}


// ΓöÇΓöÇ Cost helpers ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ

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

// ── Non-linear asymmetric fatigue model (Sec 3.6 rework) ──────────────
// psi_ij = phi_plus_ij - rho * phi_minus_ij: net fatigue contribution of
// arc (i,j), where phi_plus/phi_minus are cumulative uphill/downhill
// metres along the cost-optimal path (see python/core/pathfinding.py).
// rho in [0,1] is the downhill-recovery fraction, treated as a sensitivity
// parameter (see rho sweep in main()), not calibrated from a single value.
static inline double psi_arc(const Input& inp, int i, int j, double rho) {
    return inp.gain[i][j] - rho * inp.loss[i][j] + inp.mu * inp.dist[i][j];
}

// Which arcs are structurally admissible for the fatigue-bound graph.
// Mirrors the base-cost pre-filter in LPModel::build() but does NOT use
// the fatigue-aware filter, since that filter itself depends on Ghat --
// this would be circular. Using only the cost-only filter here is
// conservative (a superset of true survivors), which keeps Ghat valid.
static inline bool arc_survives_base(const Input& inp, int i, int j) {
    if (i == j) return false;
    if (!std::isfinite(inp.cm[i][j])) return false;
    if (!std::isfinite(inp.cm[j][0])) return false;
    if (inp.cm[i][j] + inp.cm[j][0] > inp.bud_raw) return false;
    if (inp.cm[0][i] + inp.cm[i][j] + inp.cm[j][0] > inp.bud_raw) return false;
    return true;
}

// Ghat_i: valid upper bound on the (unclipped) cumulative fatigue state
// G_i reachable at node i by any budget-admissible route from the depot.
// Bellman-Ford relaxation in the (max,+) semiring, capped at n-1 rounds so
// it stays well-defined even with positive-weight cycles in the psi graph
// (see solver_flow_torremocha.cpp for the full derivation/comment).
static std::vector<double> compute_fatigue_bounds(const Input& inp, double rho) {
    const int n = static_cast<int>(inp.pts.size());
    const double NEG_INF = -1e30;
    std::vector<double> g(n, NEG_INF);
    g[0] = 0.0;

    std::vector<std::array<double,3>> arcs;
    arcs.reserve(n * n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (arc_survives_base(inp, i, j))
                arcs.push_back({double(i), double(j), psi_arc(inp, i, j, rho)});

    for (int round = 0; round < n - 1; ++round) {
        bool changed = false;
        for (const auto& a : arcs) {
            int i = int(a[0]), j = int(a[1]);
            double psi = a[2];
            if (g[i] <= NEG_INF / 2) continue;
            double cand = g[i] + psi;
            if (cand > g[j] + 1e-12) { g[j] = cand; changed = true; }
        }
        if (!changed) break;
    }
    for (int i = 0; i < n; ++i) g[i] = std::max(g[i], 0.0);
    return g;
}

// Ghat_low_i: valid LOWER bound on the (unclipped) cumulative fatigue
// state G_i. Needed because psi_ij can be negative (net downhill
// recovery) -- see add_fatigue_flow_coupling for why g_col must be
// allowed below 0 for this to be a sound relaxation (see
// solver_flow_torremocha.cpp for the full derivation/comment). Same
// construction as compute_fatigue_bounds but minimising, capped at n-1
// rounds; no clipping at the end (only unreached nodes fall back to 0).
static std::vector<double> compute_fatigue_lower_bounds(const Input& inp, double rho) {
    const int n = static_cast<int>(inp.pts.size());
    const double POS_INF = 1e30;
    std::vector<double> g(n, POS_INF);
    g[0] = 0.0;

    std::vector<std::array<double,3>> arcs;
    arcs.reserve(n * n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (arc_survives_base(inp, i, j))
                arcs.push_back({double(i), double(j), psi_arc(inp, i, j, rho)});

    for (int round = 0; round < n - 1; ++round) {
        bool changed = false;
        for (const auto& a : arcs) {
            int i = int(a[0]), j = int(a[1]);
            double psi = a[2];
            if (g[i] >= POS_INF / 2) continue;
            double cand = g[i] + psi;
            if (cand < g[j] - 1e-12) { g[j] = cand; changed = true; }
        }
        if (!changed) break;
    }
    for (int i = 0; i < n; ++i)
        if (g[i] >= POS_INF / 2) g[i] = 0.0;
    return g;
}

// Exact (clipped, sequence-dependent) fatigue-adjusted route cost, used to
// validate a concrete candidate route (SA output or extracted B&C integer
// solution): F_0=0, F_{l+1}=max(0, F_l + phi_plus - rho*phi_minus),
// cost_leg = C_leg * (1 + fatigue_rate * F_l).
static double rcost_fatigue_asym(const Input& inp, const std::vector<int>& route, double rho) {
    if (route.empty()) return inp.cm[0][0];
    std::vector<int> seq = {0};
    seq.insert(seq.end(), route.begin(), route.end());
    seq.push_back(0);
    double total = 0.0, F = 0.0;
    for (size_t k = 0; k + 1 < seq.size(); ++k) {
        int a = seq[k], b = seq[k + 1];
        total += inp.cm[a][b] * (1.0 + inp.fatigue_rate * F);
        F = std::max(0.0, F + psi_arc(inp, a, b, rho));
    }
    return total;
}

#include "Highs.h"  // replace #include "glpk.h"

struct LPModel {
    int n = 0;
    Highs* highs = nullptr;
    std::vector<std::vector<int>> x_col;
    std::vector<int> y_col;
    std::vector<std::vector<int>> f_col;   // flow variables (cumulative time) -- legacy, kept for A-B comparison
    std::vector<std::vector<int>> g_col;   // flow variables (cumulative asymmetric fatigue state G_i)
    std::vector<double> Ghat;              // per-node upper bound on G_i, from compute_fatigue_bounds()
    std::vector<double> Glow;              // per-node lower bound on G_i, from compute_fatigue_lower_bounds()
    std::vector<double> col_ub_cache;  // cached upper bounds per column
    std::vector<double> sol_cache;      // cached primal solution after each solve()
    int last_model_status = 0;
    int n_cols_base = 0;
    int n_rows_base = 0;

    LPModel() = default;
    ~LPModel() { if (highs) Highs_destroy(highs); }
    LPModel(const LPModel&) = delete;
    LPModel& operator=(const LPModel&) = delete;

    // ΓöÇΓöÇ helpers ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ

    // Add a column: returns 0-based col index
    int add_col(double lb, double ub, double obj = 0.0) {
        Highs_addCol(highs, obj, lb, ub, 0, nullptr, nullptr);
        col_ub_cache.push_back(ub);
        return static_cast<int>(col_ub_cache.size()) - 1;
    }

    // Add a row: lhs <= sum(coeffs * cols) <= rhs
    // cols and coeffs are 0-indexed, no dummy element
    void add_row(double lhs, double rhs,
                 const std::vector<int>& cols,
                 const std::vector<double>& coeffs) {
        assert(cols.size() == coeffs.size());
        Highs_addRow(highs, lhs, rhs,
                     static_cast<int>(cols.size()),
                     cols.data(), coeffs.data());
    }

    // Fix a column to a value
    void fix_col(int col, double val) {
        Highs_changeColBounds(highs, col, val, val);
        col_ub_cache[col] = val;
    }

    double get_col_ub(int col) const {
        return col_ub_cache[col];
    }



    // ΓöÇΓöÇ build ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ

    void build(const Input& inp) {
        n = static_cast<int>(inp.pts.size());
        highs = (Highs*)Highs_create();
        Highs_setBoolOptionValue(highs, "output_flag", false);
        Highs_setStringOptionValue(highs, "presolve", "on");
        Highs_setStringOptionValue(highs, "solver", "simplex");
        Highs_setStringOptionValue(highs, "simplex_strategy", "1"); // dual simplex

        x_col.assign(n, std::vector<int>(n, -1));
        y_col.resize(n, -1);
        f_col.assign(n, std::vector<int>(n, -1));
        g_col.assign(n, std::vector<int>(n, -1));

        // x[i][j] ΓÇö pre-fix structurally infeasible arcs (base cost + fatigue-aware)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (i == j || !std::isfinite(inp.cm[i][j])) continue;
                bool infeasible = !std::isfinite(inp.cm[j][0]) ||
                                  inp.cm[i][j] + inp.cm[j][0] > inp.bud_raw ||
                                  inp.cm[0][i] + inp.cm[i][j] + inp.cm[j][0] > inp.bud_raw;
                if (!infeasible && inp.fatigue_rate > 0) {
                    // Corrected (clipped, asymmetric) model on the direct
                    // 0->i->j->0 route, replacing the old linear
                    // 1+lambda*t/B pre-filter test.
                    double fat_cost = rcost_fatigue_asym(inp, {i, j}, inp.rho);
                    if (fat_cost > inp.bud_raw) infeasible = true;
                }
                x_col[i][j] = add_col(0.0, infeasible ? 0.0 : 1.0);
            }

        // y[i]
        for (int i = 0; i < n; ++i)
            y_col[i] = add_col(0.0, 1.0, inp.pts[i]);
        fix_col(y_col[0], 1.0);

        // f[i][j] ΓÇö flow variables carrying cumulative time along each arc
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (x_col[i][j] < 0) continue;
                f_col[i][j] = add_col(0.0, inp.bud_raw);
            }

        // Ghat_i / Glow_i bounds (Bellman-Ford, max-plus / min-plus, n-1
        // rounds) -- see compute_fatigue_bounds() / compute_fatigue_lower_bounds()
        // above build().
        Ghat = compute_fatigue_bounds(inp, inp.rho);
        Glow = compute_fatigue_lower_bounds(inp, inp.rho);

        // g[i][j] -- asymmetric fatigue-flow variables, only on surviving
        // arcs. Column bound must include 0 unconditionally
        // (min(Glow[i],0)..max(Ghat[i],0)), not [Glow[i],Ghat[i]] directly:
        // when Glow[i]>0 (a node reachable only via net-uphill prefixes),
        // an unconditional Glow[i] floor forces g_ij>0 on EVERY outgoing
        // arc regardless of x_ij, which via add_fatigue_flow_coupling
        // forces x_ij>0 on all of them simultaneously -- infeasible
        // against flow conservation. See solver_flow_torremocha.cpp.
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (x_col[i][j] < 0) continue;
                g_col[i][j] = add_col(std::min(Glow[i], 0.0), std::max(Ghat[i], 0.0));
            }

        add_flow_conservation();
        add_time_flow_coupling(inp.cm, inp.bud_raw);
        add_time_flow_propagation(inp.cm, inp.bud_raw);
        // NOTE: add_fatigue_budget_flow (old, linear-in-time model) is
        // replaced by the asymmetric-state budget below. Kept defined
        // further down for A/B comparison if needed, just not called.
        add_fatigue_flow_coupling();
        add_fatigue_flow_propagation(inp, inp.rho);
        add_fatigue_budget_asym(inp.cm, inp.bud_raw, inp.fatigue_rate);

        n_rows_base = Highs_getNumRow(highs);
    }

    void clone_from(const LPModel& other,
                    const std::vector<std::pair<int,double>>& fixings) {
        assert(highs == nullptr);
        n               = other.n;
        x_col           = other.x_col;
        y_col           = other.y_col;
        f_col           = other.f_col;
        g_col           = other.g_col;
        Ghat            = other.Ghat;
        Glow            = other.Glow;
        col_ub_cache    = other.col_ub_cache;
        n_rows_base     = other.n_rows_base;
        // added_secs starts empty ΓÇö each cloned node tracks its own cuts

        // Deep-copy the HiGHS model
        highs = (Highs*)Highs_create();
        Highs_setBoolOptionValue(highs, "output_flag", false);
        Highs_setStringOptionValue(highs, "presolve", "on");
        Highs_setStringOptionValue(highs, "solver", "simplex");

        // Deep-copy via Highs_passLp ΓÇö no temp file, faster, thread-safe
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

    // ΓöÇΓöÇ constraints (flow-based) ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ

    void add_flow_conservation() {
        for (int i = 0; i < n; ++i) {
            std::vector<int> cols; std::vector<double> coeffs;
            for (int j = 0; j < n; ++j)
                if (x_col[j][i] >= 0) { cols.push_back(x_col[j][i]); coeffs.push_back(1.0); }
            cols.push_back(y_col[i]); coeffs.push_back(-1.0);
            add_row(0.0, 0.0, cols, coeffs);
        }
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
                add_row(-1e30, 0.0,
                        {f_col[i][j], x_col[i][j]},
                        {1.0,         -bud_raw});
            }
    }

    void add_time_flow_propagation(const std::vector<std::vector<double>>& cm, double bud_raw) {
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
        for (int j = 0; j < n; ++j) {
            if (f_col[0][j] >= 0)
                fix_col(f_col[0][j], 0.0);
        }
    }

    void add_fatigue_budget_flow(const std::vector<std::vector<double>>& cm,
                                 double bud_raw, double fatigue_rate) {
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

    // -- Asymmetric fatigue-state model (replaces the linear one above) --

    void add_fatigue_flow_coupling() {
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (g_col[i][j] < 0) continue;
                add_row(-1e30, 0.0,
                        {g_col[i][j], x_col[i][j]},
                        {1.0,         -Ghat[i]});
                add_row(0.0, 1e30,
                        {g_col[i][j], x_col[i][j]},
                        {1.0,         -Glow[i]});
            }
    }

    void add_fatigue_flow_propagation(const Input& inp, double rho) {
        for (int j = 1; j < n; ++j) {
            std::vector<int> cols; std::vector<double> coeffs;
            for (int i = 0; i < n; ++i) {
                if (g_col[i][j] >= 0) { cols.push_back(g_col[i][j]); coeffs.push_back(1.0); }
                if (x_col[i][j] >= 0) { cols.push_back(x_col[i][j]); coeffs.push_back(psi_arc(inp, i, j, rho)); }
            }
            for (int k = 0; k < n; ++k) {
                if (g_col[j][k] >= 0) { cols.push_back(g_col[j][k]); coeffs.push_back(-1.0); }
            }
            if (!cols.empty())
                add_row(0.0, 0.0, cols, coeffs);
        }
        for (int j = 0; j < n; ++j) {
            if (g_col[0][j] >= 0)
                fix_col(g_col[0][j], 0.0);
        }
    }

    void add_fatigue_budget_asym(const std::vector<std::vector<double>>& cm,
                                  double bud_raw, double fatigue_rate) {
        std::vector<int> cols; std::vector<double> coeffs;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (x_col[i][j] >= 0) {
                    cols.push_back(x_col[i][j]);
                    coeffs.push_back(cm[i][j]);
                }
                if (g_col[i][j] >= 0) {
                    cols.push_back(g_col[i][j]);
                    coeffs.push_back(fatigue_rate * cm[i][j]);
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

    // ΓöÇΓöÇ solve ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ

    bool solve() {
        Highs_changeObjectiveSense(highs, -1);  // -1 = kHighsObjSenseMaximize
        int status = Highs_run(highs);
        last_model_status = Highs_getModelStatus(highs);
        bool ok = status == 0 && last_model_status == 7;
        if (!ok && solve_failed_unresolved()) {
            // A warm-started re-solve after adding cuts reuses the previous
            // basis and skips presolve/rescaling -- on this model's
            // coefficient range (matrix entries span ~1e8, mixing raw
            // distances/costs with 0/1 columns) that can leave HiGHS unable
            // to certify Optimal/Infeasible even when the point is
            // numerically fine (tiny residuals). Retry once from a cleared
            // solver state so presolve+scaling run fresh, matching what the
            // (always-clean) root solve already does successfully.
            Highs_clearSolver(highs);
            status = Highs_run(highs);
            last_model_status = Highs_getModelStatus(highs);
            ok = status == 0 && last_model_status == 7;
        }
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

    // kInfeasible (8) is the only status where dropping the node is sound.
    // Any other non-Optimal status means the LP was NOT resolved -- pruning
    // it anyway silently truncates the search like the depth-limit bug.
    bool solve_failed_unresolved() const { return last_model_status != 7 && last_model_status != 8; }

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


// ΓöÇΓöÇ Subtour detection (Kosaraju SCC for asymmetric) ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ

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

// ΓöÇΓöÇ Depot-unreachable detection (asymmetric only) ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ
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

// ΓöÇΓöÇ Lifted cover cuts on the budget knapsack ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ

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

// ΓöÇΓöÇ Route extraction & validation ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ

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
    // Corrected model: clipped, asymmetric, sequence-dependent fatigue.
    return rcost_fatigue_asym(inp, route, inp.rho) <= inp.bud_raw;
}

// ΓöÇΓöÇ Greedy heuristic ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ

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

// ΓöÇΓöÇ Simulated Annealing ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ

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
        if (rcost_fatigue_asym(inp, new_route, inp.rho) > inp.bud_raw) continue;

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
           rcost_fatigue_asym(inp, best_route, inp.rho) > inp.bud_raw) {
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

// ΓöÇΓöÇ Branch-and-Cut Node ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ

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
    int max_depth = 200;  // raised from 15: x_ij branching (see integrality-check fix) needs far more decisions to fully resolve than y-only branching did
    double time_limit_s = 900.0;
    // See cpp/solver_flow_torremocha.cpp: a node dropped purely for
    // exceeding max_depth is unresolved, not proven suboptimal --
    // proved_optimal must not claim completeness once this fires.
    bool depth_limit_hit = false;
    // See LPModel::solve_failed_unresolved(): a node dropped because HiGHS
    // failed to reach kOptimal/kInfeasible is unresolved, not proven
    // suboptimal -- same completeness requirement as depth_limit_hit above.
    bool lp_solve_failure_hit = false;
    // The outer B&C loop only checks time_limit_s between nodes -- a single
    // node's cut-generation loop had no time check at all, so one slow node
    // could run past the budget indefinitely (observed: 13+ hours stuck on
    // one root node on a n=150 instance, in the sibling benchmark flow
    // solver). This flag lets process_node bail out of its own cut loop
    // once the overall time budget is spent, same completeness caveat as
    // the other two flags.
    bool time_limit_hit = false;
    std::chrono::steady_clock::time_point t_start;

    bool proved_optimal = false;
    double best_ub = std::numeric_limits<double>::infinity();
    int nodes_explored = 0;
    int nodes_unexplored = 0;

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
        if (!root_ok && root.solve_failed_unresolved()) lp_solve_failure_hit = true;
        root_node.ub = root_ok ? root.obj() : -std::numeric_limits<double>::infinity();
        std::cerr << "Root UB=" << root_node.ub << "  best_pts=" << best_pts << "\n";
        if (root_node.ub > best_pts) node_stack.push(std::move(root_node));

        int nodes = 0;
        t_start = std::chrono::steady_clock::now();
        while (!node_stack.empty() && nodes++ < 10000) {
            double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
            if (elapsed > time_limit_s) { std::cerr << "Time limit reached\n"; break; }
            BNCNode node = std::move(node_stack.top());
            node_stack.pop();
            process_node(std::move(node), node_stack);
        }
        proved_optimal = node_stack.empty() && nodes < 10000 && !depth_limit_hit && !lp_solve_failure_hit && !time_limit_hit;
        nodes_explored = nodes;
        nodes_unexplored = static_cast<int>(node_stack.size());
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
        std::cerr << "Explored: " << nodes_explored << ", Unexplored: " << nodes_unexplored
                  << ", best: " << best_pts << " pts, UB: " << best_ub
                  << ", gap: " << gap_pct << "%\n";
    }

    void process_node(BNCNode node, std::stack<BNCNode>& node_stack) {
        if (node.ub <= best_pts + 1e-6) return;
        if (node.fixings.size() > max_depth) { depth_limit_hit = true; return; }

        // Clone LP and apply fixings
        LPModel lp;
        lp.clone_from(root, node.fixings);

        // Cut loop -- accumulate cuts, do NOT delete them between iterations
        for (int cut_iter = 0; cut_iter < max_cuts; ++cut_iter) {
            double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
            if (elapsed > time_limit_s) { time_limit_hit = true; return; }
            if (!lp.solve()) { if (lp.solve_failed_unresolved()) lp_solve_failure_hit = true; return; }
            double lp_ub = lp.obj();
            std::cerr << "LP UB: " << lp_ub << " (gap: " << (lp_ub-best_pts) << ")\n";
            if (lp_ub <= best_pts + 1e-6) return;

            auto subtours = find_subtours(lp);
            auto unreachable = find_depot_unreachable(lp);
            subtours.insert(subtours.end(), unreachable.begin(), unreachable.end());

            int cover_cuts = find_and_add_cover_cuts(lp, inp);

            if (subtours.empty() && cover_cuts == 0) break;

            for (const auto& S : subtours)
                lp.add_sec(S);
        }

        if (!lp.solve()) { if (lp.solve_failed_unresolved()) lp_solve_failure_hit = true; return; }
        double lp_ub = lp.obj();
        if (lp_ub <= best_pts + 1e-6) return;

        // Check integrality (simple: all vars 0/1 within tol)
        bool integer_sol = true;
        int branch_col = -1;
        double max_frac = 0.0;
        for (int i = 1; i < root.n; ++i) {  // Branch on y[i]
            double v = lp.prim(root.y_col[i]);
            double frac = std::min(v, 1.0 - v);
            if (frac > 1e-5) {
                integer_sol = false;
                if (frac > max_frac) {
                    max_frac = frac;
                    branch_col = root.y_col[i];
                }
            }
        }

        // y integral alone does not certify an integer-feasible route: x_ij
        // can still be fractional even when every node's visit status is
        // fixed (e.g. two fractional sub-routings both consistent with the
        // same y-degree constraints, not yet separated by any SEC/
        // connectivity cut found so far). extract_route()'s x>0.5
        // threshold-and-follow-successor heuristic would silently truncate
        // or misroute in that case, undercounting the true score for a
        // genuinely achievable y-pattern -- so a leaf is only accepted as
        // resolved once x is ALSO fully integral; otherwise branch on the
        // most fractional x_ij, same as for y above.
        if (integer_sol) {
            for (int i = 0; i < root.n && integer_sol; ++i) {
                for (int j = 0; j < root.n; ++j) {
                    if (root.x_col[i][j] < 0) continue;
                    double v = lp.prim(root.x_col[i][j]);
                    double frac = std::min(v, 1.0 - v);
                    if (frac > 1e-5) {
                        integer_sol = false;
                        if (frac > max_frac) {
                            max_frac = frac;
                            branch_col = root.x_col[i][j];
                        }
                    }
                }
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
        } else if (branch_col >= 0) {
            // Branch
            BNCNode node0 = node, node1 = node;
            node0.fixings.emplace_back(branch_col, 0.0);
            node0.ub = lp_ub;
            node1.fixings.emplace_back(branch_col, 1.0);
            node1.ub = lp_ub;
            // Push in LIFO order (DFS)
            node_stack.push(std::move(node0));
            node_stack.push(std::move(node1));
        }
    }
};

// ΓöÇΓöÇ Main ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ

static void run_map(const std::string& in_path, const std::string& out_path,
                     double rho_override = -1.0, double lambda_override = -1.0) {
    std::cerr << "\n=== " << in_path << " ===\n";
    Input inp = parse_input(read_file(in_path));
    if (rho_override >= 0.0) inp.rho = rho_override;
    if (lambda_override >= 0.0) inp.fatigue_rate = lambda_override;
    std::cerr << "  rho = " << inp.rho << ", lambda = " << inp.fatigue_rate << "\n";

    auto t_sa = std::chrono::steady_clock::now();
    auto sa_route = solve_sa_iterated(inp);
    double sa_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_sa).count();
    double sa_pts     = rpts(inp.pts, sa_route);
    double sa_base    = rcost(inp.cm, sa_route);
    double sa_fatigue        = rcost_fatigue_asym(inp, sa_route, inp.rho);
    double sa_fatigue_legacy = rcost_fatigue(inp.cm, sa_route, inp.bud_raw, inp.fatigue_rate);
    std::cerr << "SA: " << sa_pts << " pts (" << sa_route.size() << " nodes) in " << sa_elapsed << "s\n";

    auto t_bnc = std::chrono::steady_clock::now();
    Solver solver(inp);
    solver.solve(sa_pts, sa_route);
    double bnc_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_bnc).count();
    double bnc_base    = rcost(inp.cm, solver.best_route);
    double bnc_fatigue        = rcost_fatigue_asym(inp, solver.best_route, inp.rho);
    double bnc_fatigue_legacy = rcost_fatigue(inp.cm, solver.best_route, inp.bud_raw, inp.fatigue_rate);
    std::cerr << "B&C(SA): " << solver.best_pts << " pts (" << solver.best_route.size() << " nodes) in " << bnc_elapsed << "s\n";

    std::ofstream out(out_path);
    out << "{\n";
    out << "  \"rho\": " << inp.rho << ",\n";
    out << "  \"lambda\": " << inp.fatigue_rate << ",\n";
    out << "  \"sa\": {\"pts\": " << sa_pts << ", \"nodes\": " << sa_route.size()
        << ", \"elapsed_s\": " << sa_elapsed << ", \"base_cost\": " << sa_base
        << ", \"fatigue_cost\": " << sa_fatigue
        << ", \"fatigue_cost_legacy\": " << sa_fatigue_legacy << ", \"route\": [";
    for (size_t i = 0; i < sa_route.size(); ++i) { if (i) out << ", "; out << sa_route[i]; }
    out << "]},\n";
    out << "  \"bnc_sa\": {\"pts\": " << solver.best_pts << ", \"nodes\": " << solver.best_route.size()
        << ", \"elapsed_s\": " << bnc_elapsed
        << ", \"proved_optimal\": " << (solver.proved_optimal ? "true" : "false")
        << ", \"best_ub\": " << solver.best_ub
        << ", \"gap_pct\": " << (solver.best_pts > 0 ? 100.0 * (solver.best_ub - solver.best_pts) / solver.best_pts : 0.0)
        << ", \"nodes_explored\": " << solver.nodes_explored
        << ", \"nodes_unexplored\": " << solver.nodes_unexplored
        << ", \"base_cost\": " << bnc_base
        << ", \"fatigue_cost\": " << bnc_fatigue
        << ", \"fatigue_cost_legacy\": " << bnc_fatigue_legacy << ", \"route\": [";
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

int main(int argc, char** argv) {
    const std::vector<std::pair<std::string,std::string>> maps = {
        {"op_input_standard.json",     "op_output_standard.json"},
        {"op_input_clustered.json",     "op_output_clustered.json"},
        {"op_input_ring.json",          "op_output_ring.json"},
        {"op_input_path_biased.json",   "op_output_path_biased.json"},
        {"op_input_elev_biased.json",   "op_output_elev_biased.json"},
        {"op_input_sparse_far.json",    "op_output_sparse_far.json"},
        {"op_input_mixed_density.json", "op_output_mixed_density.json"},

    };

    // Track whether B&C hit the time limit per map via a flag set in run_map
    // We re-read the output JSON to extract results for the summary table.
    // run_map already writes elapsed_s to the JSON, so we parse that.
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
            std::string sa_blk  = js.substr(0, bnc_pos);
            std::string bnc_blk = js.substr(bnc_pos);
            MapResult r;
            // strip "op_input_" prefix and ".json" suffix for display
            r.name      = in.substr(9, in.size() - 14);
            r.sa_pts    = json_val(sa_blk,  "pts");
            r.sa_nodes  = static_cast<int>(json_val(sa_blk,  "nodes"));
            r.sa_s      = json_val(sa_blk,  "elapsed_s");
            r.bnc_pts   = json_val(bnc_blk, "pts");
            r.bnc_nodes = static_cast<int>(json_val(bnc_blk, "nodes"));
            r.bnc_s     = json_val(bnc_blk, "elapsed_s");
            r.bnc_optimal = (bnc_blk.find("\"proved_optimal\": true") != std::string::npos);
            r.gap_pct = json_val(bnc_blk, "gap_pct");
            results.push_back(r);
        } catch (...) {}
    }

    // ΓöÇΓöÇ Summary table ΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇΓöÇ
    std::cout << "\n";
    std::cout << "+----------------------+---------------------------+----------------------------------------------+\n";
    std::cout << "| Map                  | SA                        | B&C(SA)                                      |\n";
    std::cout << "|                      |  pts  nodes     time      |  pts  nodes     time    optimal?    gap       |\n";
    std::cout << "+----------------------+---------------------------+----------------------------------------------+\n";
    for (const auto& r : results) {
        std::cout << "| " << std::left  << std::setw(20) << r.name
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
    std::cout << "+----------------------+---------------------------+----------------------------------------------+\n";

    return 0;
}
