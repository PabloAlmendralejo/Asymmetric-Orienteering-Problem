/*
 * Orienteering Problem — Branch and Cut with GLPK LP relaxation + SEC cuts
 * Input:  op_input.json  { "cm": [[...]], "pts": [...], "bud_eff": float, "bud_raw": float, "fatigue_rate": float }
 * Output: op_output.json { "route": [...], "pts": float, "base_cost": float, "fatigue_cost": float, "elapsed_s": float }
 *
 * Compile (MSVC x64):
 *   cl /O2 /EHsc /std:c++17 op_bnc.cpp /I"C:\Users\borrepa\Downloads\winglpk-4.65\glpk-4.65\src" ^
 *      /link "C:\Users\borrepa\Downloads\winglpk-4.65\glpk-4.65\w64\glpk_4_65.lib" /Fe:op_bnc.exe
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>
#include <chrono>
#include <limits>
#include <queue>
#include <functional>

#include "glpk.h"

// ── JSON parser (reused from op_bnb.cpp) ──────────────────────────────────

static std::string read_file(const std::string& path) {
    std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open " << path << "\n"; std::exit(1); }
    std::ostringstream ss; ss << f.rdbuf(); return ss.str();
}
static size_t skip_ws(const std::string& s, size_t i) {
    while (i < s.size() && (s[i]==' '||s[i]=='\t'||s[i]=='\n'||s[i]=='\r')) ++i;
    return i;
}
static double parse_number(const std::string& s, size_t& i) {
    size_t start = i;
    if (i < s.size() && (s[i]=='-'||s[i]=='+')) ++i;
    while (i < s.size() && (std::isdigit(s[i])||s[i]=='.'||s[i]=='e'||s[i]=='E'||s[i]=='+'||s[i]=='-')) ++i;
    return std::stod(s.substr(start, i - start));
}
static std::vector<double> parse_array1d(const std::string& s, size_t& i) {
    std::vector<double> v;
    i = skip_ws(s, i); assert(s[i]=='['); ++i;
    while (true) {
        i = skip_ws(s, i);
        if (s[i]==']') { ++i; break; }
        if (s[i]==',') { ++i; continue; }
        v.push_back(parse_number(s, i));
    }
    return v;
}
static std::vector<std::vector<double>> parse_array2d(const std::string& s, size_t& i) {
    std::vector<std::vector<double>> m;
    i = skip_ws(s, i); assert(s[i]=='['); ++i;
    while (true) {
        i = skip_ws(s, i);
        if (s[i]==']') { ++i; break; }
        if (s[i]==',') { ++i; continue; }
        if (s[i]=='[') m.push_back(parse_array1d(s, i));
    }
    return m;
}

struct Input {
    std::vector<std::vector<double>> cm;
    std::vector<double> pts;
    double bud_eff, bud_raw, fatigue_rate;
};

static Input parse_input(const std::string& json) {
    Input inp; size_t i = 0;
    auto find_key = [&](const std::string& key) {
        size_t pos = json.find("\"" + key + "\"", i);
        assert(pos != std::string::npos);
        i = pos + key.size() + 2;
        i = skip_ws(json, i); assert(json[i]==':'); ++i;
        i = skip_ws(json, i);
    };
    find_key("cm");           inp.cm           = parse_array2d(json, i);
    find_key("pts");          inp.pts          = parse_array1d(json, i);
    find_key("bud_eff");      inp.bud_eff      = parse_number(json, i);
    find_key("bud_raw");      inp.bud_raw      = parse_number(json, i);
    find_key("fatigue_rate"); inp.fatigue_rate = parse_number(json, i);
    return inp;
}

// ── Cost helpers ───────────────────────────────────────────────────────────

static double rcost(const std::vector<std::vector<double>>& cm,
                    const std::vector<int>& route) {
    if (route.empty()) return 0.0;
    double c = cm[0][route[0]];
    for (size_t i = 0; i+1 < route.size(); ++i) c += cm[route[i]][route[i+1]];
    return c + cm[route.back()][0];
}

static double rcost_fatigue(const std::vector<std::vector<double>>& cm,
                             const std::vector<int>& route,
                             double bud_raw, double fatigue_rate) {
    if (route.empty()) return 0.0;
    double total = 0.0, elapsed = 0.0;
    std::vector<int> seq = {0};
    for (int x : route) seq.push_back(x);
    seq.push_back(0);
    for (size_t i = 0; i+1 < seq.size(); ++i) {
        double leg = cm[seq[i]][seq[i+1]];
        total   += leg * (1.0 + fatigue_rate * (elapsed / std::max(bud_raw, 1.0)));
        elapsed += leg;
    }
    return total;
}

// ── Arc index helpers ──────────────────────────────────────────────────────
// Arcs: x[i][j] for i!=j, i,j in 0..n-1. Linearised as arc_id = i*(n-1) + (j < i ? j : j-1)
// We only create arcs where cm[i][j] is finite.

struct Arc { int i, j, col; }; // col = 1-based GLPK column index

// ── LP formulation ─────────────────────────────────────────────────────────
//
// Variables:
//   x[i][j]  ∈ [0,1]  arc i→j used          (n*(n-1) arcs, only finite ones)
//   y[i]     ∈ [0,1]  node i visited         (n nodes, y[0] fixed = 1)
//
// Objective: max  sum_i  pts[i] * y[i]
//
// Constraints:
//   (1) Flow balance:  sum_j x[j][i] = sum_j x[i][j] = y[i]   for all i
//   (2) Budget:        sum_{i,j} cm[i][j] * x[i][j] <= bud_eff
//   (3) y[0] = 1  (depot always visited)
//   (4) SECs added dynamically: sum_{i,j in S} x[i][j] <= |S|-1  for subtour S
//
// MTZ is avoided — we use pure SEC cuts found by connected-component detection
// on the fractional LP solution graph.

struct LPModel {
    int n;
    glp_prob* lp = nullptr;

    // column indices (1-based)
    std::vector<std::vector<int>> x_col; // x_col[i][j], 0 if arc doesn't exist
    std::vector<int>              y_col; // y_col[i]

    int n_arcs = 0;
    int n_rows_base = 0; // rows before any SECs

    void build(const std::vector<std::vector<double>>& cm,
               const std::vector<double>& pts,
               double bud_eff) {
        n = (int)pts.size();
        lp = glp_create_prob();
        glp_set_obj_dir(lp, GLP_MAX);
        glp_term_out(GLP_OFF); // silence GLPK output

        x_col.assign(n, std::vector<int>(n, 0));
        y_col.resize(n);

        // add arc columns x[i][j]
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (i != j && std::isfinite(cm[i][j])) {
                    int c = glp_add_cols(lp, 1);
                    glp_set_col_bnds(lp, c, GLP_DB, 0.0, 1.0);
                    glp_set_obj_coef(lp, c, 0.0); // pts come from y
                    x_col[i][j] = c;
                    ++n_arcs;
                }

        // add node columns y[i]
        for (int i = 0; i < n; ++i) {
            int c = glp_add_cols(lp, 1);
            glp_set_col_bnds(lp, c, GLP_DB, 0.0, 1.0);
            glp_set_obj_coef(lp, c, pts[i]);
            y_col[i] = c;
        }
        // depot y[0] = 1
        glp_set_col_bnds(lp, y_col[0], GLP_FX, 1.0, 1.0);

        // flow-in = y[i]: sum_j x[j][i] = y[i]
        for (int i = 0; i < n; ++i) {
            int r = glp_add_rows(lp, 1);
            glp_set_row_bnds(lp, r, GLP_FX, 0.0, 0.0);
            std::vector<int>    idx; idx.push_back(0);
            std::vector<double> val; val.push_back(0.0);
            for (int j = 0; j < n; ++j)
                if (x_col[j][i]) { idx.push_back(x_col[j][i]); val.push_back(1.0); }
            idx.push_back(y_col[i]); val.push_back(-1.0);
            glp_set_mat_row(lp, r, (int)idx.size()-1, idx.data(), val.data());
        }

        // flow-out = y[i]: sum_j x[i][j] = y[i]
        for (int i = 0; i < n; ++i) {
            int r = glp_add_rows(lp, 1);
            glp_set_row_bnds(lp, r, GLP_FX, 0.0, 0.0);
            std::vector<int>    idx; idx.push_back(0);
            std::vector<double> val; val.push_back(0.0);
            for (int j = 0; j < n; ++j)
                if (x_col[i][j]) { idx.push_back(x_col[i][j]); val.push_back(1.0); }
            idx.push_back(y_col[i]); val.push_back(-1.0);
            glp_set_mat_row(lp, r, (int)idx.size()-1, idx.data(), val.data());
        }

        // budget constraint
        {
            int r = glp_add_rows(lp, 1);
            glp_set_row_bnds(lp, r, GLP_UP, 0.0, bud_eff);
            std::vector<int>    idx; idx.push_back(0);
            std::vector<double> val; val.push_back(0.0);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    if (x_col[i][j]) { idx.push_back(x_col[i][j]); val.push_back(cm[i][j]); }
            glp_set_mat_row(lp, r, (int)idx.size()-1, idx.data(), val.data());
        }

        n_rows_base = glp_get_num_rows(lp);
    }

    // fix a column to a value (for branching)
    void fix_col(int col, double val) {
        glp_set_col_bnds(lp, col, GLP_FX, val, val);
    }
    void free_col(int col, double lo, double hi) {
        glp_set_col_bnds(lp, col, GLP_DB, lo, hi);
    }

    // solve LP relaxation, return objective or -inf if infeasible
    double solve_lp() {
        glp_smcp parm; glp_init_smcp(&parm);
        parm.msg_lev = GLP_MSG_OFF;
        parm.meth    = GLP_DUAL;
        int ret = glp_simplex(lp, &parm);
        if (ret != 0) return -std::numeric_limits<double>::infinity();
        int stat = glp_get_status(lp);
        if (stat == GLP_OPT) return glp_get_obj_val(lp);
        return -std::numeric_limits<double>::infinity();
    }

    // add SEC for a subset S (0-based node indices): sum_{i,j in S} x[i][j] <= |S|-1
    void add_sec(const std::vector<int>& S) {
        int r = glp_add_rows(lp, 1);
        glp_set_row_bnds(lp, r, GLP_UP, 0.0, (double)(S.size() - 1));
        std::vector<int>    idx; idx.push_back(0);
        std::vector<double> val; val.push_back(0.0);
        for (int i : S)
            for (int j : S)
                if (i != j && x_col[i][j]) { idx.push_back(x_col[i][j]); val.push_back(1.0); }
        glp_set_mat_row(lp, r, (int)idx.size()-1, idx.data(), val.data());
    }

    // remove all SEC rows added after base (reset for new B&C node)
    void remove_secs() {
        int cur = glp_get_num_rows(lp);
        if (cur <= n_rows_base) return;
        std::vector<int> rows(cur - n_rows_base + 1);
        for (int i = 1; i <= cur - n_rows_base; ++i) rows[i] = n_rows_base + i;
        glp_del_rows(lp, cur - n_rows_base, rows.data());
    }

    ~LPModel() { if (lp) glp_delete_prob(lp); }
};

// ── Subtour detection ──────────────────────────────────────────────────────
// Build directed graph from x[i][j] > eps, find connected components
// excluding depot (node 0). Any component not containing 0 is a subtour.

static std::vector<std::vector<int>> find_subtours(
        const LPModel& m, double eps = 0.5) {
    int n = m.n;
    // build adjacency from fractional solution
    std::vector<std::vector<int>> adj(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j && m.x_col[i][j]) {
                double v = glp_get_col_prim(m.lp, m.x_col[i][j]);
                if (v > eps) adj[i].push_back(j);
            }

    std::vector<int> comp(n, -1);
    int nc = 0;
    for (int s = 0; s < n; ++s) {
        if (comp[s] >= 0) continue;
        std::queue<int> q; q.push(s); comp[s] = nc;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int v : adj[u]) if (comp[v] < 0) { comp[v] = nc; q.push(v); }
            // also traverse reverse for undirected connectivity
            for (int v = 0; v < n; ++v)
                if (comp[v] < 0 && m.x_col[v][u]) {
                    double val = glp_get_col_prim(m.lp, m.x_col[v][u]);
                    if (val > eps) { comp[v] = nc; q.push(v); }
                }
        }
        ++nc;
    }

    // collect components not containing depot (node 0)
    int depot_comp = comp[0];
    std::vector<std::vector<int>> subtours;
    for (int c = 0; c < nc; ++c) {
        if (c == depot_comp) continue;
        std::vector<int> S;
        for (int i = 0; i < n; ++i) if (comp[i] == c) S.push_back(i);
        if (S.size() >= 2) subtours.push_back(S);
    }
    return subtours;
}

// ── Extract integer route from LP solution ─────────────────────────────────

static std::vector<int> extract_route(const LPModel& m, double eps = 0.5) {
    int n = m.n;
    // build successor map
    std::vector<int> succ(n, -1);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j && m.x_col[i][j]) {
                double v = glp_get_col_prim(m.lp, m.x_col[i][j]);
                if (v > eps) succ[i] = j;
            }
    std::vector<int> route;
    int cur = succ[0];
    while (cur > 0 && cur != -1) {
        route.push_back(cur);
        cur = succ[cur];
        if ((int)route.size() > n) break; // cycle guard
    }
    return route;
}

// ── Check if LP solution is integer ───────────────────────────────────────

static bool is_integer(const LPModel& m, double eps = 1e-4) {
    int ncols = glp_get_num_cols(m.lp);
    for (int c = 1; c <= ncols; ++c) {
        double v = glp_get_col_prim(m.lp, c);
        if (v > eps && v < 1.0 - eps) return false;
    }
    return true;
}

// ── Most fractional variable for branching ────────────────────────────────

static int most_fractional_arc(const LPModel& m) {
    int n = m.n;
    int best_col = -1;
    double best_frac = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j && m.x_col[i][j]) {
                double v = glp_get_col_prim(m.lp, m.x_col[i][j]);
                double frac = std::min(v, 1.0 - v);
                if (frac > best_frac) { best_frac = frac; best_col = m.x_col[i][j]; }
            }
    return best_col;
}

// ── Greedy warm start ──────────────────────────────────────────────────────
// Repeatedly pick best pts/cost node reachable from current position.

static std::vector<int> greedy_route(const Input& inp) {
    int n = (int)inp.pts.size();
    std::vector<bool> visited(n, false);
    visited[0] = true;
    std::vector<int> route;
    double cost = 0.0, elapsed = 0.0, fatigue = 0.0;
    int cur = 0;
    while (true) {
        int best_j = -1; double best_ratio = -1.0;
        for (int j = 1; j < n; ++j) {
            if (visited[j]) continue;
            double go = inp.cm[cur][j], back = inp.cm[j][0];
            if (!std::isfinite(go) || !std::isfinite(back)) continue;
            if (cost + go + back > inp.bud_eff) continue;
            double fm  = 1.0 + inp.fatigue_rate * (elapsed / std::max(inp.bud_raw, 1.0));
            double fat_go   = fatigue + go * fm;
            double ret_fm   = 1.0 + inp.fatigue_rate * ((elapsed + go) / std::max(inp.bud_raw, 1.0));
            double fat_back = fat_go + back * ret_fm;
            if (fat_back > inp.bud_raw) continue;
            double ratio = inp.pts[j] / std::max(go, 1e-9);
            if (ratio > best_ratio) { best_ratio = ratio; best_j = j; }
        }
        if (best_j < 0) break;
        double go = inp.cm[cur][best_j];
        double fm = 1.0 + inp.fatigue_rate * (elapsed / std::max(inp.bud_raw, 1.0));
        fatigue += go * fm;
        elapsed += go;
        cost    += go;
        visited[best_j] = true;
        route.push_back(best_j);
        cur = best_j;
    }
    return route;
}

// ── Branch and Cut ─────────────────────────────────────────────────────────
// Each B&C node clones the root LP, applies fixings, runs SEC cut loop,
// then branches on the most fractional y[i] (node variable).

struct BnC {
    const Input& inp;
    LPModel& root;
    double best_pts = 0.0;
    std::vector<int> best_route;

    struct Fix { int col; double val; };

    void solve() {
        // warm start with greedy
        auto gr = greedy_route(inp);
        double gr_pts = 0.0;
        for (int v : gr) gr_pts += inp.pts[v];
        if (gr_pts > best_pts) { best_pts = gr_pts; best_route = gr; }
        std::cerr << "  greedy warm start: " << best_pts << "pts\n";
        bnc_node({});
    }

    void bnc_node(const std::vector<Fix>& fixings) {
        glp_prob* lp = glp_create_prob();
        glp_copy_prob(lp, root.lp, GLP_ON);
        glp_term_out(GLP_OFF);

        for (auto& f : fixings)
            glp_set_col_bnds(lp, f.col, GLP_FX, f.val, f.val);

        // SEC cut loop
        double lp_obj = -1.0;
        for (int cut_iter = 0; cut_iter < 50; ++cut_iter) {
            glp_smcp parm; glp_init_smcp(&parm);
            parm.msg_lev = GLP_MSG_OFF;
            parm.meth    = GLP_DUAL;
            if (glp_simplex(lp, &parm) != 0 || glp_get_status(lp) != GLP_OPT)
                { lp_obj = -1.0; break; }
            lp_obj = glp_get_obj_val(lp);
            if (lp_obj <= best_pts + 1e-6) break;

            LPModel tmp; tmp.lp = lp; tmp.n = root.n;
            tmp.x_col = root.x_col; tmp.y_col = root.y_col;
            auto subtours = find_subtours(tmp, 0.5);
            tmp.lp = nullptr;
            if (subtours.empty()) break;
            for (auto& S : subtours) {
                int r = glp_add_rows(lp, 1);
                glp_set_row_bnds(lp, r, GLP_UP, 0.0, (double)(S.size()-1));
                std::vector<int> idx = {0}; std::vector<double> val = {0.0};
                for (int i : S) for (int j : S)
                    if (i != j && root.x_col[i][j])
                        { idx.push_back(root.x_col[i][j]); val.push_back(1.0); }
                glp_set_mat_row(lp, r, (int)idx.size()-1, idx.data(), val.data());
            }
        }

        if (lp_obj <= best_pts + 1e-6) { glp_delete_prob(lp); return; }

        // check integrality — branch on most fractional y[i]
        int branch_col = -1; double best_frac = 0.0;
        bool integer = true;
        for (int i = 1; i < root.n; ++i) {  // skip depot y[0]
            double v = glp_get_col_prim(lp, root.y_col[i]);
            double frac = std::min(v, 1.0 - v);
            if (frac > 1e-4) {
                integer = false;
                if (frac > best_frac) { best_frac = frac; branch_col = root.y_col[i]; }
            }
        }
        // also check arc variables for integrality
        if (integer) {
            for (int i = 0; i < root.n && integer; ++i)
                for (int j = 0; j < root.n && integer; ++j)
                    if (i != j && root.x_col[i][j]) {
                        double v = glp_get_col_prim(lp, root.x_col[i][j]);
                        if (v > 1e-4 && v < 1.0 - 1e-4) integer = false;
                    }
        }

        if (integer) {
            std::vector<int> succ(root.n, -1);
            for (int i = 0; i < root.n; ++i)
                for (int j = 0; j < root.n; ++j)
                    if (i != j && root.x_col[i][j])
                        if (glp_get_col_prim(lp, root.x_col[i][j]) > 0.5) succ[i] = j;
            std::vector<int> route;
            int cur = succ[0];
            while (cur > 0 && cur != -1 && (int)route.size() <= root.n)
                { route.push_back(cur); cur = succ[cur]; }
            double fc = rcost_fatigue(inp.cm, route, inp.bud_raw, inp.fatigue_rate);
            double pts_got = 0.0;
            for (int v : route) pts_got += inp.pts[v];
            if (fc <= inp.bud_raw && pts_got > best_pts)
                { best_pts = pts_got; best_route = route; }
        } else if (branch_col >= 0) {
            glp_delete_prob(lp);
            auto fix1 = fixings; fix1.push_back({branch_col, 1.0});
            bnc_node(fix1);
            auto fix0 = fixings; fix0.push_back({branch_col, 0.0});
            bnc_node(fix0);
            return;
        }

        glp_delete_prob(lp);
    }
};

// ── Main ───────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {
    std::string in_path  = (argc > 1) ? argv[1] : "op_input.json";
    std::string out_path = (argc > 2) ? argv[2] : "op_output.json";

    auto t0 = std::chrono::steady_clock::now();

    Input inp = parse_input(read_file(in_path));

    LPModel m;
    m.build(inp.cm, inp.pts, inp.bud_eff);

    BnC bnc{inp, m};
    bnc.solve();

    auto t1 = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();

    double base = rcost(inp.cm, bnc.best_route);
    double fc   = rcost_fatigue(inp.cm, bnc.best_route, inp.bud_raw, inp.fatigue_rate);

    std::cerr << "BnC C++: " << bnc.best_pts << "pts, "
              << bnc.best_route.size() << " controls, "
              << elapsed << "s\n";

    std::ofstream out(out_path);
    out << "{\n";
    out << "  \"route\": [";
    for (size_t i = 0; i < bnc.best_route.size(); ++i) { if (i) out << ", "; out << bnc.best_route[i]; }
    out << "],\n";
    out << "  \"pts\": "          << bnc.best_pts << ",\n";
    out << "  \"base_cost\": "    << base         << ",\n";
    out << "  \"fatigue_cost\": " << fc           << ",\n";
    out << "  \"elapsed_s\": "    << elapsed      << "\n";
    out << "}\n";

    return 0;
}
