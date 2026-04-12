
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

#include "glpk.h"

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

// ── LP Model ───────────────────────────────────────────────────────────────

struct LPModel {
    int n = 0;
    glp_prob* lp = nullptr;
    std::vector<std::vector<int>> x_col;  // x_col[i][j] -> col idx (0 if no arc)
    std::vector<int> y_col;               // y_col[i] -> col idx
    std::vector<int> t_col;               // t_col[i] -> arrival time at node i
    std::vector<std::vector<int>> w_col;  // w_col[i][j] -> McCormick var for x[i][j]*t[i]
    int n_rows_base = 0;

    LPModel() = default;
    ~LPModel() { if (lp) glp_delete_prob(lp); }
    LPModel(const LPModel&) = delete;
    LPModel& operator=(const LPModel&) = delete;

    void build(const Input& inp) {
        n = static_cast<int>(inp.pts.size());
        lp = glp_create_prob();
        glp_set_obj_dir(lp, GLP_MAX);
        glp_term_out(GLP_OFF);

        x_col.assign(n, std::vector<int>(n, 0));
        y_col.resize(n);
        t_col.resize(n, 0);
        w_col.assign(n, std::vector<int>(n, 0));

        // x[i][j] binary arc variables
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (i == j || !std::isfinite(inp.cm[i][j])) continue;
                int c = glp_add_cols(lp, 1);
                glp_set_col_bnds(lp, c, GLP_DB, 0.0, 1.0);
                glp_set_obj_coef(lp, c, 0.0);
                x_col[i][j] = c;
            }

        // y[i] visit indicator
        for (int i = 0; i < n; ++i) {
            int c = glp_add_cols(lp, 1);
            glp_set_col_bnds(lp, c, GLP_DB, 0.0, 1.0);
            glp_set_obj_coef(lp, c, inp.pts[i]);
            y_col[i] = c;
        }
        glp_set_col_bnds(lp, y_col[0], GLP_FX, 1.0, 1.0);

        // t[i] arrival time at node i, continuous in [0, bud_raw]
        for (int i = 0; i < n; ++i) {
            int c = glp_add_cols(lp, 1);
            glp_set_col_bnds(lp, c, GLP_DB, 0.0, inp.bud_raw);
            glp_set_obj_coef(lp, c, 0.0);
            t_col[i] = c;
        }
        // depot departs at time 0
        glp_set_col_bnds(lp, t_col[0], GLP_FX, 0.0, 0.0);

        // w[i][j] = McCormick linearisation of x[i][j] * t[i]
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (!x_col[i][j]) continue;
                int c = glp_add_cols(lp, 1);
                glp_set_col_bnds(lp, c, GLP_DB, 0.0, inp.bud_raw);
                glp_set_obj_coef(lp, c, 0.0);
                w_col[i][j] = c;
            }

        add_flow_constraints();
        add_time_propagation(inp.cm, inp.bud_raw);
        add_mccormick(inp.bud_raw);
        add_fatigue_budget(inp.cm, inp.bud_raw, inp.fatigue_rate);

        n_rows_base = glp_get_num_rows(lp);
    }

    void clone_from(const LPModel& other, const std::vector<std::pair<int, double>>& fixings) {
        assert(lp == nullptr);
        n = other.n;
        x_col = other.x_col;
        y_col = other.y_col;
        t_col = other.t_col;   // added
        w_col = other.w_col;   // added
        n_rows_base = other.n_rows_base;
        lp = glp_create_prob();
        glp_copy_prob(lp, other.lp, GLP_ON);
        for (const auto& [col, val] : fixings)
            glp_set_col_bnds(lp, col, GLP_FX, val, val);
    }

    void add_flow_constraints() {
        // In-flow
        for (int i = 0; i < n; ++i) {
            int r = glp_add_rows(lp, 1);
            glp_set_row_bnds(lp, r, GLP_FX, 0.0, 0.0);
            std::vector<int> ia = {0};
            std::vector<double> ra = {0.0};
            for (int j = 0; j < n; ++j) {
                if (x_col[j][i]) {
                    ia.push_back(x_col[j][i]);
                    ra.push_back(1.0);
                }
            }
            ia.push_back(y_col[i]);
            ra.push_back(-1.0);
            glp_set_mat_row(lp, r, static_cast<int>(ia.size() - 1), ia.data(), ra.data());
        }
        // Out-flow
        for (int i = 0; i < n; ++i) {
            int r = glp_add_rows(lp, 1);
            glp_set_row_bnds(lp, r, GLP_FX, 0.0, 0.0);
            std::vector<int> ia = {0};
            std::vector<double> ra = {0.0};
            for (int j = 0; j < n; ++j) {
                if (x_col[i][j]) {
                    ia.push_back(x_col[i][j]);
                    ra.push_back(1.0);
                }
            }
            ia.push_back(y_col[i]);
            ra.push_back(-1.0);
            glp_set_mat_row(lp, r, static_cast<int>(ia.size() - 1), ia.data(), ra.data());
        }
    }

    void add_budget_constraint(double bud_eff, const std::vector<std::vector<double>>& cm) {
        int r = glp_add_rows(lp, 1);
        glp_set_row_bnds(lp, r, GLP_UP, 0.0, bud_eff);
        std::vector<int> ia = {0};
        std::vector<double> ra = {0.0};
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (x_col[i][j]) { ia.push_back(x_col[i][j]); ra.push_back(cm[i][j]); }
        glp_set_mat_row(lp, r, static_cast<int>(ia.size()-1), ia.data(), ra.data());
    }

    // Time propagation: if arc (i,j) is used, t[j] >= t[i] + cm[i][j]
    // Linearised as: t[j] >= t[i] + cm[i][j] - bud_raw*(1 - x[i][j])
    void add_time_propagation(const std::vector<std::vector<double>>& cm, double bud_raw) {
        for (int i = 0; i < n; ++i)
            for (int j = 1; j < n; ++j) {  // skip depot as destination for propagation
                if (!x_col[i][j]) continue;
                // t[j] - t[i] - bud_raw*x[i][j] >= cm[i][j] - bud_raw
                int r = glp_add_rows(lp, 1);
                glp_set_row_bnds(lp, r, GLP_LO, cm[i][j] - bud_raw, 0.0);
                std::vector<int> ia = {0, t_col[j], t_col[i], x_col[i][j]};
                std::vector<double> ra = {0.0, 1.0, -1.0, -bud_raw};
                glp_set_mat_row(lp, r, 3, ia.data(), ra.data());
            }
    }

    // McCormick envelopes for w[i][j] = x[i][j] * t[i]
    // w >= 0                          (handled by column bound)
    // w >= t[i] - bud_raw*(1-x[i][j])  =>  w - t[i] + bud_raw*x[i][j] >= -bud_raw + bud_raw = 0... 
    // Correct form:
    //   w[i][j] >= t[i] - bud_raw*(1 - x[i][j])  =>  w - t + bud_raw*x >= 0  ... rhs = -bud_raw*(1-1)=0 when x=1
    //   w[i][j] <= bud_raw * x[i][j]
    //   w[i][j] <= t[i]
    void add_mccormick(double bud_raw) {
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (!w_col[i][j]) continue;
                // w >= t[i] - bud_raw*(1 - x[i][j])  =>  w - t[i] + bud_raw*x[i][j] >= 0
                {
                    int r = glp_add_rows(lp, 1);
                    glp_set_row_bnds(lp, r, GLP_LO, 0.0, 0.0);
                    std::vector<int> ia = {0, w_col[i][j], t_col[i], x_col[i][j]};
                    std::vector<double> ra = {0.0, 1.0, -1.0, bud_raw};
                    glp_set_mat_row(lp, r, 3, ia.data(), ra.data());
                }
                // w <= bud_raw * x[i][j]  =>  w - bud_raw*x <= 0
                {
                    int r = glp_add_rows(lp, 1);
                    glp_set_row_bnds(lp, r, GLP_UP, 0.0, 0.0);
                    std::vector<int> ia = {0, w_col[i][j], x_col[i][j]};
                    std::vector<double> ra = {0.0, 1.0, -bud_raw};
                    glp_set_mat_row(lp, r, 2, ia.data(), ra.data());
                }
                // w <= t[i]  =>  w - t[i] <= 0
                {
                    int r = glp_add_rows(lp, 1);
                    glp_set_row_bnds(lp, r, GLP_UP, 0.0, 0.0);
                    std::vector<int> ia = {0, w_col[i][j], t_col[i]};
                    std::vector<double> ra = {0.0, 1.0, -1.0};
                    glp_set_mat_row(lp, r, 2, ia.data(), ra.data());
                }
            }
    }

    // Exact fatigue budget constraint:
    // sum_{i,j} x[i][j]*cm[i][j] + (fatigue_rate/bud_raw)*sum_{i,j} cm[i][j]*w[i][j] <= bud_raw
    void add_fatigue_budget(const std::vector<std::vector<double>>& cm,
                            double bud_raw, double fatigue_rate) {
        int r = glp_add_rows(lp, 1);
        glp_set_row_bnds(lp, r, GLP_UP, 0.0, bud_raw);
        std::vector<int> ia = {0};
        std::vector<double> ra = {0.0};
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (x_col[i][j]) {
                    ia.push_back(x_col[i][j]);
                    ra.push_back(cm[i][j]);
                }
                if (w_col[i][j]) {
                    ia.push_back(w_col[i][j]);
                    ra.push_back((fatigue_rate / bud_raw) * cm[i][j]);
                }
            }
        glp_set_mat_row(lp, r, static_cast<int>(ia.size()-1), ia.data(), ra.data());
    }

    // All four asymmetric cuts for subset S (not containing depot):
    //   1. Directed SEC:     sum_{i,j in S}         x[i][j]  <= |S|-1
    //   2. Outgoing cut:     sum_{i in S, j not in S} x[i][j] >= y[k]  for each k in S
    //   3. Incoming cut:     sum_{i not in S, j in S} x[i][j] >= y[k]  for each k in S
    //   4. Combined in+out:  outgoing + incoming                        >= 2*y[k] (tighter than 2+3 separately)
    void add_sec(const std::vector<int>& S) {
        std::vector<int> Sset(S.begin(), S.end());
        std::sort(Sset.begin(), Sset.end());
        auto in_S = [&](int v) { return std::binary_search(Sset.begin(), Sset.end(), v); };

        // 1. Directed SEC
        {
            int r = glp_add_rows(lp, 1);
            glp_set_row_bnds(lp, r, GLP_UP, 0.0, static_cast<double>(S.size() - 1));
            std::vector<int> ia = {0}; std::vector<double> ra = {0.0};
            for (int ii : S)
                for (int jj : S)
                    if (ii != jj && x_col[ii][jj]) { ia.push_back(x_col[ii][jj]); ra.push_back(1.0); }
            glp_set_mat_row(lp, r, static_cast<int>(ia.size()-1), ia.data(), ra.data());
        }

        // Precompute outgoing and incoming arc lists (shared across cuts 2/3/4)
        std::vector<int> out_ia = {0}, in_ia = {0};
        std::vector<double> out_ra = {0.0}, in_ra = {0.0};
        for (int ii : S)
            for (int j = 0; j < n; ++j)
                if (!in_S(j) && x_col[ii][j]) { out_ia.push_back(x_col[ii][j]); out_ra.push_back(1.0); }
        for (int i = 0; i < n; ++i)
            if (!in_S(i))
                for (int jj : S)
                    if (x_col[i][jj]) { in_ia.push_back(x_col[i][jj]); in_ra.push_back(1.0); }

        // Add cuts once per subset (not per node k) using the most violated node
        // 2+3. Outgoing+Incoming cuts FOR ALL k in S
        for (int k : S) {
            // Outgoing for k
            {
                int r = glp_add_rows(lp, 1);
                glp_set_row_bnds(lp, r, GLP_LO, 0.0, 0.0);
                auto ia = out_ia; auto ra = out_ra;
                ia.push_back(y_col[k]); ra.push_back(-1.0);
                glp_set_mat_row(lp, r, static_cast<int>(ia.size()-1), ia.data(), ra.data());
            }
            // Incoming for k
            {
                int r = glp_add_rows(lp, 1);
                glp_set_row_bnds(lp, r, GLP_LO, 0.0, 0.0);
                auto ia = in_ia; auto ra = in_ra;
                ia.push_back(y_col[k]); ra.push_back(-1.0);
                glp_set_mat_row(lp, r, static_cast<int>(ia.size()-1), ia.data(), ra.data());
            }
        }
        // 4. Combined cut: sum_out + sum_in >= 2*y[S[0]]
        {
            int r = glp_add_rows(lp, 1);
            glp_set_row_bnds(lp, r, GLP_LO, 0.0, 0.0);
            std::vector<int> ia = {0}; std::vector<double> ra = {0.0};
            for (size_t t = 1; t < out_ia.size(); ++t) { ia.push_back(out_ia[t]); ra.push_back(1.0); }
            for (size_t t = 1; t < in_ia.size();  ++t) { ia.push_back(in_ia[t]);  ra.push_back(1.0); }
            ia.push_back(y_col[S[0]]); ra.push_back(-2.0);
            glp_set_mat_row(lp, r, static_cast<int>(ia.size()-1), ia.data(), ra.data());
        }
    }

    bool solve() {
        glp_smcp parm{};
        glp_init_smcp(&parm);
        parm.msg_lev = GLP_MSG_OFF;
        parm.meth = GLP_DUAL;
        int ret = glp_simplex(lp, &parm);
        return ret == 0 && glp_get_status(lp) == GLP_OPT;
    }

    double obj() const { return glp_get_obj_val(lp); }
    double prim(int col) const { return glp_get_col_prim(lp, col); }
    void delete_extra_rows(int base_rows) {
        int cur_rows = glp_get_num_rows(lp);
        if (cur_rows <= base_rows) return;
        std::vector<int> rows_to_del(cur_rows - base_rows + 1);
        rows_to_del[0] = 0;  // GLPK dummy element at index 0
        std::iota(rows_to_del.begin() + 1, rows_to_del.end(), base_rows + 1);
        glp_del_rows(lp, cur_rows - base_rows, rows_to_del.data());
    }
};

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
            if (i != j && model.x_col[i][j] && model.prim(model.x_col[i][j]) > eps) {
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
            if (i != j && model.x_col[i][j] && model.prim(model.x_col[i][j]) > eps)
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

// ── Route extraction & validation ──────────────────────────────────────────

std::vector<int> extract_route(const LPModel& model, double eps = 0.5) {
    const int n = model.n;
    std::vector<int> succ(n, -1);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j && model.x_col[i][j] && model.prim(model.x_col[i][j]) > eps) {
                succ[i] = j;
                break;  // Assume at most one strong successor
            }
        }
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

    while (true) {
        int best_j = -1;
        double best_ratio = -1.0;
        for (int j = 1; j < n; ++j) {
            if (visited[j]) continue;
            double go = inp.cm[cur][j], back = inp.cm[j][0];
            if (!std::isfinite(go) || !std::isfinite(back)) continue;
            if (cost + go + back > inp.bud_eff) continue;

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
        double fm = 1.0 + inp.fatigue_rate * (elapsed / std::max(inp.bud_raw, 1.0));
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
            // insert unvisited node at best position
            std::vector<int> unv;
            for (int j : all_ctrls) if (!new_visited[j]) unv.push_back(j);
            if (unv.empty()) continue;
            int j = unv[randi(0, static_cast<int>(unv.size()) - 1)];
            int bp = 0; double bc = std::numeric_limits<double>::infinity();
            for (int pos = 0; pos <= static_cast<int>(new_route.size()); ++pos) {
                new_route.insert(new_route.begin() + pos, j);
                double tc = rcost(inp.cm, new_route);
                if (tc < bc) { bc = tc; bp = pos; }
                new_route.erase(new_route.begin() + pos);
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

std::vector<int> solve_sa_iterated(const Input& inp, int n_restarts = 4,
                                    int n_iterations = 80000) {
    std::vector<int> best_route;
    double best_score = 0.0;
    std::mt19937 seed_rng(std::random_device{}());
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
    double time_limit_s = 60.0;

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
        root_node.ub = root.solve() ? root.obj() : -std::numeric_limits<double>::infinity();
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
        std::cerr << "Processed " << nodes << " nodes, best: " << best_pts << " pts\n";
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
            if (subtours.empty()) break;

            for (const auto& S : subtours)
                lp.add_sec(S);
        }

        if (!lp.solve()) return;
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
        } else if (branch_col > 0) {
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
        << ", \"elapsed_s\": " << bnc_elapsed << ", \"base_cost\": " << bnc_base
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
};

int main() {
    const std::vector<std::pair<std::string,std::string>> maps = {
        {"op_input_standard.json",     "op_output_standard.json"},
        {"op_input_clustered.json",     "op_output_clustered.json"},
        {"op_input_ring.json",          "op_output_ring.json"},
        {"op_input_path_biased.json",   "op_output_path_biased.json"},
        {"op_input_elev_biased.json",   "op_input_elev_biased.json"},
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
            r.bnc_optimal = (r.bnc_s < 59.0);
            results.push_back(r);
        } catch (...) {}
    }

    // ── Summary table ──────────────────────────────────────────────────────
    std::cout << "\n";
    std::cout << "+----------------------+---------------------------+------------------------------------+\n";
    std::cout << "| Map                  | SA                        | B&C(SA)                            |\n";
    std::cout << "|                      |  pts  nodes     time      |  pts  nodes     time    optimal?   |\n";
    std::cout << "+----------------------+---------------------------+------------------------------------+\n";
    for (const auto& r : results) {
        std::cout << "| " << std::left  << std::setw(20) << r.name
                  << " | " << std::right << std::setw(5)  << static_cast<int>(r.sa_pts)
                  << "  " << std::setw(5) << r.sa_nodes
                  << "  " << std::setw(7) << std::fixed << std::setprecision(1) << r.sa_s << "s"
                  << "   | " << std::setw(5) << static_cast<int>(r.bnc_pts)
                  << "  " << std::setw(5) << r.bnc_nodes
                  << "  " << std::setw(7) << r.bnc_s << "s"
                  << "  " << (r.bnc_optimal ? "YES (proven)" : "no  (limit) ")
                  << " |\n";
    }
    std::cout << "+----------------------+---------------------------+------------------------------------+\n";

    return 0;
}