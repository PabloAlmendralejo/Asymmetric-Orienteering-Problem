/*
 * Heuristic comparison: Greedy, GA, ACO, SA on asymmetric OP with fatigue.
 * Reads op_input_*.json files from a directory and compares all four methods.
 * Reports mean, std, and Wilcoxon signed-rank test (SA vs each).
 *
 * Usage: compare_heuristics.exe <input_dir>
 *        compare_heuristics.exe instances
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <random>
#include <iomanip>
#include <filesystem>
#include <functional>
#include <limits>

// ── JSON parser (minimal) ──
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
        return std::stod(json_str.substr(start, pos-start));
    };
    auto parse_array1d = [&](size_t& pos) -> std::vector<double> {
        std::vector<double> v; pos = skip_ws(pos);
        if (json_str[pos]!='[') throw std::runtime_error("Expected ["); ++pos;
        while (true) { pos=skip_ws(pos); if(json_str[pos]==']'){++pos;break;}
            if(json_str[pos]==','){++pos;continue;} v.push_back(parse_number(pos)); }
        return v;
    };
    auto parse_array2d = [&](size_t& pos) -> std::vector<std::vector<double>> {
        std::vector<std::vector<double>> m; pos=skip_ws(pos);
        if(json_str[pos]!='[') throw std::runtime_error("Expected ["); ++pos;
        while(true){pos=skip_ws(pos);if(json_str[pos]==']'){++pos;break;}
            if(json_str[pos]==','){++pos;continue;} if(json_str[pos]=='[') m.push_back(parse_array1d(pos));}
        return m;
    };
    auto find_key = [&](const std::string& key) {
        std::string k="\""+key+"\":"; auto p=json_str.find(k,i);
        if(p==std::string::npos) throw std::runtime_error("Missing: "+key); i=p+k.size();
    };
    find_key("cm"); inp.cm=parse_array2d(i);
    find_key("pts"); inp.pts=parse_array1d(i);
    find_key("bud_eff"); i=skip_ws(i); inp.bud_eff=parse_number(i);
    find_key("bud_raw"); i=skip_ws(i); inp.bud_raw=parse_number(i);
    find_key("fatigue_rate"); i=skip_ws(i); inp.fatigue_rate=parse_number(i);
    return inp;
}

// ── Cost helpers ──
static double rcost(const std::vector<std::vector<double>>& cm, const std::vector<int>& route) {
    if (route.empty()) return 0.0;
    double c = cm[0][route[0]];
    for (size_t i = 0; i+1 < route.size(); ++i) c += cm[route[i]][route[i+1]];
    return c + cm[route.back()][0];
}

static double rcost_fatigue(const std::vector<std::vector<double>>& cm,
                            const std::vector<int>& route, double bud_raw, double fr) {
    if (route.empty()) return 0.0;
    double total=0, elapsed=0;
    std::vector<int> seq={0}; seq.insert(seq.end(),route.begin(),route.end()); seq.push_back(0);
    for (size_t i=0;i+1<seq.size();++i) {
        double leg=cm[seq[i]][seq[i+1]];
        total+=leg*(1.0+fr*(elapsed/std::max(bud_raw,1.0))); elapsed+=leg;
    }
    return total;
}

static double rpts(const std::vector<double>& pts, const std::vector<int>& route) {
    double s=0; for(int v:route) s+=pts[v]; return s;
}

// ── 1. Greedy Insertion ──
std::vector<int> solve_greedy(const Input& inp) {
    int n = (int)inp.pts.size();
    std::vector<bool> vis(n, false); vis[0]=true;
    std::vector<int> route;
    int cur=0;
    while (true) {
        int bj=-1; double br=-1;
        for (int j=1;j<n;++j) {
            if(vis[j]) continue;
            double go=inp.cm[cur][j], back=inp.cm[j][0];
            if(!std::isfinite(go)||!std::isfinite(back)) continue;
            // Quick feasibility: try adding j
            auto trial = route; trial.push_back(j);
            if(rcost_fatigue(inp.cm, trial, inp.bud_raw, inp.fatigue_rate) > inp.bud_raw) continue;
            double ratio=inp.pts[j]/std::max(go,1e-9);
            if(ratio>br){br=ratio;bj=j;}
        }
        if(bj<0) break;
        vis[bj]=true; route.push_back(bj); cur=bj;
    }
    return route;
}

// ── 2. Genetic Algorithm ──
std::vector<int> solve_ga(const Input& inp, int pop_size=50, int generations=200, unsigned seed=42) {
    int n=(int)inp.pts.size();
    std::mt19937 rng(seed);
    auto randu=[&](){return std::uniform_real_distribution<double>(0,1)(rng);};
    auto randi=[&](int lo,int hi){return std::uniform_int_distribution<int>(lo,hi)(rng);};

    // Initialize population with random subsets
    std::vector<std::vector<int>> pop;
    for (int p=0;p<pop_size;++p) {
        std::vector<int> route;
        std::vector<int> perm;
        for(int i=1;i<n;++i) perm.push_back(i);
        std::shuffle(perm.begin(),perm.end(),rng);
        for(int j:perm) {
            auto trial=route; trial.push_back(j);
            if(rcost_fatigue(inp.cm,trial,inp.bud_raw,inp.fatigue_rate)<=inp.bud_raw)
                route.push_back(j);
        }
        pop.push_back(route);
    }

    auto fitness=[&](const std::vector<int>& r)->double{return rpts(inp.pts,r);};

    for(int gen=0;gen<generations;++gen) {
        // Sort by fitness
        std::sort(pop.begin(),pop.end(),[&](auto&a,auto&b){return fitness(a)>fitness(b);});

        // Elitism: keep top 10%
        int elite=std::max(2,pop_size/10);
        std::vector<std::vector<int>> next_pop(pop.begin(),pop.begin()+elite);

        while((int)next_pop.size()<pop_size) {
            // Tournament selection
            int i1=randi(0,pop_size/2), i2=randi(0,pop_size/2);
            auto& p1=pop[std::min(i1,i2)];
            int i3=randi(0,pop_size/2), i4=randi(0,pop_size/2);
            auto& p2=pop[std::min(i3,i4)];

            // Crossover: take nodes from p1, fill gaps from p2
            std::set<int> used;
            std::vector<int> child;
            int cut=p1.empty()?0:randi(0,(int)p1.size()-1);
            for(int k=0;k<=cut&&k<(int)p1.size();++k){child.push_back(p1[k]);used.insert(p1[k]);}
            for(int v:p2) if(!used.count(v)) child.push_back(v);

            // Repair: remove nodes until feasible
            while(!child.empty()&&rcost_fatigue(inp.cm,child,inp.bud_raw,inp.fatigue_rate)>inp.bud_raw) {
                int worst=0; double wv=inp.pts[child[0]];
                for(int k=1;k<(int)child.size();++k)
                    if(inp.pts[child[k]]<wv){wv=inp.pts[child[k]];worst=k;}
                child.erase(child.begin()+worst);
            }

            // Mutation: swap two random positions
            if(child.size()>=2&&randu()<0.3) {
                int a=randi(0,(int)child.size()-1), b=randi(0,(int)child.size()-1);
                std::swap(child[a],child[b]);
                if(rcost_fatigue(inp.cm,child,inp.bud_raw,inp.fatigue_rate)>inp.bud_raw)
                    std::swap(child[a],child[b]); // revert if infeasible
            }
            next_pop.push_back(child);
        }
        pop=next_pop;
    }

    std::sort(pop.begin(),pop.end(),[&](auto&a,auto&b){return fitness(a)>fitness(b);});
    return pop[0];
}

// ── 3. Ant Colony Optimization ──
std::vector<int> solve_aco(const Input& inp, int n_ants=30, int iterations=100, unsigned seed=42) {
    int n=(int)inp.pts.size();
    std::mt19937 rng(seed);
    auto randu=[&](){return std::uniform_real_distribution<double>(0,1)(rng);};

    // Pheromone matrix
    std::vector<std::vector<double>> tau(n, std::vector<double>(n, 1.0));
    double alpha=1.0, beta=2.0, rho=0.1;
    std::vector<int> best_route;
    double best_score=0;

    for(int iter=0;iter<iterations;++iter) {
        std::vector<std::vector<int>> ant_routes(n_ants);
        std::vector<double> ant_scores(n_ants,0);

        for(int a=0;a<n_ants;++a) {
            std::vector<bool> vis(n,false); vis[0]=true;
            std::vector<int> route;
            int cur=0;
            while(true) {
                // Compute probabilities
                std::vector<std::pair<int,double>> candidates;
                double total=0;
                for(int j=1;j<n;++j) {
                    if(vis[j]||!std::isfinite(inp.cm[cur][j])) continue;
                    auto trial=route; trial.push_back(j);
                    if(rcost_fatigue(inp.cm,trial,inp.bud_raw,inp.fatigue_rate)>inp.bud_raw) continue;
                    double attract=std::pow(tau[cur][j],alpha)*std::pow(inp.pts[j]/std::max(inp.cm[cur][j],1e-9),beta);
                    candidates.push_back({j,attract});
                    total+=attract;
                }
                if(candidates.empty()) break;
                // Roulette wheel
                double r=randu()*total, cum=0;
                int chosen=candidates.back().first;
                for(auto&[j,p]:candidates){cum+=p;if(cum>=r){chosen=j;break;}}
                vis[chosen]=true; route.push_back(chosen); cur=chosen;
            }
            ant_routes[a]=route;
            ant_scores[a]=rpts(inp.pts,route);
            if(ant_scores[a]>best_score){best_score=ant_scores[a];best_route=route;}
        }

        // Evaporate
        for(auto&row:tau) for(auto&v:row) v*=(1.0-rho);
        // Deposit
        for(int a=0;a<n_ants;++a) {
            if(ant_scores[a]<=0) continue;
            double deposit=ant_scores[a]/100.0;
            auto&r=ant_routes[a];
            if(!r.empty()) {
                tau[0][r[0]]+=deposit;
                for(size_t k=0;k+1<r.size();++k) tau[r[k]][r[k+1]]+=deposit;
                tau[r.back()][0]+=deposit;
            }
        }
    }
    return best_route;
}

// ── 4. Simulated Annealing ──
std::vector<int> solve_sa(const Input& inp, int n_iters=80000, double T0=100, double Tend=0.1, unsigned seed=42) {
    int n=(int)inp.pts.size();
    std::mt19937 rng(seed);
    auto randu=[&](){return std::uniform_real_distribution<double>(0,1)(rng);};
    auto randi=[&](int lo,int hi){return std::uniform_int_distribution<int>(lo,hi)(rng);};

    // Start from greedy
    auto route=solve_greedy(inp);
    std::vector<bool> vis(n,false); vis[0]=true;
    for(int v:route) vis[v]=true;

    auto best_route=route;
    double best_score=rpts(inp.pts,route), cur_score=best_score;
    double temp=T0, decay=std::pow(Tend/T0,1.0/std::max(n_iters,1));

    for(int it=0;it<n_iters;++it) {
        temp*=decay;
        auto nr=route; auto nv=vis;
        double mv=randu();

        if(mv<0.3&&!nr.empty()) {
            int idx=randi(0,(int)nr.size()-1); nv[nr[idx]]=false; nr.erase(nr.begin()+idx);
        } else if(mv<0.6) {
            std::vector<int> unv; for(int j=1;j<n;++j) if(!nv[j]) unv.push_back(j);
            if(unv.empty()) continue;
            int j=unv[randi(0,(int)unv.size()-1)];
            int bp=0; double bc=1e18;
            for(int pos=0;pos<=(int)nr.size();++pos) {
                int prev=pos==0?0:nr[pos-1], next=pos==(int)nr.size()?0:nr[pos];
                double d=inp.cm[prev][j]+inp.cm[j][next]-inp.cm[prev][next];
                if(d<bc){bc=d;bp=pos;}
            }
            nr.insert(nr.begin()+bp,j); nv[j]=true;
        } else if(mv<0.8&&nr.size()>=2) {
            int a=randi(0,(int)nr.size()-1),b=randi(0,(int)nr.size()-1);
            std::swap(nr[a],nr[b]);
        } else {
            std::vector<int> unv; for(int j=1;j<n;++j) if(!nv[j]) unv.push_back(j);
            if(unv.empty()||nr.empty()) continue;
            int oi=randi(0,(int)nr.size()-1), old=nr[oi], nw=unv[randi(0,(int)unv.size()-1)];
            nr[oi]=nw; nv[old]=false; nv[nw]=true;
        }

        if(rcost_fatigue(inp.cm,nr,inp.bud_raw,inp.fatigue_rate)>inp.bud_raw) continue;
        double ns=rpts(inp.pts,nr), delta=ns-cur_score;
        if(delta>0||(delta<0&&temp>1e-6&&randu()<std::exp(delta/temp))) {
            route=nr; vis=nv; cur_score=ns;
            if(ns>best_score){best_score=ns;best_route=route;}
        }
    }
    return best_route;
}

// ── Statistics ──
struct Stats { double mean, stddev, median, min, max; };

Stats compute_stats(const std::vector<double>& v) {
    Stats s;
    s.mean=std::accumulate(v.begin(),v.end(),0.0)/v.size();
    double sq=0; for(double x:v) sq+=(x-s.mean)*(x-s.mean);
    s.stddev=std::sqrt(sq/v.size());
    auto sorted=v; std::sort(sorted.begin(),sorted.end());
    s.median=sorted[sorted.size()/2];
    s.min=sorted.front(); s.max=sorted.back();
    return s;
}

// Wilcoxon signed-rank test (approximate z-score for n>20)
double wilcoxon_z(const std::vector<double>& a, const std::vector<double>& b) {
    int n=(int)a.size();
    std::vector<std::pair<double,int>> diffs;
    for(int i=0;i<n;++i) {
        double d=a[i]-b[i];
        if(std::abs(d)>1e-9) diffs.push_back({std::abs(d), d>0?1:-1});
    }
    std::sort(diffs.begin(),diffs.end());
    double W_plus=0, W_minus=0;
    for(int i=0;i<(int)diffs.size();++i) {
        double rank=i+1;
        if(diffs[i].second>0) W_plus+=rank; else W_minus+=rank;
    }
    double nn=diffs.size();
    double mu=nn*(nn+1)/4.0;
    double sigma=std::sqrt(nn*(nn+1)*(2*nn+1)/24.0);
    return sigma>0?(W_plus-mu)/sigma:0;
}

// ── Main ──
int main(int argc, char* argv[]) {
    std::string input_dir = "instances";
    if (argc > 1) input_dir = argv[1];

    // Collect all input files
    std::vector<std::string> files;
    for (const auto& entry : std::filesystem::directory_iterator(input_dir)) {
        std::string fname = entry.path().filename().string();
        if (fname.substr(0,9)=="op_input_" && fname.size()>14 &&
            fname.substr(fname.size()-5)==".json")
            files.push_back(entry.path().string());
    }
    std::sort(files.begin(), files.end());
    std::cerr << "Found " << files.size() << " instances\n";

    // Run all methods on each instance (10 random seeds for stochastic methods)
    const int N_SEEDS = 10;
    std::vector<double> all_greedy, all_ga, all_aco, all_sa;
    std::vector<std::string> names;

    // Per-instance best scores for the table
    struct InstanceResult {
        std::string name;
        double greedy, ga_mean, ga_std, aco_mean, aco_std, sa_mean, sa_std;
    };
    std::vector<InstanceResult> results;

    for (const auto& fpath : files) {
        std::string fname = std::filesystem::path(fpath).filename().string();
        std::string name = fname.substr(9, fname.size()-14);
        std::cerr << "  " << name << "...\n";

        Input inp;
        try { inp = parse_input(read_file(fpath)); }
        catch (...) { std::cerr << "    SKIP (parse error)\n"; continue; }

        // Greedy (deterministic)
        double greedy_pts = rpts(inp.pts, solve_greedy(inp));

        // GA, ACO, SA: run N_SEEDS times
        std::vector<double> ga_scores, aco_scores, sa_scores;
        for (int s=0; s<N_SEEDS; ++s) {
            unsigned seed = 1000 + s;
            ga_scores.push_back(rpts(inp.pts, solve_ga(inp, 50, 200, seed)));
            aco_scores.push_back(rpts(inp.pts, solve_aco(inp, 30, 100, seed)));
            sa_scores.push_back(rpts(inp.pts, solve_sa(inp, 80000, 100, 0.1, seed)));
        }

        auto ga_s = compute_stats(ga_scores);
        auto aco_s = compute_stats(aco_scores);
        auto sa_s = compute_stats(sa_scores);

        all_greedy.push_back(greedy_pts);
        all_ga.push_back(ga_s.mean);
        all_aco.push_back(aco_s.mean);
        all_sa.push_back(sa_s.mean);

        results.push_back({name, greedy_pts, ga_s.mean, ga_s.stddev,
                           aco_s.mean, aco_s.stddev, sa_s.mean, sa_s.stddev});
    }

    // ── Per-instance table ──
    std::cout << "\n";
    std::cout << "+----------------------------------------+----------+------------------+------------------+------------------+\n";
    std::cout << "| Instance                               |  Greedy  |  GA (mean±std)   | ACO (mean±std)   |  SA (mean±std)   |\n";
    std::cout << "+----------------------------------------+----------+------------------+------------------+------------------+\n";
    for (const auto& r : results) {
        std::cout << "| " << std::left << std::setw(38) << r.name
                  << " | " << std::right << std::setw(7) << std::fixed << std::setprecision(0) << r.greedy
                  << "  | " << std::setw(7) << r.ga_mean << "±" << std::setw(5) << std::setprecision(1) << r.ga_std
                  << "  | " << std::setw(7) << std::setprecision(0) << r.aco_mean << "±" << std::setw(5) << std::setprecision(1) << r.aco_std
                  << "  | " << std::setw(7) << std::setprecision(0) << r.sa_mean << "±" << std::setw(5) << std::setprecision(1) << r.sa_std
                  << "  |\n";
    }
    std::cout << "+----------------------------------------+----------+------------------+------------------+------------------+\n";

    // ── Aggregate statistics ──
    int n_inst = (int)results.size();
    if (n_inst == 0) { std::cout << "No instances.\n"; return 0; }

    // Compute SA gap vs each method (% below SA)
    std::vector<double> greedy_gap, ga_gap, aco_gap;
    for (int i=0; i<n_inst; ++i) {
        double sa = all_sa[i];
        if (sa > 0) {
            greedy_gap.push_back(100.0 * (sa - all_greedy[i]) / sa);
            ga_gap.push_back(100.0 * (sa - all_ga[i]) / sa);
            aco_gap.push_back(100.0 * (sa - all_aco[i]) / sa);
        }
    }

    auto gg = compute_stats(greedy_gap);
    auto gag = compute_stats(ga_gap);
    auto acg = compute_stats(aco_gap);

    std::cout << "\n── Gap vs SA (% below SA, higher = worse) ──\n";
    std::cout << std::fixed << std::setprecision(1);
    std::cout << "  Greedy: mean=" << gg.mean << "%, std=" << gg.stddev
              << "%, range=[" << gg.min << "%, " << gg.max << "%]\n";
    std::cout << "  GA:     mean=" << gag.mean << "%, std=" << gag.stddev
              << "%, range=[" << gag.min << "%, " << gag.max << "%]\n";
    std::cout << "  ACO:    mean=" << acg.mean << "%, std=" << acg.stddev
              << "%, range=[" << acg.min << "%, " << acg.max << "%]\n";

    // Wilcoxon signed-rank: SA vs each
    double z_greedy = wilcoxon_z(all_sa, all_greedy);
    double z_ga = wilcoxon_z(all_sa, all_ga);
    double z_aco = wilcoxon_z(all_sa, all_aco);

    std::cout << "\n── Wilcoxon signed-rank test (SA vs method, z-score) ──\n";
    std::cout << "  SA vs Greedy: z=" << std::setprecision(2) << z_greedy
              << (std::abs(z_greedy)>1.96?" (p<0.05, significant)":" (not significant)") << "\n";
    std::cout << "  SA vs GA:     z=" << z_ga
              << (std::abs(z_ga)>1.96?" (p<0.05, significant)":" (not significant)") << "\n";
    std::cout << "  SA vs ACO:    z=" << z_aco
              << (std::abs(z_aco)>1.96?" (p<0.05, significant)":" (not significant)") << "\n";

    return 0;
}
