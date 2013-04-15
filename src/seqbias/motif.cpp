/*
 * This file is part of Isolator.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#include "motif.hpp"
#include "logger.hpp"
#include <cmath>
#include <cstdio>

using namespace std;

/* return 4^x */
static unsigned int four_pow(unsigned int x)
{
    return 1 << (2 * x);
}


/* Compute log(exp(x) + exp(y)), avoiding overflow/underflow. */
static double logaddexp(double x, double y)
{
    double u = x - y;
    if (u > 0.0) {
        return x + log1p(exp(-u));
    }else if (u <= 0.0) {
        return y + log1p(exp(u));
    }else  {
        return x + y;
    }
}


/* Bayesian (Schwarz) Information Criterion */
static double bic(double L, double n_obs, double n_params, double c = 1.0)
{
    return L - c * n_params * log(n_obs);

    // I'm fudging this a bit. Technically, it should be:
    //       2.0 * L - c * n_params * log(n_obs)
}


/* A bunch of primative vector/matrix operations */

static void colcpy(double* dest, const double* src, size_t j, size_t n, size_t m)
{
    for(size_t i = 0; i < n; ++i) dest[i] = src[i * m + j];
}

static void vecadd(double* u, double* v, size_t n)
{
    for (size_t i = 0; i < n; ++i) u[i] += v[i];
}

static void vecsub(double* u, double* v, size_t n)
{
    for (size_t i = 0; i < n; ++i) u[i] -= v[i];
}

static void vecaddcol(double* u, double* V, size_t n, size_t m, size_t j)
{
    for (size_t i = 0; i < n; ++i) u[i] += V[i * m + j];
}

static void vecsubcol(double* u, double* V, size_t n, size_t m, size_t j)
{
    for (size_t i = 0; i < n; ++i) u[i] -= V[i * m + j];
}

static void matsetcol(double* U, double* v, size_t n, size_t m, size_t j)
{
    for (size_t i = 0; i < n; ++i) U[i * m + j] = v[i];
}




/* I've moved all of this code out of the motif class, because I think it makes
 * things a bit cleaner. */
class motif_trainer
{
    public:


    motif_trainer(const deque<twobitseq*>& seqs0,
                  const deque<twobitseq*>& seqs1,
                  size_t m,
                  size_t max_parents,
                  size_t max_distance,
                  double complexity_penalty,
                  const char* task_name)
        : col(0)
        , col_max(30)
        , m(m)
        , max_parents(max_parents)
        , max_distance(max_distance)
        , complexity_penalty(complexity_penalty)
        , task_name(task_name)
    {
        M.m = m;
        M.P0 = new kmer_matrix(m, max_parents + 1);
        M.P1 = new kmer_matrix(m, max_parents + 1);

        M.P0->set_all(0.0);
        M.P1->set_all(0.0);

        M.parents = new bool [m * m];
        memset(M.parents, 0, m * m * sizeof(bool));

        M.nparents = new size_t [m];

        sprintf(col_base, "\n%30s", "");

        n0 = seqs0.size();
        n1 = seqs1.size();
        n  = n0 + n1;

        seqs.insert(seqs.begin(), seqs1.begin(), seqs1.end());
        seqs.insert(seqs.begin(), seqs0.begin(), seqs0.end());

        reachability = new bool [m * m];

        L0 = new double [n * m];
        memset(L0, 0, n * m * sizeof(double));

        L1 = new double [n * m];
        memset(L1, 0, n * m * sizeof(double));


        ms0 = new double [n];
        memset(ms0, 0, n * sizeof(double));

        ms1 = new double [n];
        memset(ms1, 0, n * sizeof(double));


        ls0_i = new double [n];
        ls0_j = new double [n];
        ls1_i = new double [n];
        ls1_j = new double [n];


        ps0_i = new double [M.P0->ncols()];
        ps0_j = new double [M.P0->ncols()];
        ps1_i = new double [M.P1->ncols()];
        ps1_j = new double [M.P1->ncols()];

        prior = (double) n1 / (double) (n0 + n1);
    };


    ~motif_trainer()
    {
        delete [] reachability;

        delete [] L0;
        delete [] L1;

        delete [] ms0;
        delete [] ms1;

        delete [] ls0_i;
        delete [] ls0_j;
        delete [] ls1_i;
        delete [] ls1_j;

        delete [] ps0_i;
        delete [] ps0_j;
        delete [] ps1_i;
        delete [] ps1_j;
    }


    void train()
    {
        char msg[512];

        /* How good the current solution is. */
        double ll = conditional_likelihood();
        double ic = bic(ll, n, M.num_params(), complexity_penalty);


        /* Keep track of the best move found so far. */
        move_type move_best = MOVE_NA;
        int i_best = -1;
        int j_best = -1;
        double ic_best;

        int i_search, j_search;
        double ic_search;

        size_t round_num = 0;


        while (true) {
            // XXX: debugging
            if (round_num > 0) {
                break;
            }

            ++round_num;
            Logger::get_task(task_name).inc();
            compute_reachability();

            Logger::debug(col_base);
            snprintf(msg, sizeof(msg), "round %4zu (ic = %0.4e)", round_num, ic);
            Logger::debug(msg);
            Logger::debug(col_base);
            col = 0;

            ic_best = -HUGE_VAL;

            search_additions(i_search, j_search, ic_search);

            if (ic_search > ic_best) {
                ic_best = ic_search;
                i_best  = i_search;
                j_best  = j_search;
                move_best = MOVE_ADDITION;
            }

            search_subtractions(i_search, j_search, ic_search);

            if (ic_search > ic_best) {
                ic_best = ic_search;
                i_best  = i_search;
                j_best  = j_search;
                move_best = MOVE_SUBTRACTION;
            }

            search_reversals(i_search, j_search, ic_search);

            if (ic_search > ic_best) {
                ic_best = ic_search;
                i_best  = i_search;
                j_best  = j_search;
                move_best = MOVE_REVERSAL;
            }

            Logger::debug("\n");

            if (ic_best <= ic) break;

            /* make the best move */

            Logger::debug(col_base);
            switch (move_best) {
                case MOVE_ADDITION:
                    Logger::debug(" [+] ");
                    add_edge(i_best, j_best);
                    break;

                case MOVE_SUBTRACTION:
                    Logger::debug(" [-] ");
                    del_edge(i_best, j_best);
                    break;

                case MOVE_REVERSAL:
                    Logger::debug(" [.] ");
                    del_edge(i_best, j_best);
                    add_edge(j_best, i_best);

                    vecsubcol(ms0, L0, n, m, i_best);
                    vecsubcol(ms1, L1, n, m, i_best);

                    update_likelihood_column(i_best);

                    vecaddcol(ms0, L0, n, m, i_best);
                    vecaddcol(ms1, L1, n, m, i_best);
                    break;

                default:
                    Logger::debug("Inexplicable error. Please report this.");
            }

            vecsubcol(ms0, L0, n, m, j_best);
            vecsubcol(ms1, L1, n, m, j_best);

            update_likelihood_column(j_best);

            vecaddcol(ms0, L0, n, m, j_best);
            vecaddcol(ms1, L1, n, m, j_best);

            snprintf(msg, sizeof(msg), "%d->%d\n", i_best, j_best);
            Logger::debug(msg);

            ic = ic_best;
        }
    }





    private:

    enum move_type {
        MOVE_NA,
        MOVE_ADDITION,
        MOVE_SUBTRACTION,
        MOVE_REVERSAL
    };


    double evaluate_move(int i, int j, move_type move)
    {
        double ll;
        double ic;

        /* keep track of old parameters to avoid retraining */
        M.P0->get_row(j, ps0_j);
        M.P1->get_row(j, ps1_j);
        if (move == MOVE_REVERSAL) {
            M.P0->get_row(i, ps0_i);
            M.P1->get_row(i, ps1_i);
        }

        /* keep track of old likelihoods to avoid reevaluating */
        colcpy(ls0_j, L0, j, n, m);
        colcpy(ls1_j, L1, j, n, m);
        if (move == MOVE_REVERSAL) {
            colcpy(ls0_i, L0, i, n, m);
            colcpy(ls1_i, L1, i, n, m);
        }


        /* make a move */
        switch (move) {
            case MOVE_ADDITION:
                Logger::debug("+");
                add_edge(i, j);
                break;

            case MOVE_SUBTRACTION:
                Logger::debug("-");
                del_edge(i, j);
                break;

            case MOVE_REVERSAL:
                Logger::debug(".");
                rev_edge(i, j);
                break;

            default:
                Logger::abort("Inexplicable error. Please report this.");
        }

        if (++col > col_max) {
            col = 0;
            Logger::debug(col_base);
        }

        /* evaluate likelihood */
        update_likelihood_column(j);
        if (move == MOVE_REVERSAL) {
            update_likelihood_column(i);
        }

        /* update training example likelihoods */
        vecsub(ms0, ls0_j, n);
        vecaddcol(ms0, L0, n, m, j);
        vecsub(ms1, ls1_j, n);
        vecaddcol(ms1, L1, n, m, j);

        if (move == MOVE_REVERSAL) {
            vecsub(ms0, ls0_i, n);
            vecaddcol(ms0, L0, n, m, i);
            vecsub(ms1, ls1_i, n);
            vecaddcol(ms1, L1, n, m, i);
        }

        /* compute likelihood and ic */
        ll = conditional_likelihood();
        ic = bic(ll, n, M.num_params(), complexity_penalty);


        /* undo the move */
        switch (move) {
            case MOVE_ADDITION:
                M.set_edge(i, j, false);
                break;

            case MOVE_SUBTRACTION:
                M.set_edge(i, j, true);
                break;

            case MOVE_REVERSAL:
                M.set_edge(j, i, false);
                M.set_edge(i, j, true);
                break;

            default:
                Logger::abort("Inexplicable error. Please report this.");
        }


        /* restore previous parameters */
        M.P0->set_row(j, ps0_j);
        M.P1->set_row(j, ps1_j);
        if (move == MOVE_REVERSAL) {
            M.P0->set_row(i, ps0_i);
            M.P1->set_row(i, ps1_i);
        }

        /* restore previous likelihoods */
        vecsubcol(ms0, L0, n, m, j);
        vecadd(ms0, ls0_j, n);
        vecsubcol(ms1, L1, n, m, j);
        vecadd(ms1, ls1_j, n);
        if (move == MOVE_REVERSAL) {
            vecsubcol(ms0, L0, n, m, i);
            vecadd(ms0, ls0_i, n);
            vecsubcol(ms1, L1, n, m, i);
            vecadd(ms1, ls1_i, n);
        }


        matsetcol(L0, ls0_j, n, m, j);
        matsetcol(L1, ls1_j, n, m, j);
        if (move == MOVE_REVERSAL) {
            matsetcol(L0, ls0_i, n, m, i);
            matsetcol(L1, ls1_i, n, m, i);
        }

        return ic;
    }




    void search_additions(
            int& i_best,
            int& j_best,
            double& ic_best)
    {
        int i, j;

        double ic;

        i_best = 0;
        j_best = 0;
        ic_best = -HUGE_VAL;

        /* specifies the range of i's that we search through */
        int i_start, i_end;

        for (j = m - 1; j >= 0; --j) {

            /* determine the range of i's that form possible (i, j) edges */
            if (M.has_edge(j, j)) {
                if (max_distance == 0 || (size_t) j < max_distance) i_start = 0;
                else i_start = j - max_distance;

                if (max_distance == 0) i_end = m - 1;
                else i_end = min(m - 1, j + max_distance);
            }
            else i_start = i_end = j;

            for (i = i_start; i <= i_end; ++i) {
                /* skip existing edges */
                if (M.has_edge(i, j)) continue;

                /* skip edges that would introduce cycles */
                if (reachable(j, i)) continue;

                /* skip edges that would exceed the parent limit */
                if (M.num_parents(j) >= max_parents) continue;

                /* skip edges that are equivalent to one already tried:
                 * (If i and j both have in-degree 0, and i > j, then we already
                 * tried the edge (j, i).)
                 */
                if (i > j && M.num_parents(j) == 1 && M.num_parents(i) == 1) continue;


                ic = evaluate_move(i, j, MOVE_ADDITION);

                if (ic >= ic_best) {
                    i_best = i;
                    j_best = j;
                    ic_best = ic;
                }
            }
        }
    }



    void search_subtractions(
            int& i_best,
            int& j_best,
            double& ic_best)
    {
        int i, j;
        double ic;

        i_best = 0;
        j_best = 0;
        ic_best = -HUGE_VAL;

        for (j = 0; (size_t) j < m; ++j) {
            for (i = 0; (size_t) i < m; ++i) {

                /* skip non-existent edges */
                if (!M.has_edge(i, j)) continue;

                /* skip nodes with in-edges */
                if (i == j && M.num_parents(j) > 1) continue;

                ic = evaluate_move(i, j, MOVE_SUBTRACTION);

                if (ic >= ic_best) {
                    i_best = i;
                    j_best = j;
                    ic_best = ic;
                }
            }
        }
    }


    void search_reversals(
            int& i_best,
            int& j_best,
            double& ic_best)
    {
        int i, j;
        double ic;
        bool has_ij_path;

        i_best = 0;
        j_best = 0;
        ic_best = -HUGE_VAL;

        for (j = 0; (size_t) j < m; ++j) {
            for (i = 0; (size_t) i < m; ++i) {

                /* skip nodes */
                if (i == j) continue;

                /* skip non-existent edges */
                if (!M.has_edge(i, j)) continue;

                /* skip reversals that add parameters */
                if (!M.has_edge(i, i) || !M.has_edge(j, j)) continue;

                /* skip reversals that would introduce cycles
                 * (this is determined by the cutting the (i, j) edge,
                 * recomputing reachability and checking if i, j remains) */
                M.set_edge(i, j, false);
                compute_reachability();

                has_ij_path = reachable(i, j);

                M.set_edge(i, j, true);
                compute_reachability();

                if (has_ij_path) continue;

                ic = evaluate_move(i, j, MOVE_REVERSAL);

                if (ic >= ic_best) {
                    i_best = i;
                    j_best = j;
                    ic_best = ic;
                }
            }
        }
    }


    void compute_reachability()
    {
        /* find the transitive closure of the parents matrix via flyod-warshall */

        memcpy(reachability, M.parents, m * m * sizeof(bool));

        size_t k, i, j;
        for (k = 0; k < m; ++k) {
            for (i = 0; i < m; ++i) {
                for (j = 0; j < m; ++j) {
                    reachability[j * m + i] =
                        reachability[j * m + i] ||
                        (reachability[k * m + i] && reachability[j * m + k]);
                }
            }
        }
    }

    double conditional_likelihood()
    {
        double ll = 0.0;
        size_t i;

        double p = log(prior);
        double q = log(1.0 - prior);

        /* background sequence likelihood */
        for (i = 0; i < n0; ++i) {
            ll += (q + ms0[i]) - logaddexp(q + ms0[i], p + ms1[i]);
        }

        /* foreground sequnece likelihood */
        for (i = n0; i < n; ++i) {
            ll += (p + ms1[i]) - logaddexp(q + ms0[i], p + ms1[i]);
        }

        return ll;
    }


    void update_likelihood_column(int j)
    {
        size_t i;
        kmer K;

        // zero columns
        for (i = 0; i < n; ++i) {
            L0[i * m + j] = 0.0;
            L1[i * m + j] = 0.0;
        }

        // compute likelihood
        deque<twobitseq*>::const_iterator seq;
        for (i = 0, seq = seqs.begin(); seq != seqs.end(); ++i, ++seq) {
            if ((*seq)->make_kmer(K, 0, M.parents + j * m, m) > 0) {
                L0[i * m + j] = (*M.P0)(j, K);
                L1[i * m + j] = (*M.P1)(j, K);
            }
        }
    }



    void calc_row_params(int j)
    {
        /* initialize */

        M.P0->set_row(j, 0.0);
        M.P1->set_row(j, 0.0);

        size_t np = M.num_parents(j);

        /* If this node is now removed, do nothing further. */
        if (np == 0) return;

        size_t K_max = four_pow(np);
        kmer K;

        for (K = 0; K < K_max; ++K) {
            (*M.P0)(j, K) = pseudocount;
            (*M.P1)(j, K) = pseudocount;
        }


        /* tabulate */

        size_t c;
        deque<twobitseq*>::const_iterator seq;
        for (c = 0, seq = seqs.begin(); seq != seqs.end(); ++c, ++seq) {

            if ((*seq)->make_kmer(K, 0, M.parents + j * m, m) == 0) continue;

            /* background sequence */
            if (c < n0) (*M.P0)(j, K) += 1;

            /* foreground sequence */
            else (*M.P1)(j, K) += 1;
        }


        /* normalize */

        size_t np_pred = 0; // predecessor parents
        int u;
        for (u = 0; u < j; ++u) {
            if (M.has_edge(u, j)) ++np_pred;
        }

        M.P0->make_conditional_distribution(j, np_pred, np);
        M.P1->make_conditional_distribution(j, np_pred, np);

        M.P0->dist_log_transform_row(j, np);
        M.P1->dist_log_transform_row(j, np);
    }


    void add_edge(int i, int j)
    {
        M.set_edge(i, j, true);
        calc_row_params(j);
    }


    void del_edge(int i, int j)
    {
        M.set_edge(i, j, false);
        calc_row_params(j);
    }


    void rev_edge(int i, int j)
    {
        del_edge(i, j);
        add_edge(j, i);
    }


    bool reachable(int i, int j)
    {
        return reachability[j * m + i];
    }

    public:
    motif M;

    private:
    deque<twobitseq*> seqs;


    /* for printing pretty output */
    int col;
    int col_max;
    char col_base[80];

    /** a n*n 0-1 matrix marking if node i is reachable from node j */
    bool* reachability;

    /* number of training examples (foreground, background, and total, resp.) */
    size_t n0, n1, n;

    /* number of positions */
    size_t m;

    double prior;
    
    size_t max_parents;
    size_t max_distance;
    double complexity_penalty;


    double* L0;
    double* L1;

    double* ms0;
    double* ms1;

    double* ls0_i;
    double* ls0_j;
    double* ls1_i;
    double* ls1_j;

    double* ps0_i;
    double* ps0_j;
    double* ps1_i;
    double* ps1_j;

    const char* task_name;

    static const double pseudocount;
};


const double motif_trainer::pseudocount = 1.0;


motif::motif()
    : m(0)
    , P0(NULL)
    , P1(NULL)
    , parents(NULL)
    , nparents(NULL)
{

}


motif::motif(const motif& other)
{
    m = other.m;
    P0 = new kmer_matrix(*other.P0);
    P1 = new kmer_matrix(*other.P0);
    parents = new bool [m * m];
    memcpy(parents, other.parents, m * m * sizeof(bool));

    nparents = new size_t [m];
    memcpy(nparents, other.nparents, m * sizeof(size_t));
}


motif::motif(const deque<twobitseq*>& seqs0,
             const deque<twobitseq*>& seqs1,
             size_t m,
             size_t max_parents,
             size_t max_distance,
             double complexity_penalty,
             const char* task_name)
    : nparents(NULL)
{
    /* We are essentially outsourcing construction here. The motif_trainer
     * object does all the construction, because it is complicated, and I want
     * to encapsulate it a bit. */

    motif_trainer trainer(seqs0, seqs1,
                          m,
                          max_parents,
                          max_distance,
                          complexity_penalty,
                          task_name);

    trainer.train();

    this->m = trainer.M.m;
    P0 = new kmer_matrix(*trainer.M.P0);
    P1 = new kmer_matrix(*trainer.M.P1);
    parents = new bool [m * m];
    memcpy(parents, trainer.M.parents, m * m * sizeof(bool));

    nparents = new size_t [m];
    for (size_t i = 0; i < m; ++i) nparents[i] = num_parents(i);
}




#if 0
motif::motif(const YAML::Node& node)
{
    /* m */
    unsigned int m_;
    node["m"] >> m_;
    m = (size_t) m_;


    /* parents */
    parents = new bool [m * m];
    memset(parents, 0, m * m * sizeof(bool));
    const YAML::Node& parents_node = node["parents"];

    size_t i;
    int b;
    for (i = 0; i < m * m; ++i) {
        parents_node[i] >> b;
        parents[i] = (bool) b;
    }

    nparents = new size_t [m];
    for (size_t i = 0; i < m; ++i) nparents[i] = num_parents(i);


    /* P0 */
    const YAML::Node& P0_node = node["P0"];
    P0 = new kmer_matrix(P0_node);


    /* P1 */
    const YAML::Node& P1_node = node["P1"];
    P1 = new kmer_matrix(P1_node);
}
#endif


#if 0
void motif::to_yaml(YAML::Emitter& out) const
{
    out << YAML::BeginMap;

    /* m */
    out << YAML::Key   << "m";
    out << YAML::Value << (unsigned int) m;


    /* parents */
    out << YAML::Key << "parents";
    out << YAML::Value;
    out << YAML::Flow << YAML::BeginSeq;

    size_t i;
    for (i = 0; i < m * m; ++i) out << (parents[i] ? 1 : 0);

    out << YAML::EndSeq;


    /* background parameters */
    out << YAML::Key << "P0";
    out << YAML::Value;
    P0->to_yaml(out);


    /* foreground parameters */
    out << YAML::Key << "P1";
    out << YAML::Value;
    P1->to_yaml(out);

    out << YAML::EndMap;
}
#endif


string motif::model_graph(int offset) const
{
    string graph_str;
    char strbuf[512];

    graph_str += "digraph {\n";
    graph_str += "splines=\"true\";\n";
    graph_str += "node [shape=\"box\"];\n";

    /* print nodes */
    size_t i, j;
    for (j = 0; j < m; j++) {
        snprintf(strbuf, sizeof(strbuf), 
                 "n%d [label=\"%d\",pos=\"%d,0\",style=\"%s\"];\n",
                 (int) j, (int) j - offset, (int) j * 100,
                 parents[j * m + j] ? "solid" : "dotted");

        graph_str += strbuf;
    }

    /* print edges */
    for (j = 0; j < m; ++j) {
        if (!parents[j * m + j]) continue;

        for (i = 0; i < m; ++i) {
            if (i == j) continue;
            if (parents[j * m + i]) {
                snprintf(strbuf, sizeof(strbuf),
                         "n%zu -> n%zu;\n", i, j);
                graph_str += strbuf;
            }
        }
    }

    graph_str += "}\n";
    return graph_str;
}


double motif::eval(const twobitseq& seq, size_t offset) const
{
    double ll0 = 0.0;
    double ll1 = 0.0;
    kmer K;

    size_t n = P0->nrows();
    size_t i;
    for (i = 0; i < n; ++i) {
        if (nparents != NULL && nparents[i] == 0) continue;
        if (seq.make_kmer(K, offset, parents + i * m, m) == 0) continue;

        ll0 += (*P0)(i, K);
        ll1 += (*P1)(i, K);
    }

    return exp(ll1 - ll0);
}



void motif::set_edge(int i, int j, bool x)
{
    parents[j * m + i] = x;
}


bool motif::has_edge(int i, int j) const
{
    return parents[j * m + i];
}

size_t motif::num_parents(size_t j) const
{
    size_t i;
    size_t N = 0;
    for (i = 0; i < m; ++i) if (has_edge(i, j)) ++N;
    
    return N;
}

size_t motif::num_params() const
{
    size_t N =0;
    int i;
    for (i = 0; (size_t) i < m; ++i) {
        N += four_pow(num_parents(i)) - 1;
    }

    N *= 2; /* foreground and background parameters */

    return N;
}

motif::~motif()
{
    delete P0;
    delete P1;

    delete [] parents;
    delete [] nparents;
}



