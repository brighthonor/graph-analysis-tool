#include "eigencentrality.h"

#include <cmath>

namespace snu {
std::string EigenCentrality::statName() {
    return "EigenCentrality";
}

bool EigenCentrality::calculateStat(Graph &graph, bool verify) {
    std::vector<double> probv;
    probv.resize(graph.V);

    if ((max_eigenvalue = calcEigenCentrality(graph, probv)) >= 1E-12) {
        normalizeProb((double) graph.V, probv);
        eigencentrality_converged = true;
        max_eigencentrality = 0;
        for (size_t i = 0; i < graph.V; i++) {
            if (probv[i] > max_eigencentrality) {
                max_eigencentrality = probv[i];
                max_eigencentrality_id = graph.vid_order[i];
            }
        }
    }

    if (calcPageRank(graph, probv)) {
        normalizeProb((double) graph.V, probv);
        pagerank_converged = true;
        max_pagerank = 0;
        for (size_t i = 0; i < graph.V; i++) {
            if (probv[i] > max_pagerank) {
                max_pagerank = probv[i];
                max_pagerank_id = graph.vid_order[i];
            }
        }
    }

    calcKatzCentrality(graph, probv, max_eigenvalue);
    normalizeProb((double) graph.V, probv);
    katz_centrality_computed = true;
    max_katz_centrality = 0;
    for (size_t i = 0; i < graph.V; i++) {
        if (probv[i] > max_katz_centrality) {
            max_katz_centrality = probv[i];
            max_katz_centrality_id = graph.vid_order[i];
        }
    }

    return true;
}

void EigenCentrality::writeToHTMLStat(FILE *fp, bool directed) {
    fprintf(fp,
            "\
            <h2>\
                Eigenvector Centrality Statistics\
            </h2>\
            <h3>\
            <p> max eigenvalue = %lf</p>",
            max_eigenvalue);
    if (eigencentrality_converged) {
        fprintf(fp, "<p> max eigenvector centrality value = %lf at ID = %lld </p>",
                max_eigencentrality, max_eigencentrality_id);
    } else {
        fprintf(fp, "<p> eigenvector centrality failed to converge. </p>");
    }
    if (pagerank_converged) {
        fprintf(fp, "<p> max PageRank value = %lf at ID = %lld </p>",
                max_pagerank, max_pagerank_id);
    } else {
        fprintf(fp, "<p> PageRank failed to converge. </p>");
    }
    if (katz_centrality_computed) {
        fprintf(fp, "<p> max Katz centrality value = %lf at ID = %lld </p>",
                max_katz_centrality, max_katz_centrality_id);
    }
    fprintf(fp, "</h3>");
}

void EigenCentrality::normalizeProb(double n, std::vector<double> &probv) {
    double sum = 0;
    for (double p: probv) {
        sum += p;
    }
    for (double& p: probv) {
        p *= n / sum;
    }
}

bool EigenCentrality::calcPageRank(Graph &graph, std::vector<double> &probv) {
    std::vector<double> mulv;
    mulv.resize(graph.V);

    std::fill(probv.begin(), probv.end(), 1);

    bool not_finished = true;
    int iterations = 0;
    double damper = 1. - PAGERANK_DAMPING_FACTOR;
    while (not_finished && iterations < MAX_ITERATIONS) {
        iterations++;
        not_finished = false;
        std::fill(mulv.begin(), mulv.end(), 0); // prepare matrix multiplication result
        for (size_t i = 0; i < graph.V; i++) {
            auto g = graph.vertices[i];
            if (!g->edges.empty()) {
                double w = 1. / (double) g->edges.size();
                for (auto& edge: g->edges) {
                    auto to = edge->to == g ? edge->from : edge->to;
                    mulv[to->order] += w * probv[g->order];
                }
            }
        }
        for (size_t i = 0; i < graph.V; i++) {
            double next = damper + PAGERANK_DAMPING_FACTOR * mulv[i];
            if (std::abs(next - probv[i]) > CONVERGENCE_TEST) {
                not_finished = true;
            }
            probv[i] = next;
        }
    }

    return !not_finished;
}

double EigenCentrality::calcEigenCentrality(Graph &graph, std::vector<double> &probv) {
    std::vector<double> mulv;
    mulv.resize(graph.V);
    double eigenvalue;

    std::fill(probv.begin(), probv.end(), 1);

    bool not_finished = true;
    int iterations = 0;
    while (not_finished && iterations < MAX_ITERATIONS) {
        iterations++;
        not_finished = false;
        std::fill(mulv.begin(), mulv.end(), 0); // prepare matrix multiplication result
        for (size_t i = 0; i < graph.V; i++) {
            auto g = graph.vertices[i];
            if (!g->edges.empty()) {
                for (auto& edge: g->edges) {
                    auto to = edge->to == g ? edge->from : edge->to;
                    mulv[to->order] += probv[g->order];
                }
            }
        }
        eigenvalue = 0;
        for (double m: mulv) {
            eigenvalue += m * m;
        }
        eigenvalue = std::sqrt(eigenvalue);
        if (eigenvalue <= 1E-12) {
            eigenvalue = 0;
            break;
        }
        for (size_t i = 0; i < graph.V; i++) {
            double next = mulv[i] / eigenvalue;
            if (std::abs(probv[i] - next) > CONVERGENCE_TEST) {
                not_finished = true;
            }
            probv[i] = next;
        }
    }

    return eigenvalue;
}

void EigenCentrality::calcKatzCentrality(Graph &graph, std::vector<double> &probv, double eigenvalue) {
    std::vector<double> powv, tmpv;
    powv.resize(graph.V, 1);
    tmpv.resize(graph.V);

    double att = eigenvalue <= 0 ? 0.8 : std::min(1. / eigenvalue, 1.) * 0.8;

    for (int it = 0; it < KATZ_ITERATIONS; it++) {
        std::fill(tmpv.begin(), tmpv.end(), 0);
        for (size_t i = 0; i < graph.V; i++) {
            auto g = graph.vertices[i];
            if (!g->edges.empty()) {
                for (auto& edge: g->edges) {
                    auto to = edge->to == g ? edge->from : edge->to;
                    tmpv[to->order] += powv[g->order] * att;
                }
            }
        }
        std::swap(powv, tmpv);

        for (size_t i = 0; i < graph.V; i++) {
            probv[i] += powv[i];
        }
    }
}
}  // namespace snu
