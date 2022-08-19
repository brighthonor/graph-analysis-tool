/**
    @file   betweennesscentrality.cpp
    @author Gyeonghyo Min
    @date   7/8/2022
*/
/*
    Betweenness Centrality:
    Ulrik Brandes. (2001).
    A Faster Algorithm for Betweenness Centrality

    Robert Geisberger, Peter Sanders, Dominik Schultes (2008)
    Better Approximation of Betweenness Centrality
*/

#include "betweennesscentrality.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <random>
#include <set>

namespace snu {
std::string BetweennessCentrality::statName() {
    return "BetweennessCentrality";
}

bool BetweennessCentrality::calculateStat(Graph &graph, bool verify) {
    const auto &vertices = graph.id_to_vertex;
    int V = vertices.size();
    int sample_sz = std::min(V, MAX_BETWEENNESS_SAMPLE_SZ);

    std::vector<std::pair<Graph::Vid, Graph::Vertex *>> samples;
    std::sample(vertices.begin(), vertices.end(), std::back_inserter(samples),
                sample_sz, std::mt19937{std::random_device{}()});

    bool reversed = false;
    for (auto p : samples) {
        bool suc = applyPartialValue(graph, p.second, betweennessValue, reversed);
        reversed = !reversed;
        if (!suc) {
            std::cout << "-- skipping betweenness centrality for disconnected graph\n";
            return false;
        }
    }
    // adjust constant
    for (auto [vid, b_val] : betweennessValue) {
        betweennessValue[vid] *= (float)V / sample_sz;
    }

    max_betweenness_centrality_id = 0;
    max_betweenness_centrality = 0.0;
    for (auto p : betweennessValue) {
        if (p.second > max_betweenness_centrality) {
            max_betweenness_centrality_id = p.first;
            max_betweenness_centrality = p.second;
        }
    }

    return true;
}

bool BetweennessCentrality::writeToFileStat(std::string graph_name, bool directed) {
    std::string fName = graph_name + "_Betweennes.txt";
    std::ofstream fout(fName.data());
    for (auto [nodeId, between_val] : betweennessValue) {
        fout << nodeId << ' ' << between_val << '\n';
    }
    return true;
}

void BetweennessCentrality::writeToHTMLStat(FILE *fp, bool directed) {
    fprintf(fp,
            "\
            <h2>\
                Betweenness Centrality Statistics\
            </h2>\
            <h3>\
                <p> max betweenness centrality value = %lf at ID = %lld </p>\
            </h3>",
            max_betweenness_centrality, max_betweenness_centrality_id);
}

// runs augmented dijkstra
void BetweennessCentrality::augmented_dijkstra(const Graph &graph, Graph::Vertex *start,
                                               std::unordered_map<Graph::Vid, std::set<Graph::Vid>> &out_predecessors,
                                               std::unordered_map<Graph::Vid, int64_t> &out_pathcount,
                                               std::unordered_map<Graph::Vid, int64_t> &dist,
                                               bool reversed) {
    // initialize output
    dist.clear();
    out_predecessors.clear();
    out_pathcount.clear();
    out_pathcount[start->id] = 1;

    typedef std::pair<Graph::Vertex *, Graph::Vertex *> vertex_pair;
    typedef std::pair<int64_t, vertex_pair> elem;
    std::priority_queue<elem, std::vector<elem>, std::greater<elem>> pq;
    pq.emplace(0, std::make_pair(nullptr, start));

    while (!pq.empty()) {
        elem el = pq.top();
        pq.pop();

        auto cur_dist = el.first;
        auto prev = el.second.first;
        auto cur = el.second.second;

        if (!dist.count(cur->id)) {
            dist[cur->id] = cur_dist;
            for (const auto &edge : cur->edges) {
                auto next_dist = cur_dist + edge->weight;
                bool nextIsTo = edge->to != cur;
                auto nex = nextIsTo ? edge->to : edge->from;

                // prevent going to parent
                if (dist.count(nex->id) > 0 && dist[nex->id] < cur_dist)
                    continue;
                pq.emplace(next_dist, std::make_pair(cur, nex));
            }
        }
        if (prev && dist[cur->id] == cur_dist) {
            out_predecessors[cur->id].insert(prev->id);
            out_pathcount[cur->id] += out_pathcount[prev->id];
        }
    }
}

bool BetweennessCentrality::applyPartialValue(const Graph &graph, Graph::Vertex *A, std::unordered_map<Graph::Vid, float> &betweennessValue, bool reversed) {
    std::unordered_map<Graph::Vid, std::set<Graph::Vid>> predecessors;
    std::unordered_map<Graph::Vid, int64_t> pathcount;
    std::unordered_map<Graph::Vid, int64_t> dist;

    int V = graph.id_to_vertex.size();
    augmented_dijkstra(graph, A, predecessors, pathcount, dist, reversed);
    if ((int)pathcount.size() < V) {
        return false;  // disconnected graph
    }

    typedef std::pair<int64_t, Graph::Vid> dist_info;
    std::vector<dist_info> dist_decreasing;
    dist_decreasing.reserve(V);
    for (auto p : dist) {
        dist_decreasing.emplace_back(p.second, p.first);
    }
    std::sort(dist_decreasing.begin(), dist_decreasing.end(), std::greater<dist_info>());

    std::unordered_map<Graph::Vid, float> P;
    for (auto p : dist_decreasing) {
        // auto d = p.first;
        auto v = p.second;

        P[v] += 1;
        for (auto w : predecessors[v]) {
            if (dist[v] == 0 || pathcount[v] == 0)
                continue;

            float sigmaW = pathcount[w];
            float sigmaV = pathcount[v];
            float dist_ratio = (float)dist[w] / dist[v];
            P[w] += dist_ratio * sigmaW / sigmaV * P[v];
        }
        P[v] -= 1;
    }

    // apply partial value to betweenness centrality
    for (auto p : P) {
        betweennessValue[p.first] += p.second;
    }
    return true;
}

}  // namespace snu