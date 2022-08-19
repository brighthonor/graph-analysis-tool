/**
    @file   closenesscentrality.cpp
    @author Gyeonghyo Min
    @date   7/8/2022
*/
/*
    Closeness Centrality:
    David Eppstein, Joseph Wang. (2000).
    Fast Approximation of Centrality.
*/
#include "closenesscentrality.h"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <random>
#include <sstream>

namespace snu {
std::string ClosenessCentrality::statName() {
    return "ClosenessCentrality";
}

bool ClosenessCentrality::calculate(Graph &graph, bool outbound,
                                    std::map<Graph::Vid, double> &cc, double &max, long long &max_id) {
    const auto &vertices = graph.id_to_vertex;
    int V = vertices.size();
    int sample_sz = std::min(V, MAX_CLOSENESS_SAMPLE_SZ);

    std::vector<std::pair<Graph::Vid, Graph::Vertex *>> samples;
    std::sample(vertices.begin(), vertices.end(), std::back_inserter(samples),
                sample_sz, std::mt19937{std::random_device{}()});

    std::map<Graph::Vid, int64_t> dist_sum;
    for (auto vid : samples) {
        auto dist = dijkstra(graph, vid.second, outbound);
        // if ((int)dist.size() != V) return false; // disconnected graph

        for (auto p : dist) {
            dist_sum[p.first] += p.second;
        }
    }

    // Consider nodes that do not connect anywhere to have closeness centrality of 1.0
    // This is to match SNAP's definition of closeness centrality.

    for (auto [vid, vert] : vertices)
        cc[vid] = 1.0;

    for (auto p : dist_sum)
        cc[p.first] = (p.second != 0) ? ((double)(std::min(sample_sz - 1, (int)dist_sum.size() - 1)) / p.second) : 1.0;
    double total_inv_sum = 0;
    for (auto p : cc)
        total_inv_sum += p.second;
    if (total_inv_sum == 0.0)
        return false;

    max_id = 0;
    max = 0.0;
    for (auto [nodeId, closeness_val] : cc) {
        if (closeness_val > max) {
            max_id = nodeId;
            max = closeness_val;
        }
    }
    return true;
}

bool ClosenessCentrality::calculateStat(Graph &graph, bool verify) {
    return calculate(graph, false, closeness_centrality, max_closeness_centrality, max_closeness_centrality_id) &&
           calculate(graph, true, outbound_closeness_centrality, max_outbound_closeness_centrality,
                     max_outbound_closeness_centrality_id);
}

// calculates single-source shortest-path for all nodes connected to start.
std::map<Graph::Vid, int64_t> ClosenessCentrality::dijkstra(const Graph &graph, Graph::Vertex *start, bool outbound) {
    std::map<Graph::Vid, int64_t> dist;

    typedef std::pair<int64_t, Graph::Vertex *> elem;
    std::priority_queue<elem, std::vector<elem>, std::greater<elem>> pq;
    pq.emplace(0, start);

    while (!pq.empty()) {
        elem cur = pq.top();
        pq.pop();

        auto cur_dist = cur.first;
        auto V = cur.second;
        if (!dist.count(V->id)) {
            dist[V->id] = cur_dist;
            // this may be confusing, but note that the sampling is done in *reverse*
            // (from the viewpoint of the destination, not the source)
            for (const auto &edge : (outbound ? V->inbound_edges : V->edges)) {
                auto next_dist = cur_dist + edge->weight;
                auto to = edge->to == V ? edge->from : edge->to;
                pq.emplace(next_dist, to);
            }
        }
    }
    return dist;
}

void ClosenessCentrality::writeMap(std::string fname, std::map<Graph::Vid, double> map) {
    std::ofstream fout(fname.data());
    for (auto [nodeId, val] : map) {
        fout << nodeId << ' ' << val << '\n';
    }
}

bool ClosenessCentrality::writeToFileStat(std::string graph_name, bool directed) {
    writeMap(graph_name + "_Closeness.txt", closeness_centrality);
    if (directed) {
        writeMap(graph_name + "_OutboundCloseness.txt", outbound_closeness_centrality);
    }
    return true;
}

void ClosenessCentrality::writeToHTMLStat(FILE *fp, bool directed) {
    fprintf(fp,
            "\
            <h2>\
                Closeness Centrality Statistics\
            </h2>\
            <h3>\
                <p> max closeness centrality value = %lf at ID = %lld </p>",
            max_closeness_centrality, max_closeness_centrality_id);
    if (directed) {
        fprintf(fp,
                "<p> max outbound closeness centrality value = %lf at ID = %lld </p>",
                max_outbound_closeness_centrality, max_outbound_closeness_centrality_id);
    }
    fprintf(fp, "</h3>");
}
}  // namespace snu