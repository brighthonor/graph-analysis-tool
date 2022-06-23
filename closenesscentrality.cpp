#include "closenesscentrality.h"
#include <map>
#include <queue>
#include <random>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <iostream>

namespace snu {
    // calculates single-source shortest-path for all nodes connected to start.
    static std::map<Graph::Vid, int64_t> dijkstra(const Graph & graph, Graph::Vertex* start) {
        std::map<Graph::Vid, int64_t> dist;

        typedef std::pair<int64_t, Graph::Vertex*> elem;
        std::priority_queue<elem, std::vector<elem>, std::greater<elem>> pq;
        pq.emplace(0, start);

        while(!pq.empty()) {
            elem cur = pq.top();
            pq.pop();

            auto cur_dist = cur.first;
            auto V = cur.second;
            if (!dist.count(V->id)) {
                dist[V->id] = cur_dist;
                for (const auto& edge : V->edges) {
                    auto next_dist = cur_dist + edge->weight;
                    auto to = edge->to == V ? edge->from : edge->to;
                    pq.emplace(next_dist, to);
                }
            }
        }
        return dist;
    }

    static bool calcClosenessCentrality(const Graph &graph, StatResult &result, bool file_output, std::string graph_name) {
        const auto& vertices = graph.id_to_vertex;
        int V = vertices.size();
        int sample_sz = std::min(V, MAX_CLOSENESS_SAMPLE_SZ);

        std::vector<std::pair<Graph::Vid, Graph::Vertex*>> samples;
        std::sample(vertices.begin(), vertices.end(), std::back_inserter(samples), 
                    sample_sz, std::mt19937{std::random_device{}()});

        std::map<Graph::Vid, int64_t> dist_sum;
        for (auto vid : samples) {
            auto dist = dijkstra(graph, vid.second);
            //if ((int)dist.size() != V) return false; // disconnected graph

            for (auto p : dist) {
                dist_sum[p.first] += p.second;
            }
        }

        std::map<Graph::Vid, double> inv_dist_sum;

        // Consider nodes that do not connect anywhere to have closeness centrality of 1.0
        // This is to match SNAP's definition of closeness centrality.

        for (auto [vid, vert] : vertices)
            inv_dist_sum[vid] = 1.0; 

        for (auto p : dist_sum)
            inv_dist_sum[p.first] = (p.second != 0) ? (100.0 / p.second) : 1.0;
        double total_inv_sum = 0;
        for (auto p : inv_dist_sum)
            total_inv_sum += p.second;
        if (total_inv_sum == 0.0) return false;

        Graph::Vid centralNodeId = 0;
        double max_closeness_centrality = 0.0;
        for (auto [nodeId, closeness_val] : inv_dist_sum) {
            if (closeness_val > max_closeness_centrality) {
                centralNodeId = nodeId;
                max_closeness_centrality = closeness_val;
            }
        }


        // output all values in a file
        if (file_output) {
            std::string fName = graph_name + "_Closeness.txt";
            std::ofstream fout(fName.data());
            for (auto [nodeId, closeness_val] : inv_dist_sum) {
                fout << nodeId << ' ' << closeness_val << '\n';
            }
        }
        
        
        // apply results
        result.closenessstat = true;
        result.max_closeness_centrality = max_closeness_centrality / total_inv_sum;
        result.max_closeness_centrality_id = centralNodeId;
        return true;
    }

    bool closenessCentrality(const Graph &graph, StatResult &result, bool file_output, std::string graph_name) {
        //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        bool ret = calcClosenessCentrality(graph, result, file_output, graph_name);
        //std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

        //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
        //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << std::endl;
        return ret;
    }
}