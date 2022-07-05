#include "biconnected.h"
#include <stack>

namespace snu {
    typedef struct VertexMetadata {
        bool visited = false;
        long long num_connected_bcc = 1; // number of connected bcc's (1 for non-articulation points)
        long long bcc_size = 1; // size of the bcc the vertex is located in
        long long depth = -1; // the depth in the dfs
        long long low_point = -1; // depth of topmost ancestor reachable
        long long child_count = 0; // number of children in dfs tree
        Graph::Vertex *parent = nullptr; // the parent in the dfs tree
    } VertexMetadata;

    static VertexMetadata *getMetadata(Graph::Vertex *v) {
        return (VertexMetadata *) v->temp;
    }

    // This is the algorithm presented by Hopcroft and Tarjan (1973).
    // It has been modified to avoid recursion (to prevent stack overflows)
    static void countBcc(Graph &graph, StatResult &result) {
        std::stack<Graph::Vertex *> dfs_stack;

        // Attempt dfs starting from every vertex:
        // this accounts for disconnected graphs.
        for (auto& pair: graph.id_to_vertex) {
            Graph::Vertex *root = pair.second;
            VertexMetadata *root_meta = getMetadata(root);
            if (!root_meta->visited) { // begin dfs with this root!
                root_meta->depth = root_meta->low_point = 0;
                root_meta->num_connected_bcc = 0;
                dfs_stack.push(root);

                while (!dfs_stack.empty()) {
                    bool escape = false; // used to escape nested loop
                    Graph::Vertex *v = dfs_stack.top();
                    VertexMetadata *meta = getMetadata(v);
                    meta->visited = true;

                    for (auto& e: v->edges) {
                        Graph::Vertex *to = e->to == v ? e->from : e->to;
                        VertexMetadata *to_meta = getMetadata(to);

                        if (!to_meta->visited) { // this is a child!
                            meta->child_count++;

                            // initialize child metadata
                            to_meta->parent = v;
                            to_meta->depth = to_meta->low_point = meta->depth + 1;
                            dfs_stack.push(to);

                            escape = true; // escape loop, begin processing child immediately
                            break;
                        } else if (to != meta->parent) {
                            if (meta->low_point > to_meta->depth) { // this is a back edge!
                                meta->low_point = to_meta->depth;
                            }
                        }
                    }
                    if (escape) {
                        continue;
                    }
                    if (v == root) { // if i'm root
                        break; // finished processing root, stack should be empty
                    }
                    bool parent_is_arp = meta->depth == 1
                                      || meta->low_point == meta->depth - 1
                                      || meta->low_point == meta->depth;
                    if (parent_is_arp) {
                        getMetadata(meta->parent)->num_connected_bcc++;
                        result.num_bcc++;
                        if (meta->bcc_size + 1 > result.size_lbcc) {
                            result.size_lbcc = meta->bcc_size + 1;
                        }
                    } else { // propagate bcc size to parent
                        getMetadata(meta->parent)->bcc_size += meta->bcc_size;
                    }
                    // update low point of parent before popping
                    if (getMetadata(meta->parent)->low_point > meta->low_point) {
                        getMetadata(meta->parent)->low_point = meta->low_point;
                    }
                    dfs_stack.pop();
                }
            }
        }
    }

    void biconnectedComponents(Graph& graph, StatResult& result) {
        for (auto& pair: graph.id_to_vertex) {
            pair.second->temp = new VertexMetadata();
        }

        result.num_bcc = 0;
        result.size_lbcc = 0;
        countBcc(graph, result);

        // count articulation points
        result.num_arp = 0;
        result.max_conn_bcc = 0;
        for (auto& pair: graph.id_to_vertex) {
            long long cbcc = getMetadata(pair.second)->num_connected_bcc;
            if (cbcc > 1) {
                result.num_arp++;
                if (cbcc > result.max_conn_bcc) {
                    result.max_conn_bcc = cbcc;
                }
            }
        }

        result.biconnectedstat = true;

        for (auto& pair: graph.id_to_vertex) {
            delete getMetadata(pair.second);
        }
    }
}

