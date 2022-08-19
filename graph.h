#ifndef GRAPH_H
#define GRAPH_H

#include <list>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "plot.h"

namespace snu {

class Graph;    // graph having common aspects
class DSGraph;  // directed simple graph
class USGraph;  // undirected simple graph

class Graph {
public:
    typedef long long Vid;       // vertex id
    typedef long long Eid;       // edge id
    typedef std::string Vlabel;  // vertex label
    typedef std::string Elabel;  // edge label
    typedef int Weight;          // integer type weight, other types of weight can be added

    class Vertex;
    class Edge;
    class LabelOfVertices;
    class LabelOfEdges;

    class Vertex {
    public:
        Vid id;
        size_t order;
        std::list<LabelOfVertices *> labels;
        std::list<Edge *> edges;
        std::list<Edge *> inbound_edges;
        long long indegree;
        void *temp;
    };

    class Edge {
    public:
        Eid id;
        std::list<LabelOfEdges *> labels;
        Vertex *from;  // In case of undirected graphs, there is no difference between from and to.
        Vertex *to;
        Weight weight;
    };

    class LabelOfVertices {
    public:
        Vlabel label;
        std::list<Vertex *> vertices;
    };

    class LabelOfEdges {
    public:
        Elabel label;
        std::list<Edge *> edges;
    };

    std::unordered_map<Vid, Vertex *> id_to_vertex;
    std::unordered_map<Eid, Edge *> id_to_edge;

    std::unordered_map<Vlabel, LabelOfVertices *> vlabel_to_class;
    std::unordered_map<Elabel, LabelOfEdges *> elabel_to_class;

    // faster access data structures
    std::vector<Vid> vid_order;
    std::vector<Vertex *> vertices;

    size_t V, E; // number of vertices / edges

    struct pair_hash {
        inline unsigned long long operator()(const std::pair<Vid, Vid> &v) const {
            return ((unsigned long long)v.first << 32) + v.second;
        }
    };

    std::unordered_set<std::pair<Vid, Vid>, pair_hash> is_connected;

    long long negative_edge_num = 0;
    char visit = 0;

    Graph();   // prevent from creating Graph class and allow creating graph in subclass
    ~Graph();  // destructor

    virtual int parseGraph(std::string file_path) = 0;

    // add vertex
    // if vertex is created normally, then return 0
    // if error occurs, then return 1
    int addVertex(Vid id, long long num, Vlabel label[]);  // array version
    int addVertex(Vid id, std::vector<Vlabel> *label);     // vector version

    // directed form of addEdge
    // if edge is created normally, then return 0
    // if error occurs, then return 1
    virtual int addEdge(Eid id, long long num, Elabel label[], Vid from, Vid to, Weight weight) = 0;
    int addDirectedEdge(Eid id, long long num, Elabel label[], Vid from, Vid to, Weight weight);  // array version

    // finalize graph
    // call only once after graph is generated
    void finalize();
};

// directed simple graph
class DSGraph : public Graph {
public:
    virtual int parseGraph(std::string file_path) override;
    virtual int addEdge(Eid id, long long num, Elabel label[], Vid from, Vid to, Weight weight) override;  // array version
    int addEdge(Eid id, std::vector<Elabel> *label, Vid from, Vid to, Weight weight);                      // vector version
};

// undirected simple graph
class USGraph : public Graph {
public:
    virtual int parseGraph(std::string file_path) override;
    virtual int addEdge(Eid id, long long num, Elabel label[], Vid from, Vid to, Weight weight) override;  // array version
    int addEdge(Eid id, std::vector<Elabel> *label, Vid from, Vid to, Weight weight);                      // vector version
};
}  // namespace snu

#endif  // GRAPH_H
