#ifndef NUCLEUSDECOMPOSITION_H
#define NUCLEUSDECOMPOSITION_H

#include <vector>
#include <set>
#include "kclist.h"
#include "stat.h"

namespace snu {
class NucleusDecomposition : public UndirectedStat {
   public:
    virtual std::string statName() override;

   protected:
    virtual bool calculateUndirectedStat(USGraph &graph, bool verify) override;
    virtual void writeToHTMLStat(FILE *fp, bool directed) override;
    virtual bool writeToFileStat(std::string graph_name, bool directed) override;

   private:
    void nd();
    void findRSCliques(USGraph &graph, int r, int s);

    // nd tree
    struct nd_tree_node {
        int id_;
        int parent_;
        std::vector<int> children_;
        std::vector<int> vertices_;
        uint32_t k_;
        uint32_t r_;
        uint32_t s_;
        uint32_t num_edges_;
        double density_;
    };

    // for nBucket (MinHeap)
    struct KeyValue {
        unsigned key;
        unsigned value;
    };

    struct BHeap {
        unsigned n_max;
        unsigned n;
        std::vector<unsigned> pt;
        std::vector<KeyValue> kv;
    };

    struct BHeap *construct(unsigned n_max);
    void swap(struct BHeap *heap, unsigned i, unsigned j);
    void bubbleUp(struct BHeap *heap, unsigned i);
    void bubbleDown(struct BHeap *heap);
    void insert(struct BHeap *heap, struct KeyValue kv);
    void update(struct BHeap *heap, unsigned key);
    struct KeyValue popmin(struct BHeap *heap);
    struct BHeap *makeHeap();
    void freeHeap(struct BHeap *heap);

    struct subcore {
        bool visible;
        int rank;
        int K;
        int parent;
        int root;
        std::vector<int> children;
        int size;
        int nEdge;
        double ed;

        subcore(int k) {
            K = k;
            rank = 0;
            parent = -1;
            root = -1;
            visible = true;
            size = 0;
            nEdge = 0;
            ed = -1;
        }
    };

    std::vector<std::set<Graph::Vid>> rcliques;
    std::vector<std::set<Graph::Vid>> scliques;
    std::vector<std::set<int>> scsHasr;
    std::vector<std::set<int>> rcsIns;

    int cid;
    std::vector<subcore> skeleton;
    std::vector<int> component;
    std::vector<std::pair<int, int>> relations;
    std::vector<int> unassigned;
    int nSubcores;

    std::vector<int> K;
    std::vector<nd_tree_node> nd_tree;
};
}

#endif //NUCLEUSDECOMPOSITION_H
