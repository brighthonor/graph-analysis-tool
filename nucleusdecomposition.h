#ifndef NUCLEUSDECOMPOSITION_H
#define NUCLEUSDECOMPOSITION_H

#include <vector>
#include <set>
#include <stack>
#include <queue>
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
    void nd(USGraph &graph, int *max);
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
    struct Naive_Bucket_element{
        Naive_Bucket_element *prev;
        Naive_Bucket_element *next;
        Naive_Bucket_element();
    };

    struct Naive_Bucket {
        Naive_Bucket_element **buckets;
        Naive_Bucket_element *elements;
        int nb_elements;
        int max_value;
        int *values;
        int current_min_value;

        Naive_Bucket();
        ~Naive_Bucket();
        void Initialize(int max_v, int nb_element);
        void Free ();
        void Insert(int id, int value);
        void DecVal(int id);
        int PopMin(int* ret_id, int* ret_value);
        int CurrentValue(int id);
    };

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

    void assignToRoot(int *ch);
    void assignToRepresentative(int *ch);
    void store(int uComp, int vComp);
    void merge(int u, int v);
    void createSkeleton(int u, std::set<int> neighbors);
    void updateUnassigned(int t);
    void buildHierarchy(int cn, int nEdge, int nVtx);

    inline int commons(std::vector<int> &a, std::list<Graph::Edge *> edges, int u);
    void rearrange();
    void reportSubgraph(int r, int s, int index, USGraph &graph, std::unordered_map<int, int> &skeleton_to_nd_tree, std::vector<bool> visited);
    void bfsHierarchy(std::stack<int> &scs);
    inline void findRepresentative(int *child);
    void presentNuclei(int r, int s, USGraph &graph);

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
