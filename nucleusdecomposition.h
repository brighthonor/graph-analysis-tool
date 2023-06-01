#ifndef NUCLEUSDECOMPOSITION_H
#define NUCLEUSDECOMPOSITION_H

#include <vector>
#include <set>
#include <map>
#include <stack>
#include <queue>
#include <algorithm>
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
    void findRSCliques(USGraph &graph, int r, int s, bool inadv);

    // for all
    void fill_rcs_to_id();

    // for Improsive
    void ndWing(snu::USGraph &graph, int *max, int r, int s);
    void fill_conn_deg(USGraph &graph);
    void fill_edge_umap(USGraph &graph);
    int count_scliques(USGraph &graph, std::vector<Graph::Vid> &clique, int r, int s, Graph::Vid maxNew);
    void get_scliques_rnow(USGraph &graph, std::vector<Graph::Vid> &clique, int r, int s, std::vector<std::vector<Graph::Vid>> &result, Graph::Vid maxNew);
    void find_rc_in_sc(std::vector<Graph::Vid> &sclique, std::vector<Graph::Vid> &chosen, int r, int s, int idx, std::set<int> &result);
    int calculate_combination(int s, int r);

    // for inAdv
    void ndInadv(USGraph &graph, int *max, int r, int s);
    void combination(int sidx, std::vector<Graph::Vid> &chosen, int s, int r, int idx);

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

    struct VectorHasher {
        int operator()(const std::vector<Graph::Vid> &V) const;
    };

    void assignToRoot(int* ch, std::vector<subcore>& skeleton);
    void assignToRepresentative(int* ch, std::vector<subcore>& skeleton);
    void store(int uComp, int vComp, std::vector<int>& unassigned, std::vector<std::pair<int, int> >& relations);
    void merge(int u, int v, std::vector<int>& component, std::vector<subcore>& skeleton, int* nSubcores);
    void createSkeleton(int u, const std::set<int> neighbors, int* nSubcores, std::vector<int>& K, std::vector<subcore>& skeleton,
		std::vector<int>& component, std::vector<int>& unassigned, std::vector<std::pair<int, int> >& relations);
    void updateUnassigned(int t, std::vector<int>& component, int* cid, std::vector<std::pair<int, int> >& relations, std::vector<int>& unassigned);
    void buildHierarchy(int cn, std::vector<std::pair<int, int> >& relations, std::vector<subcore>& skeleton, int* nSubcores, int nEdge, int nVtx);

    inline int commons(std::vector<int>& a, std::vector<int>& b);
    void rearrange(std::vector<subcore>& skeleton);
    void reportSubgraph(int r, int s, int index, std::vector<int>& component, std::vector<subcore>& skeleton, Graph& graph, int nEdge, std::vector<nd_tree_node> &nd_tree,
					std::unordered_map<int, int> &skeleton_to_nd_tree, std::vector<bool> visited);
    inline void removeChild(int i, std::vector<subcore> &backup);
    void bfsHierarchy(std::vector<subcore>& skeleton, std::stack<int>& scs);
    inline void findRepresentative(int* child, std::vector<subcore>& skeleton);
    void presentNuclei(int r, int s, std::vector<subcore>& skeleton, std::vector<int>& component, USGraph& graph, int nEdge, std::vector<nd_tree_node> &nd_tree);

    std::vector<std::vector<Graph::Vid>> rcliques;
    std::vector<std::vector<Graph::Vid>> scliques;

    // for improsive
    std::unordered_map<Graph::Vid, Graph::Vid> degree;
    std::unordered_map<Graph::Vid, std::unordered_set<Graph::Vid>> edge_umap;
    std::unordered_map<Graph::Vid, std::unordered_map<Graph::Vid, bool>> conn;

    std::vector<unsigned> r_idx_to_vid;
    std::unordered_map<Graph::Vid, unsigned> r_vid_to_idx;

    // for inAdv
    std::vector<std::set<int>> scsHasr;
    std::vector<std::set<int>> rcsIns;
    std::unordered_map<std::vector<Graph::Vid>, int, VectorHasher> rcs_to_id;

    std::vector<int> K;
    std::vector<nd_tree_node> nd_tree;
};
}

#endif //NUCLEUSDECOMPOSITION_H
