#ifndef GRAPH_ANALYSIS_TOOL_KCLIST_H
#define GRAPH_ANALYSIS_TOOL_KCLIST_H

#include <vector>
#include <map>
#include <queue>

#include "graph.h"
#include "stat.h"

namespace snu {
class Kclist : public UndirectedStat {
   public:
    virtual std::string statName() override;
    std::vector<std::vector<Graph::Vid>> getAllKCliques(USGraph &graph, unsigned k);

   protected:
    virtual bool calculateUndirectedStat(USGraph &graph, bool verify) override;
    virtual void writeToHTMLStat(FILE *fp, bool directed) override;
    virtual bool writeToFileStat(std::string graph_name, bool directed) override;

   private:
    unsigned cliquek;
    unsigned long long cliqueNumber;

    std::vector<std::vector<Graph::Vid>> cliques;

    typedef struct {
        unsigned key;
        unsigned value;
    } KeyValue;

    typedef struct {
        unsigned n_max;
        unsigned n;
        std::vector<unsigned> pt;
        std::vector<KeyValue> kv;
    } BHeap;

    BHeap *construct(unsigned n_max);
    void swap(BHeap *heap, unsigned i, unsigned j);
    void bubbleUp(BHeap *heap, unsigned i);
    void bubbleDown(BHeap *heap);
    void insert(BHeap *heap, KeyValue kv);
    void update(BHeap *heap, unsigned key);
    KeyValue popmin(BHeap *heap);
    BHeap *makeHeap(unsigned n, std::vector<unsigned> &v);
    void freeHeap(BHeap *heap);

    void createVidIndex(USGraph &graph);
    void ordCore(USGraph &graph);
    void reLabel(USGraph &graph);
    void makeSpecial(USGraph &graph, unsigned char k);
    void kclique(unsigned l, unsigned long long *n, std::vector<Graph::Vid> &clq);

    std::vector<unsigned> idx_to_vid;
    std::unordered_map<Graph::Vid, unsigned> vid_to_idx;

    std::vector<unsigned> ns;
    std::vector<std::vector<unsigned>> d;
    std::vector<unsigned> cd;
    std::vector<unsigned> adj;
    std::vector<unsigned> rank;
    std::vector<unsigned char> lab;
    std::vector<std::vector<unsigned>> sub;
};
}


#endif //GRAPH_ANALYSIS_TOOL_KCLIST_H
