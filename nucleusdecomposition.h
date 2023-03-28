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
};
}

#endif //NUCLEUSDECOMPOSITION_H
