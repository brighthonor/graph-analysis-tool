#ifndef CONNECTSTAT_H
#define CONNECTSTAT_H

#include <algorithm>

#include "graph.h"
#include "stat.h"

namespace snu {

class ConnectStat : public CommonStat {
   public:
    virtual std::string statName() override;

   protected:
    virtual bool calculateStat(Graph &graph, bool verify) override;
    virtual void writeToHTMLStat(FILE *fp, bool directed) override;

   private:
    bool calcDS(DSGraph &graph);
    bool calcUS(USGraph &graph);

    void htmlDS(FILE *fp);
    void htmlUS(FILE *fp);

    long long num_cc;     // number of connected components
    long long num_scc;    // number of strongly connected components
    long long num_wcc;    // number of weakly connected components
    long long size_lcc;   // largest size of connected components
    long long size_lscc;  // largest size of strongly connected components
    long long size_lwcc;  // largest size of weakly connected components
    long long diameter;   // longest shortest paths
};
}  // namespace snu

#endif  // CONNECTSTAT_H
