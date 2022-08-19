#ifndef BICONNECTED_H
#define BICONNECTED_H

#include <map>
#include "graph.h"
#include "stat.h"

namespace snu {
class BiconnectedComponents : public UndirectedStat {
   public:
    virtual std::string statName() override;

   protected:
    virtual bool calculateUndirectedStat(USGraph &graph, bool verify) override;
    virtual void writeToHTMLStat(FILE *fp, bool directed) override;
    virtual bool writeToFileStat(std::string graph_name, bool directed) override;

   private:
    void countBcc(USGraph &graph);
    void resetVisited(USGraph &graph); // reset visited
    bool verifyArticulationPoints(USGraph &graph); // for verification, O(V^2)

    long long num_arp;       // number of articulation points
    long long num_bcc;       // number of biconnected components
    long long max_conn_bcc;  // maximum number of bcc's connected to a single arp
    long long size_lbcc;     // largest size of biconnected components

    std::map<Graph::Vid, long long> connected_bcc;
};
}  // namespace snu

#endif  // BICONNECTED_H
