#ifndef BASICSTAT_H
#define BASICSTAT_H

#include <algorithm>

#include "graph.h"
#include "stat.h"

namespace snu {
class BasicStat : public CommonStat {
   public:
    virtual std::string statName() override;

   protected:
    virtual bool calculateStat(Graph &graph, bool verify) override;
    virtual void writeToHTMLStat(FILE *fp, bool directed) override;

   private:
    bool calcDS(const DSGraph &graph);
    bool calcUS(const USGraph &graph);

    void htmlDS(FILE *fp);
    void htmlUS(FILE *fp);

    long long size;           // number of vertices
    long long volume;         // number of edges
    long long num_vlabel;     // number of vertex label
    long long num_elabel;     // number of edge label
    double density;           // m / n^2, directed graph with loops
    double avg_degree;        // average degree, 2m / n
    long long max_degree;     // maximum degree
    long long max_indegree;   // maximum indegree
    long long max_outdegree;  // maximum outdegree
    // size of vertex label is defined as number of vertices having that label
    long long max_vlabel_size;  // maximum size of vertex label
    long long max_elabel_size;  // maximum size of edge label
    double reciprocity;         // only defined in directed graph
    double negativity;          // ratio of nagative edges
};
}  // namespace snu

#endif  // BASICSTAT_H
