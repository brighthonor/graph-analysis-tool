#ifndef COUNTSTAT_H
#define COINTSTAT_H

#include "graph.h"
#include "stat.h"

namespace snu {
class CountStat : public UndirectedStat {
   public:
    virtual std::string statName() override;

   protected:
    virtual bool calculateUndirectedStat(USGraph &graph, bool verify) override;
    virtual void writeToHTMLStat(FILE *fp, bool directed) override;

   private:
    long long TriangleCount(std::vector<std::pair<long long, Graph::Vertex *>> &vec);

    unsigned long long wedge_count;     // number of wedges
    unsigned long long claw_count;      // number of claws
    unsigned long long triangle_count;  // number of triangles
};
}  // namespace snu

#endif  // COUNTSTAT_H
