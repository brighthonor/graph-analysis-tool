#ifndef NUCLEUSDECOMPOSITION_H
#define NUCLEUSDECOMPOSITION_H

#include <vector>
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
    void test();

    std::vector<std::vector<Graph::Vid>> rcliques;
    std::vector<std::vector<Graph::Vid>> scliques;
};
}

#endif //NUCLEUSDECOMPOSITION_H
