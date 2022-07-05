#ifndef BETWEENNESSCENTRALITY_H
#define BETWEENNESSCENTRALITY_H

#include "graph.h"
#include "stat.h"

#define MAX_BETWEENNESS_SAMPLE_SZ 100

namespace snu {
    bool betweennessCentrality(const Graph &graph, StatResult &result);
}

#endif //BETWEENNESSCENTRALITY_H
