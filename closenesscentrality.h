#ifndef CLOSENESSCENTRALITY_H
#define CLOSENESSCENTRALITY_H

#include "graph.h"
#include "stat.h"

#define MAX_CLOSENESS_SAMPLE_SZ 100

namespace snu {
    bool closenessCentrality(const Graph &graph, StatResult &result);
}

#endif //CLOSENESSCENTRALITY_H
