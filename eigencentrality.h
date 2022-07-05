#ifndef EIGENCENTRALITY_H
#define EIGENCENTRALITY_H

#include "graph.h"
#include "stat.h"

#define PAGERANK_DAMPING_FACTOR 0.85
#define CONVERGENCE_TEST 1E-5
#define MAX_ITERATIONS 100
#define KATZ_ITERATIONS 100
#define KATZ_ATTENUATION_FACTOR 0.5

namespace snu {
    void eigenCentrality(Graph &graph, StatResult &result);
}

#endif //EIGENCENTRALITY_H
