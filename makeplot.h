#ifndef MAKEPLOT_H
#define MAKEPLOT_H

#include "graph.h"
#include "plot.h"

namespace snu {
void makePlot(Graph &graph, Plot &plot);
void makeDSPlot(DSGraph &graph, Plot &plot);
void makeUSPlot(USGraph &graph, Plot &plot);
}  // namespace snu

#endif  // MAKEPLOT_H
