#include "stat.h"

namespace snu {

bool Stat::getSuccess() {
    return success;
}

void Stat::writeToHTML(FILE *fp, bool directed) {
    if (success)
        writeToHTMLStat(fp, directed);
    else
        fprintf(stderr, "writeToHTML called on failed stat %s\n", statName().c_str());
}

bool Stat::writeToFile(std::string graph_name, bool directed) {
    if (success)
        return writeToFileStat(graph_name, directed);
    else
        fprintf(stderr, "writeToFile called on failed stat %s\n", statName().c_str());
    return false;
}

bool CommonStat::calculate(Graph &graph, bool verify) {
    success = calculateStat(graph, verify);
    return success;
}

bool DirectedStat::calculateDirected(DSGraph &graph, bool verify) {
    success = calculateDirectedStat(graph, verify);
    return success;
}

bool UndirectedStat::calculateUndirected(USGraph &graph, bool verify) {
    success = calculateUndirectedStat(graph, verify);
    return success;
}
}  // namespace snu
