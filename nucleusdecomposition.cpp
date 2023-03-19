#include <fstream>
#include "nucleusdecomposition.h"

namespace snu {
std::string NucleusDecomposition::statName() {
    return "nucleusdecomposition";
}

bool NucleusDecomposition::calculateUndirectedStat(USGraph &graph, bool verify) {
    bool success = true;
    return success;
}

void NucleusDecomposition::writeToHTMLStat(FILE *fp, bool directed) {
    fprintf(fp,
            "\
            <h2>\
                Nucleus Decomposition Statistics\
            </h2>\
            <h3>\
                <p> r = %lld </p>\
                <p> s = %lld </p>\
            </h3>\
            ",
            3, 4);
}

bool NucleusDecomposition::writeToFileStat(std::string graph_name, bool directed) {
    std::ofstream fout(graph_name+"_Nucleusdecompostion.txt");
    fout << "test\n";
    return true;
}

void NucleusDecomposition::test() {
}
} // nmaespace snu
