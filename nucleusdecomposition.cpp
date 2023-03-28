#include <fstream>
#include "nucleusdecomposition.h"

namespace snu {
std::string NucleusDecomposition::statName() {
    return "nucleusdecomposition";
}

bool NucleusDecomposition::calculateUndirectedStat(USGraph &graph, bool verify) {
    bool success = true;

    findRSCliques(graph, 3, 4);
    int index = 0;
    for(auto &w: scsHasr){
        printf("%d. vertices: ", index);
        for(auto ww: rcliques[index]){
            printf("%d ", ww);
        } printf("\n");
        for(auto ww: w) {
            for(auto www: scliques[ww]){
                printf("%d ", www);
            } printf("\n");
        } printf("\n");
        index++;
    }

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

void NucleusDecomposition::findRSCliques(USGraph &graph, int r, int s) { // s should be always larger than r
    if(s < r) {
        int tmp = s;
        s = r;
        r = tmp;
    }

    Kclist *rclist = new Kclist();
    Kclist *sclist = new Kclist();

    rcliques = rclist->getAllKCliquesSet(graph, r);
    scliques = sclist->getAllKCliquesSet(graph, s);

    scsHasr.resize(rcliques.size());
    rcsIns.resize(scliques.size());

    for(int i=0; i<rcliques.size(); i++) {
        for(int j=0; j<scliques.size(); j++) {
            auto rc = rcliques[i];
            auto sc = scliques[j];
            bool isin = true;
            for(auto &t: rc) {
                if(!sc.count(t)) {
                    isin = false;
                    break;
                }
            }
            if(isin){
                scsHasr[i].insert(j);
                rcsIns[j].insert(i);
            }
        }
    }
}

void NucleusDecomposition::nd() {

}
} // nmaespace snu
