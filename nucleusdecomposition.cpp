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

////////////// BHeap Construct //////////////
struct NucleusDecomposition::BHeap *NucleusDecomposition::construct(unsigned n_max) {
    unsigned i;
    auto *heap = new BHeap();
    heap->n_max = n_max;
    heap->n = 0;
    heap->pt.resize(n_max);
    for(i=0; i<n_max; i++) heap->pt[i] = -1;
    heap->kv.resize(n_max);
    return heap;
}

void NucleusDecomposition::swap(struct BHeap *heap, unsigned i, unsigned j) {
    KeyValue kv_tmp = heap->kv[i];
    unsigned pt_tmp = heap->pt[kv_tmp.key];
    heap->pt[heap->kv[i].key] = heap->pt[heap->kv[j].key];
    heap->kv[i] = heap->kv[j];
    heap->pt[heap->kv[j].key] = pt_tmp;
    heap->kv[j] = kv_tmp;
}

void NucleusDecomposition::bubbleUp(struct BHeap *heap, unsigned i) {
    unsigned j = (i-1)/2;
    while(i>0) {
        if(heap->kv[j].value > heap->kv[i].value) {
            swap(heap, i, j);
            i = j;
            j = (i-1)/2;
        } else break;
    }
}

void NucleusDecomposition::bubbleDown(struct BHeap *heap) {
    unsigned i=0, j1=1, j2=2, j;
    while(j1 < heap->n) {
        j = ((j2 < heap->n) && (heap->kv[j2].value < heap->kv[j1].value)) ? j2 : j1;
        if(heap->kv[j].value < heap->kv[i].value) {
            swap(heap, i, j);
            i = j;
            j1 = 2*i+1;
            j2 = j1+1;
            continue;
        }
        break;
    }
}

void NucleusDecomposition::insert(struct BHeap *heap, struct KeyValue kv) {
    heap->pt[kv.key] = (heap->n)++;
    heap->kv[heap->n-1] = kv;
    bubbleUp(heap, heap->n-1);
}

void NucleusDecomposition::update(struct BHeap *heap, unsigned key) {
    unsigned i = heap->pt[key];
    if(i!=-1) {
        ((heap->kv[i]).value)--;
        bubbleUp(heap, i);
    }
}

struct NucleusDecomposition::KeyValue NucleusDecomposition::popmin(struct BHeap *heap) {
    KeyValue min = heap->kv[0];
    heap->pt[min.key] = -1;
    heap->kv[0] = heap->kv[--(heap->n)];
    heap->pt[heap->kv[0].key] = 0;
    bubbleDown(heap);
    return min;
}

/////////// Building BHeap ///////////////
struct NucleusDecomposition::BHeap* NucleusDecomposition::makeHeap() {
    unsigned i;
    struct KeyValue kv;
    struct BHeap *heap = construct(scsHasr.size());

    K.resize(scsHasr.size(), -1);

    for(i=0; i<scsHasr.size(); i++) {
        if(scsHasr[i].empty()){
            K[i] = 0;
            continue;
        }
        kv.key = i;
        kv.value = scsHasr[i].size();
        insert(heap, kv);
    }
    return heap;
}

void NucleusDecomposition::freeHeap(struct BHeap *heap) {
    free(heap);
}

//////////// Nucleus Decomposition ////////////
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

void NucleusDecomposition::nd(USGraph &graph, int *max) {
    int fc_t = 0;

    cid = 0;
    nSubcores = 0;
    component.resize(rcliques.size(), -1);

    struct BHeap *heap = makeHeap();

    for(int i=0; i<rcliques.size(); i++) {
        int t;
        int val;
        struct KeyValue kv = popmin(heap);
        t = kv.key;
        val = kv.value;

        unassigned.clear();
        subcore sc (val);
        skeleton.push_back(sc);

        fc_t = K[t] = val;

        for(auto sc: scsHasr[t]) {
            bool nv = false;
            for(auto rc: rcsIns[sc]) {
                if(K[rc]!=-1){
                    nv = true;
                    break;
                }
            }

            if(nv) {
                for(auto rc: rcsIns[sc]) {
                    if(heap->pt)
                }
            }
        }
    }
}

} // nmaespace snu
