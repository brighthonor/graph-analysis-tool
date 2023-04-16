/*
 *
 * Author: Sanghyun Yi
 * Date: 2023/03/19
 * reference: https://github.com/maxdan94/kClist/blob/master/kClist.c
 *
 */

#include <fstream>
#include "kclist.h"

namespace snu {
std::string Kclist::statName() {
    return "kclist";
}

std::vector<std::vector<Graph::Vid>> Kclist::getAllKCliques(USGraph &graph, unsigned k, bool vid) {
    createVidIndex(graph);
    ordCore(graph);
    reLabel(graph);

    cliquek = k;
    makeSpecial(graph, cliquek);

    cliqueNumber = 0;
    std::vector<Graph::Vid> tmp;
    std::set<Graph::Vid> tmpset;
    kclique(cliquek, &cliqueNumber, tmp, tmpset, vid);

    for(auto &v: cliques) {
        std::sort(v.begin(), v.end());
    }

    return cliques;
}

std::vector<std::set<Graph::Vid>> Kclist::getAllKCliquesSet(USGraph &graph, unsigned int k) {
    createVidIndex(graph);
    ordCore(graph);
    reLabel(graph);

    cliquek = k;
    makeSpecial(graph, cliquek);

    cliqueNumber = 0;
    std::vector<Graph::Vid> tmp;
    std::set<Graph::Vid> tmpset;
    kclique(cliquek, &cliqueNumber, tmp, tmpset, true);

    return cliques_set;
}

bool Kclist::calculateUndirectedStat(USGraph &graph, bool verify) {
    bool success = true;

    createVidIndex(graph);
    ordCore(graph);
    reLabel(graph);

    cliquek = 3;
    makeSpecial(graph, cliquek);

    cliqueNumber = 0;
    std::vector<Graph::Vid> tmp;
    std::set<Graph::Vid> tmpset;
    kclique(cliquek, &cliqueNumber, tmp, tmpset, true);

    return success;
}

void Kclist::rel_idx_and_vid(std::vector<unsigned int> &itv, std::unordered_map<Graph::Vid, unsigned int> &vti) {
    for(auto &n: idx_to_vid) {
        itv.push_back(n);
    }
    for(auto &p: vid_to_idx) {
        vti.insert(p);
    }
}

void Kclist::writeToHTMLStat(FILE *fp, bool directed) {
    fprintf(fp,
            "\
            <h2>\
                KCList (listing all k-cliquese) Statistics\
            </h2>\
            <h3>\
                <p> k = %d </p>\
                <p> num of k-cliques = %lld </p>\
            </h3>\
            ",
            cliquek, cliqueNumber);
}

bool Kclist::writeToFileStat(std::string graph_name, bool directed) {
    std::ofstream fout(graph_name+"_Kclist.txt");
    fout << "kclist\n";
    return true;
}

////////////// BHeap Construct //////////////

Kclist::BHeap *Kclist::construct(unsigned n_max) {
    unsigned i;
    auto *heap = new BHeap();
    heap->n_max = n_max;
    heap->n = 0;
    heap->pt.resize(n_max);
    for(i=0; i<n_max; i++) heap->pt[i] = -1;
    heap->kv.resize(n_max);
    return heap;
}

void Kclist::swap(BHeap *heap, unsigned i, unsigned j) {
    KeyValue kv_tmp = heap->kv[i];
    unsigned pt_tmp = heap->pt[kv_tmp.key];
    heap->pt[heap->kv[i].key] = heap->pt[heap->kv[j].key];
    heap->kv[i] = heap->kv[j];
    heap->pt[heap->kv[j].key] = pt_tmp;
    heap->kv[j] = kv_tmp;
}

void Kclist::bubbleUp(BHeap *heap, unsigned i) {
    unsigned j = (i-1)/2;
    while(i>0) {
        if(heap->kv[j].value > heap->kv[i].value) {
            swap(heap, i, j);
            i = j;
            j = (i-1)/2;
        } else break;
    }
}

void Kclist::bubbleDown(BHeap *heap) {
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

void Kclist::insert(BHeap *heap, KeyValue kv) {
    heap->pt[kv.key] = (heap->n)++;
    heap->kv[heap->n-1] = kv;
    bubbleUp(heap, heap->n-1);
}

void Kclist::update(BHeap *heap, unsigned key) {
    int i = (int) heap->pt[key];
    if(i!=-1) {
        ((heap->kv[i]).value)--;
        bubbleUp(heap, i);
    }
}

Kclist::KeyValue Kclist::popmin(BHeap *heap) {
    KeyValue min = heap->kv[0];
    heap->pt[min.key] = -1;
    heap->kv[0] = heap->kv[--(heap->n)];
    heap->pt[heap->kv[0].key] = 0;
    bubbleDown(heap);
    return min;
}

/////////// Building BHeap ///////////////
Kclist::BHeap *Kclist::makeHeap(unsigned n, std::vector<unsigned> &v) {
    unsigned i;
    KeyValue kv;
    BHeap *heap = construct(n);

    for(i=0; i<n; i++) {
        kv.key = i;
        kv.value = v[i];
        insert(heap, kv);
    }
    return heap;
}

void Kclist::freeHeap(BHeap *heap) {
    free(heap);
}

/////////// Vertex Indexing //////////
void Kclist::createVidIndex(USGraph &graph) {
    unsigned n=graph.V;

    Graph::Vid index = 0;

    idx_to_vid.resize(n);
    for(auto &p: graph.id_to_vertex) {
        idx_to_vid[index] = p.first;
        vid_to_idx[p.first] = index;
        index += 1;
    }
}

//////////// CORE ORDERING //////////////
void Kclist::ordCore(USGraph &graph) {
    unsigned i, j, r=0;
    unsigned n=graph.V, e=graph.E;
    KeyValue kv;
    BHeap *heap;

    std::vector<unsigned> d0(n);
    std::vector<unsigned> cd0(n+1);
    std::vector<unsigned> adj0(2*e);

    // degree
    for(auto &p: graph.id_to_edge) {
        unsigned f = vid_to_idx[p.second->from->id];
        unsigned t = vid_to_idx[p.second->to->id];
        d0[f]++;
        d0[t]++;
    }

    // cumulate degree
    cd0[0] = 0;
    for(i=1; i<n+1; i++) {
        cd0[i] = cd0[i-1] + d0[i-1];
        d0[i-1] = 0;
    }

    // truncated list of neighbors
    for(auto &p: graph.id_to_edge) {
        unsigned f = vid_to_idx[p.second->from->id];
        unsigned t = vid_to_idx[p.second->to->id];
        adj0[cd0[f] + d0[f]++] = t;
        adj0[cd0[t] + d0[t]++] = f;
    }

    // create BHeap
    heap = makeHeap(n, d0);

    // fill Rank - Ranking of the nodes according to degeneracy ordering
    rank.resize(n);
    for(i=0; i<n; i++) {
        kv = popmin(heap);
        rank[kv.key] = n - (++r);
        for(j=cd0[kv.key]; j<cd0[kv.key+1]; j++) {
            update(heap, adj0[j]);
        }
    }

    freeHeap(heap);
}

///////// RELABEL EDGES TO MAKE DAG ///////////
void Kclist::reLabel(USGraph &graph) {
    unsigned source, target;

    for(auto &p: graph.id_to_edge) {
        unsigned f = vid_to_idx[p.second->from->id];
        unsigned t = vid_to_idx[p.second->to->id];
        source = rank[f];
        target = rank[t];
        if(source < target) {
            auto tmp = p.second->from;
            p.second->from = p.second->to;
            p.second->to = tmp;
        }
    }
}

void Kclist::makeSpecial(USGraph &graph, unsigned char k) {
    unsigned V = graph.V, E = graph.E;

    unsigned i, lns, max;
    std::vector<unsigned> ld(V), lsub(V);
    std::vector<unsigned char> llab(V);

    for(auto &p: graph.id_to_edge) {
        unsigned f = vid_to_idx[p.second->from->id];
        ld[f]++;
    }

    cd.resize(V+1);
    lns = 0;
    cd[0] = 0;
    max = 0;
    lsub.resize(V);
    llab.resize(V);
    for(i=1; i<V+1; i++) {
        cd[i] = cd[i-1] + ld[i-1];
        max = (max > ld[i-1]) ? max : ld[i-1];
        lsub[lns++] = i-1;
        ld[i-1] = 0;
        llab[i-1] = k;
    }

    adj.resize(E);
    for(auto &p: graph.id_to_edge) {
        unsigned f = vid_to_idx[p.second->from->id];
        unsigned t = vid_to_idx[p.second->to->id];
        adj[cd[f] + ld[f]++] = t;
    }

    ns.resize(k+1);
    ns[k] = lns;
    d.resize(k+1);
    sub.resize(k+1);
    for(i=2; i<k; i++) {
        d[i].resize(V);
        sub[i].resize(max);
    }
    d[k] = ld;
    sub[k] = lsub;
    lab = llab;
}

void Kclist::kclique(unsigned int l, unsigned long long *n, std::vector<Graph::Vid> &clq, std::set<Graph::Vid> &clqset, bool vid) {
    unsigned i, j, k, end, u, v, w;

    if(l==2) {
        for(i=0; i<ns[2]; i++) {
            u = sub[2][i];

            if(vid) clq.push_back(idx_to_vid[u]);
            else clq.push_back(u);
            clqset.insert(idx_to_vid[u]);

            end = cd[u] + d[2][u];
            for(j=cd[u]; j<end; j++) {

                if(vid) clq.push_back(idx_to_vid[adj[j]]);
                else clq.push_back(adj[j]);
                clqset.insert(idx_to_vid[adj[j]]);

                cliques.push_back(clq);
                cliques_set.push_back(clqset);

                clq.pop_back();
                clqset.erase(idx_to_vid[adj[j]]);

                (*n)++; // listing here!
            }
            clq.pop_back();
            clqset.erase(idx_to_vid[u]);
        }
        return;
    }

    for(i=0; i<ns[l]; i++) {
        u = sub[l][i];
        ns[l-1] = 0;
        end = cd[u] + d[l][u];
        for(j=cd[u]; j<end; j++) { // forming U'
            v = adj[j];
            lab[v] = l-1;
            sub[l-1][ns[l-1]++] = v;
            d[l-1][v] = 0;
        }
        for(j=0; j<ns[l-1]; j++) { // reordering adjacency list and computing new degrees
            v = sub[l-1][j];
            end = cd[v] + d[l][v];
            for(k=cd[v]; k<end; k++) {
                w = adj[k];
                if(lab[w] == l-1) {
                    d[l-1][v]++;
                } else {
                    adj[k--] = adj[--end];
                    adj[end] = w;
                }
            }
        }

        if(vid) clq.push_back(idx_to_vid[u]);
        else clq.push_back(u);
        clqset.insert(idx_to_vid[u]);

        kclique(l-1, n, clq, clqset, vid);

        clq.pop_back();
        clqset.erase(idx_to_vid[u]);

        for(j=0; j<ns[l-1]; j++) { // restoring labels
            v = sub[l-1][j];
            lab[v] = l;
        }
    }
}
}