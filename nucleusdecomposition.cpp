#include <fstream>
#include "nucleusdecomposition.h"

namespace snu {
std::string NucleusDecomposition::statName() {
    return "nucleusdecomposition";
}

bool NucleusDecomposition::calculateUndirectedStat(USGraph &graph, bool verify) {
    bool success = true;

    int max;
    findRSCliques(graph, 3, 4);
    nd(graph, &max);

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

//////////// Building nBucket ////////////////
NucleusDecomposition::Naive_Bucket_element::Naive_Bucket_element()
: next(NULL), prev(NULL)
{}

NucleusDecomposition::Naive_Bucket::Naive_Bucket()
: max_value(0), current_min_value(1), elements(NULL), buckets(NULL), values(NULL), nb_elements(0)
{}

NucleusDecomposition::Naive_Bucket::~Naive_Bucket() {
    if(buckets!=NULL)
        free(buckets);
    if(elements!=NULL)
        free(elements);
    if(values!=NULL)
        free(values);
}

void NucleusDecomposition::Naive_Bucket::Free() {
    free(buckets);
    buckets = NULL;
    free(elements);
    elements = NULL;
    free(values);
    values = NULL;
}

void NucleusDecomposition::Naive_Bucket::Initialize(int max_v, int nb_element) {
    int i;
    max_value = max_v;
    buckets = (Naive_Bucket_element **) malloc(sizeof(Naive_Bucket_element *) * (max_value+1));
    elements = (Naive_Bucket_element *) malloc(sizeof(Naive_Bucket_element) * nb_element);
    values = (int *) malloc(sizeof(int) * nb_element);
    nb_elements = nb_element;
    if (buckets == NULL || elements == NULL || values == NULL) {
        free(values);
        free(buckets);
        free(elements);
    }
    else {
        for (i = 0; i <= max_value; i++)
            buckets[i] = NULL;
        for (i = 0; i < nb_element; i++) {
            elements[i].prev = NULL;
            elements[i].next = NULL;
        }
    }
    current_min_value = max_value + 1;
}

void NucleusDecomposition::Naive_Bucket::Insert(int id, int value) {
    values[id] = value;
    elements[id].prev = NULL;
    elements[id].next = buckets[value];
    if (buckets[value] != NULL)
        buckets[value]->prev = &(elements[id]);
    else if (current_min_value > value)
        current_min_value = value;
    buckets[value] = &(elements[id]);
}

int NucleusDecomposition::Naive_Bucket::PopMin(int *id, int *ret_value) {
    for (; current_min_value <= max_value; current_min_value++) {
        if (buckets[current_min_value] != NULL) {
            *id = (int) (buckets[current_min_value] - elements); // pointer arithmetic. finds the index of element that buckets[current_min_value] points to
            buckets[current_min_value] = buckets[current_min_value]->next; // adjust the pointer to the new head of the list that has same degree elements
            if (buckets[current_min_value] != NULL)
                buckets[current_min_value]->prev = NULL;
            *ret_value = values[*id];
            values[*id] = -1;
            return 0;
        }
    }
    return -1; // if the bucket is empty
}

int NucleusDecomposition::Naive_Bucket::CurrentValue(int id) {
    return values[id];
}

void NucleusDecomposition::Naive_Bucket::DecVal(int id) {
    int old_value = values[id];
    // adjust the prev and next pointers
    if (elements[id].prev == NULL)
        buckets[old_value] = elements[id].next;
    else
        elements[id].prev->next = elements[id].next;
    if (elements[id].next != NULL)
        elements[id].next->prev = elements[id].prev;
    Naive_Bucket::Insert(id, values[id]-1);
    return;
}

///////////////////////////////////////////////
//////////// Nucleus Decomposition ////////////
///////////////////////////////////////////////

void NucleusDecomposition::findRSCliques(USGraph &graph, int r, int s) { // s should be always larger than r

    auto *rclist = new Kclist();
    auto *sclist = new Kclist();

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

    K.resize(rcliques.size(), -1);
    Naive_Bucket nBucket;
    nBucket.Initialize(graph.V*graph.V, rcliques.size());

    int id = 0;
    for(int i=0; i<rcliques.size(); i++) {
        if(scsHasr[i].empty()) {
            K[i] = 0;
            continue;
        }
        nBucket.Insert(i, scsHasr[i].size());
    }

    while(true) {
        int t;
        int val;
        if(nBucket.PopMin(&t, &val))
            break;

        unassigned.clear();
        subcore sc (val);
        skeleton.push_back(sc);

        fc_t = K[t] = val;

        for(auto sc: scsHasr[t]) {
            bool nv = true;
            for(auto rc: rcsIns[sc]) {
                if(K[rc]!=-1){
                    nv = false;
                    break;
                }
            }

            if(nv) {
                for(auto rc: rcsIns[sc]) {
                    if(rc==t) continue;
                    if(nBucket.CurrentValue(rc) > fc_t)
                        nBucket.DecVal(rc);
                }
            }
            else {
                rcsIns[sc].erase(t);
                createSkeleton(t, rcsIns[sc]);
                rcsIns[sc].insert(t);
            }
        }

        updateUnassigned(t);
    }

    nBucket.Free();
    *max = fc_t;

    buildHierarchy(*max, graph.E, graph.V);
}

////////////////////////////////////////////////////////////
//////////////// Build Hierarchy Functions /////////////////
////////////////////////////////////////////////////////////

void NucleusDecomposition::assignToRoot(int *ch) {
    std::vector<int> acc;
    int s = *ch;
    while(skeleton[s].root != -1) {
        acc.push_back(s);
        s = skeleton[s].root;
    }
    for(int i: acc)
        skeleton[i].root = s;
    *ch = s;
}

void NucleusDecomposition::assignToRepresentative(int *ch) {
    int u = *ch;
    std::vector<int> vs;
    while(skeleton[u].parent != -1) {
        int n = skeleton[u].parent;
        if(skeleton[n].K == skeleton[u].K) {
            vs.push_back(u);
            u = n;
        }
        else break;
    }
    *ch = u;
    for(int i: vs) {
        if(i != u)
            skeleton[i].parent = u;
    }
}

void NucleusDecomposition::store(int uComp, int vComp) {
    std::pair<int, int> c (vComp, uComp);
    if(uComp == -1)
        unassigned.push_back(relations.size());
    relations.push_back(c);
}

void NucleusDecomposition::merge(int u, int v) {
    if (component[u] == -1) {
        component[u] = component[v];
        skeleton.erase (skeleton.end() - 1);
    }
    else { // merge component[u] and component[v] nodes
        int child = component[u];
        int parent = component[v];
        assignToRepresentative (&child);
        assignToRepresentative (&parent);
        if (child != parent) {
            if (skeleton[child].rank > skeleton[parent].rank)
                std::swap(child, parent);
            skeleton[child].parent = parent;
            skeleton[child].visible = false;
            if (skeleton[parent].rank == skeleton[child].rank)
                skeleton[parent].rank++;
            nSubcores--;
        }
    }
}

void NucleusDecomposition::createSkeleton(int u, std::set<int> neighbors) {
    int smallest = -1, minK = INT32_MAX;
    for(auto i: neighbors) {
        if(K[i] != -1 && K[i] < minK) {
            smallest = i;
            minK = K[i];
        }
    }
    if(smallest==-1) return;

    if(K[smallest] == K[u])
        merge(u, smallest);
    else
        store(component[u], component[smallest]);
}

void NucleusDecomposition::updateUnassigned(int t) {
    if (component[t] == -1) { // if e didn't get a component, give her a new one
        component[t] = cid;
        cid++;
    }

    // update the unassigned components that are in the relations
    for (int i : unassigned)
        relations[i] = std::make_pair(relations[i].first, component[t]);
}

void NucleusDecomposition::buildHierarchy(int cn, int nEdge, int nVtx) {
    // bin the relations w.r.t. first's K
    std::vector<std::vector<std::pair<int, int>>> binnedRelations (cn + 1);

    for (auto r : relations) {
        int a = r.first;
        int b = r.second;
        assignToRepresentative (&a);
        assignToRepresentative (&b);
        if (a == b)
            continue;
        std::pair<int, int> c (a, b);
        binnedRelations[skeleton[a].K].push_back (c);
    }

    // process binnedRelations in reverse order
    for (int i = binnedRelations.size() - 1; i >= 0; i--) {
        std::vector<std::pair<int, int>> mergeList;
        for (int j = 0; j < binnedRelations[i].size(); j++) { // each binnedRelations[i] has K of skeleton[b].K
            int a = binnedRelations[i][j].first;
            int root = binnedRelations[i][j].second;
            assignToRoot (&root);
            if (a != root) {
                if (skeleton[a].K < skeleton[root].K) {
                    skeleton[root].parent = a;
                    skeleton[root].root = a;
                }
                else { // skeleton[root].K == skeleton[a].K
                    std::pair<int, int> c = (root < a) ? std::make_pair(root, a) : std::make_pair(a, root);
                    mergeList.push_back (c);
                }
            }
        }

        // handle merges
        for (auto sc : mergeList) {
            int child = sc.first;
            int parent = sc.second;
            assignToRepresentative (&child);
            assignToRepresentative (&parent);
            if (child != parent) {
                if (skeleton[child].rank > skeleton[parent].rank)
                    std::swap(child, parent);
                skeleton[child].parent = parent;
                skeleton[child].root = parent;
                skeleton[child].visible = false;
                if (skeleton[parent].rank == skeleton[child].rank)
                    skeleton[parent].rank++;
            }
        }
    }

    nSubcores += skeleton.size();

    // root core
    int nid = skeleton.size();
    subcore sc (0);
    for (size_t i = 0; i < skeleton.size(); i++)
        if (skeleton[i].parent == -1)
            skeleton[i].parent = nid;

    sc.size = nVtx;
    sc.nEdge = nEdge;
    sc.ed = nEdge / double (nVtx * (nVtx - 1) / 2);
    skeleton.push_back (sc);
}

} // nmaespace snu
