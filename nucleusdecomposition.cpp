#include <fstream>
#include <string>
#include "nucleusdecomposition.h"

namespace snu {
std::string NucleusDecomposition::statName() {
    return "nucleusdecomposition";
}

bool NucleusDecomposition::calculateUndirectedStat(USGraph &graph, bool verify) {
    bool success = true;

    int max;
    FILE *fp = fopen("test.txt", "r");
    int r = 3, s = 4;
    bool inadv = true;
    char buffer[3];
    char fileName[200];
    fread(buffer, sizeof(char), 2, fp);
    r = atoi(&buffer[0])/10;
    s = atoi(&buffer[1]);
    printf("[DEBUG] r: %d, s: %d", r, s);
    printf("[DEBUG] START...\n");
    findRSCliques(graph, r, s, inadv);
    printf("[DEBUG] findRSCliques Completed...\n");
    if(inadv) ndInadv(graph, &max, r, s);
    else ndWing(graph, &max, r, s);
    printf("[DEBUG] nd Completed...\n");

    ///////////////////////DEBUG/////////////////////////
    sprintf(fileName, "nd_result_%d%d.txt", r, s);
    fp = fopen(fileName, "w");
    fprintf(fp, "[INFO] nd_tree\n");
    for(auto nd: nd_tree) {
        fprintf(fp, "[#%d nd_tree]\n\tparent: %d\n", nd.id_, nd.parent_);
        fprintf(fp, "\tk: %d, r: %d, s: %d\n", nd.k_, nd.r_, nd.s_);
        // printf("\tvertices: ");
        // for(auto v: nd.vertices_) printf("%d ", v);
        // printf("\n");
        // printf("\tchildren: ");
        // for(auto c: nd.children_) printf("%d ", c);
        // printf("\n");
        fprintf(fp, "\t# of edges: %d, edge density: %f\n", nd.num_edges_, nd.density_);
    } fprintf(fp, "\n");
    ///////////////////////DEBUG/////////////////////////

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
            3LL, 4LL);
}

bool NucleusDecomposition::writeToFileStat(std::string graph_name, bool directed) {
    std::ofstream fout(graph_name+"_Nucleusdecompostion.txt");
    fout << "test\n";
    return true;
}

//////////// Building nBucket ////////////////
NucleusDecomposition::Naive_Bucket_element::Naive_Bucket_element()
: prev(NULL), next(NULL)
{}

NucleusDecomposition::Naive_Bucket::Naive_Bucket()
: buckets(NULL), elements(NULL), nb_elements(0), max_value(0), values(NULL), current_min_value(1)
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

int NucleusDecomposition::VectorHasher::operator()(const std::vector<Graph::Vid> &V) const {
    int hash = (int) V.size();
    for(auto &i: V) {
        hash ^= (int) ((unsigned) i + 0x9e3779b9 + (hash << 6) + (hash >> 2));
    }
    return hash;
}

void NucleusDecomposition::fill_rcs_to_id() {
    for(int i=0; i<(int) rcliques.size(); i++) {
        rcs_to_id.insert({rcliques[i], i});
    }
}

void NucleusDecomposition::combination(int sidx, std::vector<Graph::Vid> &chosen, int s, int r, int idx) {
    if(r==0) {
        std::vector<Graph::Vid> sortedChosen;
        sortedChosen.reserve(chosen.size());
        for(auto &c: chosen) {
            sortedChosen.push_back(c);
        }
        std::sort(sortedChosen.begin(), sortedChosen.end());
        if(!rcs_to_id.count(sortedChosen)) return;
        int id = rcs_to_id[sortedChosen];
        scsHasr[id].insert(sidx);
        rcsIns[sidx].insert(id);
        return;
    }

    if(idx==s) return;
    chosen.push_back(scliques[sidx][idx]);
    combination(sidx, chosen, s, r-1, idx+1);
    chosen.pop_back();
    combination(sidx, chosen, s, r, idx+1);
}

void NucleusDecomposition::findRSCliques(USGraph &graph, int r, int s, bool inadv) { // s should be always larger than r

    auto *rclist = new Kclist();
    auto *sclist = new Kclist();

    rcliques = rclist->getAllKCliques(graph, r, false);
    rclist->rel_idx_and_vid(r_idx_to_vid, r_vid_to_idx);

    fill_rcs_to_id();

    printf("[DEBUG] rcliques: %zu...\n", rcliques.size());

    if(!inadv) return;

    scliques = sclist->getAllKCliques(graph, s, false);
    printf("[DEBUG] scliques: %zu...\n", scliques.size());

    scsHasr.resize(rcliques.size());
    rcsIns.resize(scliques.size());

    for(int i=0; i<(int) scliques.size(); i++) {
        std::vector<Graph::Vid> chosen;
        combination(i, chosen, s, r, 0);
    }
}

void NucleusDecomposition::ndInadv(snu::USGraph &graph, int *max, int r, int s) {
    int fc_t = 0;

    // required for hierarchy
	int cid; // subcore id number
	std::vector<subcore> skeleton; // equal K valued cores
	std::vector<int> component; // subcore ids for each vertex
	std::vector<std::pair<int, int>> relations;
	std::vector<int> unassigned;
	int nSubcores;

    cid = 0;
    nSubcores = 0;
    component.resize(rcliques.size(), -1);

    printf("component resize finished\n");

    K.resize(rcliques.size(), -1);
    Naive_Bucket nBucket;

    int sCliqueMaxNum = 0;
    for(int i=0; i<(int) rcliques.size(); i++) {
        if(sCliqueMaxNum < (int) scsHasr[i].size())
            sCliqueMaxNum = (int) scsHasr[i].size();
    }

    nBucket.Initialize(sCliqueMaxNum, rcliques.size());

    printf("nBucket initialize finished\n");

    for(int i=0; i<(int) rcliques.size(); i++) {
        if(scsHasr[i].empty()) {
            K[i] = 0;
            continue;
        }
        nBucket.Insert(i, scsHasr[i].size());
    }

    printf("nBucket first insertion completed\n");

    while(true) {
        int t;
        int val;
        if(nBucket.PopMin(&t, &val) == -1)
            break;

        unassigned.clear();
        subcore sc (val);
        skeleton.push_back(sc);

        fc_t = K[t] = val;

        for(auto scliq: scsHasr[t]) {
            bool nv = true;
            for(auto rc: rcsIns[scliq]) {
                if(rc==t) continue;
                if(K[rc]!=-1){
                    nv = false;
                    break;
                }
            }

            if(nv) {
                for(auto rc: rcsIns[scliq]) {
                    if(rc==t) continue;
                    if(nBucket.CurrentValue(rc) > fc_t)
                        nBucket.DecVal(rc);
                }
            }
            else {
                rcsIns[scliq].erase(t);
                createSkeleton(t, rcsIns[scliq], &nSubcores, K, skeleton, component, unassigned, relations);
                rcsIns[scliq].insert(t);
            }
        }

        updateUnassigned(t, component, &cid, relations, unassigned);
    }

    printf("[DEBUG] while loop completed...\n");

    nBucket.Free();
    *max = fc_t;

    printf("[DEBUG] K values: %d\n", *max);

    buildHierarchy(*max, relations, skeleton, &nSubcores, graph.E, graph.V);

    printf("[DEBUG] bulidHierarchy Completed...\n");

    printf ("[DEBUG] # subcores: %d\n# subsubcores: %ld\n|V|: %ld\n", nSubcores, skeleton.size(), graph.V);

    presentNuclei(r, s, skeleton, component, graph, graph.E, nd_tree);

    printf("[DEBUG] presentNuclei Completed...\n");
}

void NucleusDecomposition::fill_conn_deg(snu::USGraph &graph) {
    for(auto &v: graph.vertices) {
        degree[r_vid_to_idx[v->id]] = (int)v->edges.size();
        for(auto &w: graph.vertices) {
            if(v==w) continue;
            conn[r_vid_to_idx[v->id]].insert({r_vid_to_idx[w->id], false});
        }
        for(auto &w: v->edges) {
            int from = (int)w->from->id;
            int to = (int)w->to->id;
            if(v->id == from) conn[r_vid_to_idx[v->id]][r_vid_to_idx[to]] = true;
            else conn[r_vid_to_idx[v->id]][r_vid_to_idx[from]] = true;
        }
    }
}

/* fixing... */
void NucleusDecomposition::fill_edge_umap(snu::USGraph &graph) {

    for(auto &v: graph.vertices) {
        degree[v->id] = (int)v->edges.size();
        for(auto &w: v->edges) {
            int from = (int)w->from->id;
            int to = (int)w->to->id;
            if(v->id == from) edge_umap[v->id].insert(to);
            else edge_umap[v->id].insert(from);
        }
    }
}

int NucleusDecomposition::count_scliques(USGraph &graph, std::vector<Graph::Vid> &clique, int r, int s, Graph::Vid maxNew) {

    if(r==s) return 1;
    int ret = 0;
    int smallest = -1, minD = INT32_MAX;
    for(auto &v: clique) {
        if(degree[v] < minD) {
            smallest = v;
            minD = degree[v];
        }
    }

    smallest = r_idx_to_vid[smallest];

    for(auto &v: graph.id_to_vertex[smallest]->edges) {
        int from = (int) v->from->id;
        int to = (int) v->to->id;
        int check = (smallest==from) ? to : from;
        check = r_vid_to_idx[check];

        if(std::find(clique.begin(), clique.end(), check) != clique.end()) continue;
        if(check < maxNew) continue;

        bool canadd = true;
        for(auto &u: clique) {
            if(!conn[check][u]){
                canadd = false;
                break;
            }
//            if(!edge_umap[check].count(u)) {
//                canadd = false;
//                break;
//            }
        }
        if(canadd){
            clique.push_back(check);
            ret += count_scliques(graph, clique, r+1, s, check);
            clique.pop_back();
        }
    }
    return ret;
}

void NucleusDecomposition::get_scliques_rnow(snu::USGraph &graph, std::vector<Graph::Vid> &clique, int r, int s, std::vector<std::vector<Graph::Vid>> &result, Graph::Vid maxNew) {
    if(r==s) {
        std::vector<Graph::Vid> sorted;
        sorted.reserve(clique.size());
        for(auto c: clique) {
            sorted.push_back(c);
        }
        std::sort(sorted.begin(), sorted.end());
        result.push_back(sorted);
        return;
    }
    int smallest = -1, minD = INT32_MAX;
    for(auto &v: clique) {
        if(degree[v] < minD) {
            smallest = v;
            minD = degree[v];
        }
    }
    smallest = r_idx_to_vid[smallest];
    for(auto &v: graph.id_to_vertex[smallest]->edges) {
        int from = (int) v->from->id;
        int to = (int) v->to->id;
        int check = (smallest==from) ? to : from;
        check = r_vid_to_idx[check];

        if(std::find(clique.begin(), clique.end(), check) != clique.end()) continue;
        if(check < maxNew) continue;

        bool canadd = true;
        for(auto &u: clique) {
            if(!conn[check][u]){
                canadd = false;
                break;
            }
//            if(!edge_umap[check].count(u)) {
//                canadd = false;
//                break;
//            }
        }
        if(canadd){
            clique.push_back(check);
            get_scliques_rnow(graph, clique, r+1, s, result, check);
            clique.pop_back();
        }
    }
}

void NucleusDecomposition::find_rc_in_sc(std::vector<Graph::Vid> &sclique, std::vector<Graph::Vid> &chosen, int r, int s, int idx, std::set<int> &result) {
    if(r==0) {
        std::vector<Graph::Vid> sortedChosen;
        sortedChosen.reserve(chosen.size());
        for(auto &c: chosen) {
            sortedChosen.push_back(c);
        }
        std::sort(sortedChosen.begin(), sortedChosen.end());
        if(!rcs_to_id.count(sortedChosen)) return;
        int id = rcs_to_id[sortedChosen];
        result.insert(id);
        return;
    }

    if(idx==s) return;
    chosen.push_back(sclique[idx]);
    find_rc_in_sc(sclique, chosen, r-1, s, idx+1, result);
    chosen.pop_back();
    find_rc_in_sc(sclique, chosen, r, s, idx+1, result);
}

int NucleusDecomposition::calculate_combination(int s, int r) {
    int result = 1;
    for(int i=1; i<=s; i++) result *= i;
    for(int i=1; i<=r; i++) result /= i;
    for(int i=1; i<=(s-r); i++) result /= i;
    return result;
}

void NucleusDecomposition::ndWing(snu::USGraph &graph, int *max, int r, int s) {
    int fc_t = 0;

    // required for hierarchy
	int cid; // subcore id number
	std::vector<subcore> skeleton; // equal K valued cores
	std::vector<int> component; // subcore ids for each vertex
	std::vector<std::pair<int, int>> relations;
	std::vector<int> unassigned;
	int nSubcores;

    cid = 0;
    nSubcores = 0;
    component.resize(rcliques.size(), -1);

    K.resize(rcliques.size(), -1);
    Naive_Bucket nBucket;
    nBucket.Initialize(graph.V * graph.V, rcliques.size());

    printf("[DEBUG] Bucket Initialize Completed...\n");

    fill_conn_deg(graph);
//    fill_edge_umap(graph);

    printf("[DEBUG] Fill Edge Unordered map Completed...\n");

    std::vector<int> sc_count;
    sc_count.reserve(rcliques.size());
    for(auto &rc: rcliques) {
        std::unordered_set<Graph::Vid> save;
        sc_count.push_back(count_scliques(graph, rc, r, s, -1));
    }

    for(int i=0; i<(int) rcliques.size(); i++) {
        if(sc_count[i] == 0) {
            K[i] = 0;
            continue;
        }
        nBucket.Insert(i, sc_count[i]);
    }

    printf("[DEBUG] Bucket Insertion Completed...\n");

    while(true) {
        int t;
        int val;
        if(nBucket.PopMin(&t, &val) == -1)
            break;

        unassigned.clear();
        subcore sc (val);
        skeleton.push_back(sc);

        fc_t = K[t] = val;

        std::vector<std::vector<Graph::Vid>> scs_has_r_now;
        std::unordered_set<Graph::Vid> save;
        get_scliques_rnow(graph, rcliques[t], r, s, scs_has_r_now, -1);

        for(auto scliq: scs_has_r_now) {
            std::set<int> rcs;
            std::vector<Graph::Vid> chosen;
            find_rc_in_sc(scliq, chosen, r, s, 0, rcs);

            bool nv = true;
            for(auto rc: rcs) {
                if(rc==t) continue;
                if(K[rc]!=-1){
                    nv = false;
                    break;
                }
            }

            if(nv) {
                for(auto rc: rcs) {
                    if(rc==t) continue;
                    if(nBucket.CurrentValue(rc) > fc_t)
                        nBucket.DecVal(rc);
                }
            }
            else {
                rcs.erase(t);
                createSkeleton(t, rcs, &nSubcores, K, skeleton, component, unassigned, relations);
                rcs.insert(t);
            }
        }

        updateUnassigned(t, component, &cid, relations, unassigned);
    }

    printf("[DEBUG] while loop completed...\n");

    nBucket.Free();
    *max = fc_t;

    printf("[DEBUG] K values: %d\n", *max);

    buildHierarchy(*max, relations, skeleton, &nSubcores, graph.E, graph.V);

    printf("[DEBUG] bulidHierarchy Completed...\n");

    printf ("[DEBUG] # subcores: %d\n# subsubcores: %ld\n|V|: %ld\n", nSubcores, skeleton.size(), graph.V);

    presentNuclei(r, s, skeleton, component, graph, graph.E, nd_tree);

    printf("[DEBUG] presentNuclei Completed...\n");
}

////////////////////////////////////////////////////////////
//////////////// Build Hierarchy Functions /////////////////
////////////////////////////////////////////////////////////

void NucleusDecomposition::assignToRoot(int* ch, std::vector<subcore>& skeleton) {
    std::vector<int> acc;
	int s = *ch;
	while (skeleton[s].root != -1) {
		acc.push_back (s);
		s = skeleton[s].root;
	}
	for (int i : acc)
		skeleton[i].root = s;
	*ch = s;
}

void NucleusDecomposition::assignToRepresentative(int* ch, std::vector<subcore>& skeleton) {
    int u = *ch;
	std::vector<int> vs;
	while (skeleton[u].parent != -1) {
		int n = skeleton[u].parent;
		if (skeleton[n].K == skeleton[u].K) {
			vs.push_back (u);
			u = n;
		}
		else
			break;
	}
	*ch = u;
	for (int i : vs) {
		if (i != u)
			skeleton[i].parent = u;
	}
}

void NucleusDecomposition::store(int uComp, int vComp, std::vector<int>& unassigned, std::vector<std::pair<int, int> >& relations) {
    std::pair<int, int> c (vComp, uComp);
	if (uComp == -1) // it is possible that u didn't get an id yet
		unassigned.push_back (relations.size()); // keep those indices to process after the loop, below
	relations.push_back (c);
}

inline void NucleusDecomposition::merge(int u, int v, std::vector<int>& component, std::vector<subcore>& skeleton, int* nSubcores) {
    if (component[u] == -1) {
		component[u] = component[v];
		skeleton.erase (skeleton.end() - 1);
	}
	else { // merge component[u] and component[v] nodes
		int child = component[u];
		int parent = component[v];
		assignToRepresentative (&child, skeleton);
		assignToRepresentative (&parent, skeleton);
		if (child != parent) {
			if (skeleton[child].rank > skeleton[parent].rank)
				std::swap (child, parent);
			skeleton[child].parent = parent;
			skeleton[child].visible = false;
			if (skeleton[parent].rank == skeleton[child].rank)
				skeleton[parent].rank++;
			*nSubcores--;
		}
	}
}

void NucleusDecomposition::createSkeleton(int u, const std::set<int> neighbors, int* nSubcores, std::vector<int>& K, std::vector<subcore>& skeleton,
		std::vector<int>& component, std::vector<int>& unassigned, std::vector<std::pair<int, int> >& relations) {
    int smallest = -1, minK = INT32_MAX;
	for (auto i : neighbors)
		if (K[i] != -1 && K[i] < minK) {
			smallest = i;
			minK = K[i];
		}
	if (smallest == -1)
		return;

	if (K[smallest] == K[u])
		merge (u, smallest, component, skeleton, nSubcores);
	else
		store (component[u], component[smallest], unassigned, relations);
}

void NucleusDecomposition::updateUnassigned(int t, std::vector<int>& component, int* cid, std::vector<std::pair<int, int> >& relations, std::vector<int>& unassigned) {
    if (component[t] == -1) { // if e didn't get a component, give her a new one
		component[t] = *cid;
		++(*cid);
	}

	// update the unassigned components that are in the relations
	for (int i : unassigned)
		relations[i] = std::make_pair (relations[i].first, component[t]);
}

void NucleusDecomposition::buildHierarchy(int cn, std::vector<std::pair<int, int> >& relations, std::vector<subcore>& skeleton, int* nSubcores, int nEdge, int nVtx) {
    // bin the relations w.r.t. first's K
	std::vector<std::vector<std::pair<int, int>>> binnedRelations (cn + 1);

	for (auto r : relations) {
		int a = r.first;
		int b = r.second;
		assignToRepresentative (&a, skeleton);
		assignToRepresentative (&b, skeleton);
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
			assignToRoot (&root, skeleton);
			if (a != root) {
				if (skeleton[a].K < skeleton[root].K) {
					skeleton[root].parent = a;
					skeleton[root].root = a;
				}
				else { // skeleton[root].K == skeleton[a].K
					std::pair<int, int> c = (root < a) ? std::make_pair (root, a) : std::make_pair (a, root);
					mergeList.push_back (c);
				}
			}
		}

		// handle merges
		for (auto sc : mergeList) {
			int child = sc.first;
			int parent = sc.second;
			assignToRepresentative (&child, skeleton);
			assignToRepresentative (&parent, skeleton);
			if (child != parent) {
				if (skeleton[child].rank > skeleton[parent].rank)
					std::swap (child, parent);
				skeleton[child].parent = parent;
				skeleton[child].root = parent;
				skeleton[child].visible = false;
				if (skeleton[parent].rank == skeleton[child].rank)
					skeleton[parent].rank++;
			}
		}
	}

	*nSubcores += skeleton.size();

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

////////////////////////////////////////////////////////////
//////////////// Building Nucleus Functions ////////////////
////////////////////////////////////////////////////////////

inline int NucleusDecomposition::commons(std::vector<int>& a, std::vector<int>& b) {
    int i = 0, j = 0;
	int count = 0;
	while (i < a.size() && j < b.size()) {
		if (a[i] < b[j])
			i++;
		else if (b[j] < a[i])
			j++;
		else {
			count++;
			i++;
			j++;
		}
	}
	return count;
}

void NucleusDecomposition::rearrange(std::vector<subcore>& skeleton) {
    for (size_t i = 0; i < skeleton.size(); i++)
		skeleton[i].children.clear();
	for (size_t i = 0; i < skeleton.size() - 1; i++) {
		if (skeleton[i].visible) {
			int pr = skeleton[i].parent;
			while (!skeleton[pr].visible)
				pr = skeleton[pr].parent;
			skeleton[i].parent = pr;
			skeleton[pr].children.push_back (i);
		}
	}
}

void NucleusDecomposition::reportSubgraph(int r, int s, int index, std::vector<int>& component, std::vector<subcore>& skeleton, Graph& graph, int nEdge, std::vector<nd_tree_node> &nd_tree,
					std::unordered_map<int, int> &skeleton_to_nd_tree, std::vector<bool> visited) {
    if (skeleton[index].parent == -1 && skeleton[index].K == 0) {
		skeleton[index].size = graph.V;
		skeleton[index].nEdge = nEdge;
		skeleton[index].ed = 0;
		return;
	}

	std::fill(visited.begin(), visited.end(), false);

	std::vector<int> vset;
	for(int i=0; i<component.size(); i++) {
        if(component[i] == index) {
            for(auto j: rcliques[i]) {
                if(!visited[j]) {
                    vset.push_back(j);
                    visited[j] = true;
                }
            }
        }
    }

	for (auto child : skeleton[index].children) {
		int nd_node_id = skeleton_to_nd_tree[child];

		for (auto u : nd_tree[nd_node_id].vertices_) {
			if (!visited[u]) {
				vset.push_back(u);
				visited[u] = true;
			}
		}
	}

	std::sort(vset.begin(), vset.end());

	// edge density
	int edge_count = 0;
	for (size_t i = 0; i < vset.size(); i++) {
        std::vector<int> tmp;
        for(auto t: graph.id_to_vertex[r_idx_to_vid[vset[i]]]->edges) {
            int from = r_vid_to_idx[t->from->id];
            int to = r_vid_to_idx[t->to->id];
            if(vset[i]==from) tmp.push_back(to);
            else tmp.push_back(from);
        }
        edge_count += commons (vset, tmp);
    }

	edge_count /= 2;

	double density = 0;
	if (vset.size() > 1)
		density = edge_count /(vset.size() * (vset.size() - 1) / 2.0);

	nd_tree_node node;
	node.parent_ = -1;
	node.k_ = skeleton[index].K;
	node.r_ = r;
	node.s_ = s;
	node.num_edges_ = edge_count;
	node.density_ = density;
	node.vertices_ = vset;
	node.id_ = (int)nd_tree.size();
	skeleton_to_nd_tree[index] = node.id_;
	for (auto child : skeleton[index].children) {
		int node_id = skeleton_to_nd_tree[child];
		node.children_.push_back(node_id);
		nd_tree[node_id].parent_ = node.id_;
	}
	nd_tree.push_back(node);
}

void NucleusDecomposition::bfsHierarchy(std::vector<subcore>& skeleton, std::stack<int>& scs) {
    rearrange (skeleton);
	std::queue<int> bfsorder; // we are doing bfs on the hierarchy tree and push the dequeued nodes to the stack
	bfsorder.push(skeleton.size() - 1);
	while (!bfsorder.empty()) {
		int s = bfsorder.front();
		bfsorder.pop();
		scs.push (s);
		for (int r : skeleton[s].children)
			bfsorder.push (r);
	}
}

inline void NucleusDecomposition::findRepresentative(int* child, std::vector<subcore>& skeleton) {
    int u = *child;
	if (skeleton[u].parent != -1) {
		int pr = skeleton[u].parent;
		while (skeleton[u].K == skeleton[pr].K) {
			u = pr;
			if (skeleton[u].parent != -1)
				pr = skeleton[u].parent;
			else
				break;
		}
	}
	*child = u;
}

inline void NucleusDecomposition::removeChild(int i, std::vector<subcore> &backup) {
    if (backup[i].parent != -1) {
		int pr = backup[i].parent;
		for (int j = 0; j < (int) backup[pr].children.size(); j++)
			if (backup[pr].children[j] == i) {
				backup[pr].children.erase (backup[pr].children.begin() + j);
				break;
			}
	}
}

void NucleusDecomposition::presentNuclei(int r, int s, std::vector<subcore>& skeleton, std::vector<int>& component, USGraph& graph, int nEdge, std::vector<nd_tree_node> &nd_tree) {
    // assign unassigned items to top subcore
	for (int i = 0; i < component.size(); i++)
		if (component[i] == -1)
			component[i] = skeleton.size() - 1;

	// match each component with its representative
	std::unordered_map<int, int> replace;
	for (int i = 0; i < skeleton.size(); i++) {
		int sc = i;
		int original = sc;
		findRepresentative (&sc, skeleton);
		if (original != sc)
			skeleton[original].visible = false;
		replace[original] = sc;
	}

	// each component takes its representative's component number
	for (int i = 0; i < component.size(); i++)
		if (!replace.count(component[i]))
			component[i] = replace[component[i]];

	std::stack<int> subcoreStack;
	bfsHierarchy (skeleton, subcoreStack);

	std::vector<bool> visited(graph.V);
	std::unordered_map<int, int> skeleton_to_nd_tree;

	while (!subcoreStack.empty()) {
		int i = subcoreStack.top();
		subcoreStack.pop();
		if (skeleton[i].visible) {
			reportSubgraph (r, s, i, component, skeleton, graph, nEdge, nd_tree, skeleton_to_nd_tree, visited);
		}
	}

	printf("[DEBUG] skeleton length: %d\n", skeleton.size());
}

} // namespace snu
