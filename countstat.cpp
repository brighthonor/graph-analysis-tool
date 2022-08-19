#include "countstat.h"

#include <algorithm>

namespace snu {
std::string CountStat::statName() {
    return "CountStat";
}

bool CountStat::calculateUndirectedStat(USGraph &graph, bool verify) {
    unsigned long long s = 0;
    unsigned long long z = 0;

    std::vector<std::pair<long long, Graph::Vertex *>> vec;

    for (auto it = graph.id_to_vertex.begin(); it != graph.id_to_vertex.end(); it++) {
        long long degree = it->second->indegree;
        s += (unsigned long long)degree * (degree - 1) / 2;
        z += (unsigned long long)degree * (degree - 1) * (degree - 2) / 6;
        vec.push_back(std::make_pair(degree, it->second));
    }

    wedge_count = s;
    claw_count = z;
    triangle_count = TriangleCount(vec);
    return true;
}

void CountStat::writeToHTMLStat(FILE *fp, bool directed) {
    fprintf(fp,
            "\
			<h2>\
					Count statistics\
			</h2>\
			<h3>\
					<p> number of wedges = %llu </p>\
					<p> number of claws = %llu </p>\
					<p> number of triangles = %llu </p>\
			</h3>\
		",
            wedge_count, claw_count, triangle_count);
}

// O(m^(3/2)) time implementation, where m is the number of edges
// http://theory.stanford.edu/~tim/s14/l/l1.pdf
long long CountStat::TriangleCount(std::vector<std::pair<long long, Graph::Vertex *>> &vec) {
    unsigned long long t = 0;

    std::sort(vec.begin(), vec.end());

    int size = vec.size();
    for (int i = 0; i < size; i++)
        vec[i].second->temp = new int(i);  // rank

    bool *check = new bool[size];
    for (int i = 0; i < size; ++i)
        check[i] = false;

    for (int i = 0; i < size; i++) {
        Graph::Vertex *v = vec[i].second;
        int vrank = *(int *)v->temp;

        for (auto it = v->edges.begin(); it != v->edges.end(); it++) {
            Graph::Vertex *u = v == (*it)->to ? (*it)->from : (*it)->to;
            int urank = *(int *)u->temp;
            check[urank] = 1;
        }

        for (auto it = v->edges.begin(); it != v->edges.end(); it++) {
            Graph::Vertex *u = v == (*it)->to ? (*it)->from : (*it)->to;
            int urank = *(int *)u->temp;
            if (vrank < urank)
                for (auto it2 = u->edges.begin(); it2 != u->edges.end(); it2++) {
                    Graph::Vertex *w = u == (*it2)->to ? (*it2)->from : (*it2)->to;
                    int wrank = *(int *)w->temp;
                    if (check[wrank])
                        t++;
                }
        }

        for (auto it = v->edges.begin(); it != v->edges.end(); it++) {
            Graph::Vertex *u = v == (*it)->to ? (*it)->from : (*it)->to;
            int urank = *(int *)u->temp;
            check[urank] = 0;
        }
    }

    for (int i = 0; i < size; i++)
        delete (int *)vec[i].second->temp;

    delete[] check;

    return t / 3;
}
}  // namespace snu
