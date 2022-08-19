#include "basicstat.h"

namespace snu {
std::string BasicStat::statName() {
    return "BasicStat";
}

bool BasicStat::calculateStat(Graph &graph, bool verify) {
    const DSGraph *dsgraph = dynamic_cast<const DSGraph *>(&graph);
    const USGraph *usgraph = dynamic_cast<const USGraph *>(&graph);

    if (dsgraph)
        return calcDS(*dsgraph);
    else if (usgraph)
        return calcUS(*usgraph);
    return false;
}

bool BasicStat::calcDS(const DSGraph &graph) {
    long long n = graph.id_to_vertex.size();
    long long m = graph.id_to_edge.size();

    size = n;
    volume = m;
    num_vlabel = graph.vlabel_to_class.size();
    num_elabel = graph.elabel_to_class.size();
    density = (double)m / (n * n);  // assume there is loop

    max_indegree = 0;
    max_outdegree = 0;
    for (auto it = graph.id_to_vertex.begin(); it != graph.id_to_vertex.end(); it++) {
        max_indegree = std::max(max_indegree, it->second->indegree);
        max_outdegree = std::max(max_outdegree, (long long)it->second->edges.size());
    }
    max_indegree = max_indegree;
    max_outdegree = max_outdegree;

    max_vlabel_size = 0;
    for (auto it = graph.vlabel_to_class.begin(); it != graph.vlabel_to_class.end(); it++)
        max_vlabel_size = std::max(max_vlabel_size, (long long)it->second->vertices.size());
    max_vlabel_size = max_vlabel_size;

    max_elabel_size = 0;
    for (auto it = graph.elabel_to_class.begin(); it != graph.elabel_to_class.end(); it++)
        max_elabel_size = std::max(max_elabel_size, (long long)it->second->edges.size());
    max_elabel_size = max_elabel_size;

    long long cycle_count = 0;
    for (auto it = graph.id_to_edge.begin(); it != graph.id_to_edge.end(); it++) {
        Graph::Vid from = it->second->from->id;
        Graph::Vid to = it->second->to->id;
        if (graph.is_connected.count(std::make_pair(to, from)))
            cycle_count++;
    }
    reciprocity = (double)cycle_count / m;

    negativity = (double)graph.negative_edge_num / m;
    return true;
}

bool BasicStat::calcUS(const USGraph &graph) {
    long long n = graph.id_to_vertex.size();
    long long m = graph.id_to_edge.size();

    size = n;    // n
    volume = m;  // m
    num_vlabel = graph.vlabel_to_class.size();
    num_elabel = graph.elabel_to_class.size();
    density = 2.0 * m / (n * (n + 1));  // assume there is loop
    avg_degree = 2.0 * m / n;           // 2 * m / n

    max_degree = 0;
    for (auto it = graph.id_to_vertex.begin(); it != graph.id_to_vertex.end(); it++)
        max_degree = std::max(max_degree, it->second->indegree);  // indegree == outdegree
    max_degree = max_degree;

    max_vlabel_size = 0;
    for (auto it = graph.vlabel_to_class.begin(); it != graph.vlabel_to_class.end(); it++)
        max_vlabel_size = std::max(max_vlabel_size, (long long)it->second->vertices.size());
    max_vlabel_size = max_vlabel_size;

    max_elabel_size = 0;
    for (auto it = graph.elabel_to_class.begin(); it != graph.elabel_to_class.end(); it++)
        max_elabel_size = std::max(max_elabel_size, (long long)it->second->edges.size());
    max_elabel_size = max_elabel_size;

    negativity = (double)graph.negative_edge_num / m;

    return true;
}

void BasicStat::writeToHTMLStat(FILE *fp, bool directed) {
    if (directed)
        htmlDS(fp);
    else
        htmlUS(fp);
}

void BasicStat::htmlDS(FILE *fp) {
    fprintf(fp,
            "\
				<h2>\
						Basic Statistics\
				</h2>\
				<h3>\
						<p> size = %lld </p>\
						<p> volume = %lld </p>\
						<p> number of vertex labels = %lld </p>\
						<p> number of edge labels = %lld </p>\
						<p> density = %lf </p>\
						<p> maximum indegree = %lld </p>\
						<p> maximum outdegree = %lld </p>\
						<p> maximum size of vertex label = %lld </p>\
						<p> maximum size of edge label = %lld </p>\
						<p> reciprocity = %lf </p>\
						<p> negativity = %lf </p>\
				</h3>\
		",
            size, volume, num_vlabel, num_elabel, density, max_indegree, max_outdegree, max_vlabel_size, max_elabel_size, reciprocity, negativity);
}
void BasicStat::htmlUS(FILE *fp) {
    fprintf(fp,
            "\
				<h2>\
						Basic Statistics\
				</h2>\
				<h3>\
						<p> size = %lld </p>\
						<p> volume = %lld </p>\
						<p> number of vertex labels = %lld </p>\
						<p> number of edge labels = %lld </p>\
						<p> density = %lf </p>\
						<p> average degree = %lf </p>\
						<p> maximum degree = %lld </p>\
						<p> maximum size of vertex label = %lld </p>\
						<p> maximum size of edge label = %lld </p>\
						<p> negativity = %lf </p>\
				</h3>\
		",
            size, volume, num_vlabel, num_elabel, density, avg_degree, max_degree, max_vlabel_size, max_elabel_size, negativity);
}
}  // namespace snu
