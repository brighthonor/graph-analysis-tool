#ifndef STAT_H
#define STAT_H

namespace snu {

	typedef struct StatResult {
		long long size; // number of vertices
		long long volume; // number of edges
		long long num_vlabel; // number of vertex label
		long long num_elabel; // number of edge label
		double density; // m / n^2, directed graph with loops
		double avg_degree; // average degree, 2m / n
		long long max_degree; // maximum degree
		long long max_indegree; // maximum indegree
		long long max_outdegree; // maximum outdegree
		// size of vertex label is defined as number of vertices having that label
		long long max_vlabel_size; // maximum size of vertex label
		long long max_elabel_size; // maximum size of edge label
		double reciprocity; // only defined in directed graph
		double negativity; // ratio of nagative edges

		// connectivity statistics 
		long long num_cc; // number of connected components
		long long num_scc; // number of strongly connected components
		long long num_wcc; // number of weakly connected components
		long long size_lcc; // largest size of connected components
		long long size_lscc; // largest size of strongly connected components
		long long size_lwcc; // largest size of weakly connected components
		long long diameter; // longest shortest paths

        // eigencentrality
        bool pagerank_converged = false; // whether pagerank converged
        double max_pagerank; // maximum pagerank value
        long long max_pagerank_id; // id of max pagerank value vertex
        bool eigencentrality_converged = false; // whether eigencentrality converged
        double max_eigencentrality; // maximum eigencentrality value
        long long max_eigencentrality_id; // id of max eigencentrality vertex
        bool katz_centrality_computed = false; // whether katz centrality was computed
        double max_katz_centrality; // maximum katz centrality value
        long long max_katz_centrality_id; // id of max katz centrality vertex

        // biconnected components and articulation points
        long long num_arp; // number of articulation points
        long long num_bcc; // number of biconnected components
        long long max_conn_bcc; // maximum number of bcc's connected to a single arp
        long long size_lbcc; // largest size of biconnected components

		// closeness centrality
		double max_closeness_centrality; // maximum closeness centrality value
		long long max_closeness_centrality_id; // id of max closeness centrality vertex

		// betweenness centrality
		double max_betweenness_centrality; // maximum betweenness centrality value
		long long max_betweenness_centrality_id; // id of max betweenness centrality vertex

		// counting statistics 
		unsigned long long wedge_count; // number of wedges
		unsigned long long claw_count; // number of claws
		unsigned long long triangle_count; // number of triangles

		// statistic check
		bool basicstat = false;
		bool connectstat = false;
        bool eigencentralitystat = false;
        bool biconnectedstat = false;
		bool countstat = false;
		bool closenessstat = false;
		bool betweennessstat = false;
	} StatResult;

	void initStat(StatResult& result);
}

#endif // STAT_H
