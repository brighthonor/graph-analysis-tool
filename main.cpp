#include <getopt.h>
#include <stdio.h>

#include <chrono>

#include "snugal.h"
#include "statanalyzer.h"

template <typename T>
static void printElapsedTime(const char *msg, T &t) {
    using namespace std::chrono;
    printf("%10ldms: %s\n", (duration_cast<milliseconds>(high_resolution_clock::now() - t)).count(), msg);
    t = high_resolution_clock::now();
}

void usage(void);

int main(int argc, char *argv[]) {
    if (argc < 2) {
        usage();
        return 1;
    }

    bool directed = false;     // suppose given graph is undirected.
    bool file_output = false;  // default: skip output of detailed files
    bool display_time = true;  // always: show time statistics
    bool verify = false;       // verify results (default: skip verification)

    while (true) {
        static struct option long_options[] = {
            {"directed", 0, NULL, 'd'},
            {"undirected", 0, NULL, 'u'},
            {"help", 0, NULL, 'h'},
            {"file_output", 0, NULL, 'f'},
            {"verify", 0, NULL, 'v'},
            {0, 0, 0, 0}
        };

        int option = getopt_long(argc, argv, "duhfv", long_options, NULL);

        if (option < 0) {
            break;
        }

        switch (option) {
            case 'd':
                directed = true;
                break;

            case 'u':
                directed = false;
                break;

            case 'f':
                file_output = true;
                break;

            case 'v':
                verify = true;
                break;

            case 'h':
                usage();
                return 0;

            default:
                usage();
                return 1;
        }
    }

    // There exists no non-option argument.
    // That is, input file path is not given.
    if (argc <= optind) {
        fprintf(stderr, "input file is not given.\n");
        return 1;
    }

    // Suppose input file path is the first non-option argument.
    std::string input_path(argv[optind]);

    // Extract name of given graph from input path.
    std::size_t slash = input_path.rfind("/");
    if (slash == std::string::npos) {
        slash = -1;
    }

    std::size_t dot = input_path.rfind(".");
    if (dot == std::string::npos) {
        dot = input_path.length();
    }

    std::string graph_name = input_path.substr(slash + 1, dot - (slash + 1));

    // parse graph from file
    std::shared_ptr<snu::Graph> graph_shptr = snu::parseFile(input_path, directed);
    if (!graph_shptr)
        return 1;
    snu::Graph &graph = *graph_shptr.get();

    std::vector<snu::Stat *> all_stats;
    std::vector<snu::CommonStat *> common_stats;
    std::vector<snu::DirectedStat *> directed_only_stats;
    std::vector<snu::UndirectedStat *> undirected_only_stats;

    /* Add Stats Here */
    auto basicStat = snu::BasicStat();
    auto connectStat = snu::ConnectStat();
    auto closenessCentrality = snu::ClosenessCentrality();
    auto betweennessCentrality = snu::BetweennessCentrality();
    auto eigenCentrality = snu::EigenCentrality();

    auto countStat = snu::CountStat();
    auto biconnectedComponents = snu::BiconnectedComponents();

    auto kclist = snu::Kclist();
    auto nucleusDecomposition = snu::NucleusDecomposition();

    /* -- Add Common Stats To Run (support both directed & undirected) -- */
    common_stats = {
        &basicStat,
//        &connectStat,
//        &eigenCentrality,
//        &closenessCentrality,
//        &betweennessCentrality,
    };

    /* -- Add Directed Stats To Run -- */
    directed_only_stats = {};

    /* -- Add Undirected Stats To Run -- */
    undirected_only_stats = {
//        &countStat,
//        &biconnectedComponents,
        &kclist,
        &nucleusDecomposition,
    };

    // run all stats using the graph
    auto analyzer = snu::StatAnalyzer(graph_name, graph_shptr, file_output, display_time, verify);
    analyzer.addCommonStats(common_stats);
    analyzer.addDirectedStats(directed_only_stats);
    analyzer.addUndirectedStats(undirected_only_stats);
    analyzer.run();

    // plot
    auto t = std::chrono::high_resolution_clock::now();
    snu::Plot plot(graph_name);
    snu::makePlot(graph, plot);
    printElapsedTime("makePlot()", t);

    // write to HTML
    FILE *fp = snu::openHtml(graph_name.c_str(), directed);
    std::vector<snu::Stat *> html_order;
    if (directed) {
        /* -- Set Directed Html Order Here -- */
        html_order = {
            &basicStat,
            &connectStat,
            &eigenCentrality,
            &closenessCentrality,
            &betweennessCentrality,
        };
    } else {
        /* -- Set Undirected Html Order Here -- */
        html_order = {
            &basicStat,
//            &connectStat,
//            &eigenCentrality,
//            &biconnectedComponents,
//            &closenessCentrality,
//            &betweennessCentrality,
//            &countStat,
            &kclist,
            &nucleusDecomposition,
        };
    }
    for (auto stat : html_order) {
        snu::addStatToHtml(fp, stat, directed);
    }
    snu::closeHtml(fp, plot);
    printElapsedTime("Make HTML", t);

    return 0;
}

void usage(void) {
    printf(" usage: main <input file> [options]\n");
    printf("\n");
    printf("  options:\n");
    printf("   -d | --directed       set graph directed\n");
    printf("   -u | --undirected     set graph undirected (default)\n");
    printf("   -f | --file_output    output centrality details as separate files\n");
    printf("   -v | --verify         verify correctness of results (slower)\n");
    printf("   -h | --help           print this list of help\n");
    printf("\n");
    printf("   e.g. ./main some_path/some_graph.snap --directed\n");

    return;
}
