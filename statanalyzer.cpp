#include "statanalyzer.h"

#include <getopt.h>
#include <stdio.h>

#include <chrono>

#include "parse.h"

using namespace std::chrono;

namespace snu {
template <typename T>
static void printElapsedTime(const char *msg, T &t) {
    printf("%10ldms: %s\n", (duration_cast<milliseconds>(high_resolution_clock::now() - t)).count(), msg);
    t = high_resolution_clock::now();
}

StatAnalyzer::StatAnalyzer(std::string _graph_name,
                           std::shared_ptr<Graph> _graph_shptr,
                           bool _file_output,
                           bool _display_time,
                           bool _verify)
    : graph_name(_graph_name), graph_shptr(_graph_shptr), file_output(_file_output), display_time(_display_time), verify(_verify) {
}

bool StatAnalyzer::run() {
    std::vector<snu::Stat *> all_stats;

    USGraph *usgraph_ptr = dynamic_cast<USGraph *>(graph_shptr.get());
    DSGraph *dsgraph_ptr = dynamic_cast<DSGraph *>(graph_shptr.get());
    bool directed = dsgraph_ptr != nullptr;

    auto &graph = *graph_shptr.get();

    for (auto p : common_stats)
        all_stats.push_back(p);
    if (directed) {
        for (auto p : directed_only_stats)
            all_stats.push_back(p);
    } else {
        for (auto p : undirected_only_stats)
            all_stats.push_back(p);
    }

    auto t = high_resolution_clock::now();

    // calculate common stats
    for (auto commonStat : common_stats) {
        if (!commonStat->calculate(graph, verify)) {
            printf("%s failed\n", commonStat->statName().c_str());
        }
        printElapsedTime(commonStat->statName().c_str(), t);
    }

    // calculate directed stats
    if (directed) {
        for (auto dirStat : directed_only_stats) {
            if (!dirStat->calculateDirected(*dsgraph_ptr, verify)) {
                printf("%s failed\n", dirStat->statName().c_str());
            }
            printElapsedTime(dirStat->statName().c_str(), t);
        }
    }
    // calculate undirected stats
    else {
        for (auto undirStat : undirected_only_stats) {
            if (!undirStat->calculateUndirected(*usgraph_ptr, verify)) {
                printf("%s failed\n", undirStat->statName().c_str());
            }
            printElapsedTime(undirStat->statName().c_str(), t);
        }
    }

    // file output
    if (file_output) {
        for (auto stat : all_stats) {
            if (stat->getSuccess() && stat->writeToFile(graph_name, directed)) {
                printElapsedTime(("Writing to file: " + stat->statName()).c_str(), t);
            }
        }
    }

    return true;
}

void StatAnalyzer::addCommonStats(std::vector<CommonStat *> vec) {
    common_stats = vec;
}

void StatAnalyzer::addDirectedStats(std::vector<DirectedStat *> vec) {
    directed_only_stats = vec;
}

void StatAnalyzer::addUndirectedStats(std::vector<UndirectedStat *> vec) {
    undirected_only_stats = vec;
}
}  // namespace snu
