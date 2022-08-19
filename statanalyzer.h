#ifndef STATANALYZER_H
#define STATANALYZER_H

#include <memory>
#include <vector>

#include "stat.h"

/*
 * Calculate stats that were added
 */

namespace snu {
class StatAnalyzer {
public:
    StatAnalyzer(std::string _graph_name,
                 std::shared_ptr<Graph> _graph_shptr,
                 bool _file_output,
                 bool _display_time,
                 bool _verify);
    bool run();
    void addCommonStats(std::vector<CommonStat *> vec);
    void addDirectedStats(std::vector<DirectedStat *> vec);
    void addUndirectedStats(std::vector<UndirectedStat *> vec);

private:
    std::vector<CommonStat *> common_stats;
    std::vector<snu::DirectedStat *> directed_only_stats;
    std::vector<snu::UndirectedStat *> undirected_only_stats;

    std::string graph_name;
    std::shared_ptr<Graph> graph_shptr;
    bool file_output;
    bool display_time;
    bool verify;
};
}  // namespace snu

#endif  // STATANALYZER_H
