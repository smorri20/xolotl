#include "mpi.h"
#include <vector>
#include "xolotlMemUsage/summaryproc/MemUsageStats.h"

namespace xolotlMemUsage {

void
MemUsageStats::Aggregate(std::shared_ptr<MemUsageStats> localValue,
                    int myRank,
                    bool localIsValid,
                    uint32_t processCount) {

    std::vector<MemUsagePageStats MemUsageStats::*> allMetricsToAggregate{
        &MemUsageStats::vmSize,
        &MemUsageStats::vmRSS,
        &MemUsageStats::rss,
        &MemUsageStats::text,
        &MemUsageStats::dataAndStack};

    // Find global statistics about all metrics we know about.
    for(auto pCurrMetric : allMetricsToAggregate) {
        MinReduceKnown<uint64_t>(pCurrMetric, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
        MaxReduceKnown<uint64_t>(pCurrMetric, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
        AvgReduceKnown<double>(pCurrMetric, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
        StdevReduceKnown<double>(pCurrMetric, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
    }
}

} // namespace xolotlMemUsage

