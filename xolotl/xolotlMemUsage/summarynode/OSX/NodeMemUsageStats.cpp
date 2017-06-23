#include "mpi.h"
#include <vector>
#include "xolotlMemUsage/summarynode/OSX/NodeMemUsageStats.h"

namespace xolotlMemUsage {

void
NodeMemUsageStats::Aggregate(std::shared_ptr<NodeMemUsageStats> localValue,
                    MPI_Comm aggComm,
                    int myRank,
                    bool localIsValid,
                    uint32_t processCount) {

    std::vector<NodeMemUsageMetricStats NodeMemUsageStats::*> allMetricsToAggregate{
            &NodeMemUsageStats::free_count,
            &NodeMemUsageStats::active_count,
            &NodeMemUsageStats::inactive_count,
            &NodeMemUsageStats::wire_count,
            &NodeMemUsageStats::zero_fill_count,
            &NodeMemUsageStats::reactivations,
            &NodeMemUsageStats::pageins,
            &NodeMemUsageStats::pageouts,
            &NodeMemUsageStats::faults,
            &NodeMemUsageStats::cow_faults,
            &NodeMemUsageStats::lookups,
            &NodeMemUsageStats::hits,
            &NodeMemUsageStats::purges,
            &NodeMemUsageStats::purgeable_count
        };

    // Find global statistics about all metrics we know about.
    for(auto pCurrMetric : allMetricsToAggregate) {

        MinReduceKnown<uint64_t>((localValue.get()->*pCurrMetric).min,
                                    &((this->*pCurrMetric).min),
                                    aggComm,
                                    myRank,
                                    MPI_UNSIGNED_LONG,
                                    localIsValid);
        MaxReduceKnown<uint64_t>((localValue.get()->*pCurrMetric).max,
                                    &((this->*pCurrMetric).max),
                                    aggComm,
                                    myRank,
                                    MPI_UNSIGNED_LONG,
                                    localIsValid);

        AvgReduceKnown<double>((localValue.get()->*pCurrMetric).average,
                                    &((this->*pCurrMetric).average),
                                    aggComm,
                                    myRank,
                                    MPI_DOUBLE,
                                    localIsValid,
                                    processCount);
        StdevReduceKnown<double>((localValue.get()->*pCurrMetric).average,
                                    (this->*pCurrMetric).average,
                                    &((this->*pCurrMetric).stdev),
                                    aggComm,
                                    myRank,
                                    MPI_DOUBLE,
                                    localIsValid,
                                    processCount);
    }
}

} // namespace xolotlMemUsage

