#include "mpi.h"
#include <vector>
#include "xolotlMemUsage/summaryproc/Statm/MemUsageStats.h"

namespace xolotlMemUsage {

void
MemUsageStats::Aggregate(std::shared_ptr<MemUsageStats> localValue,
                    MPI_Comm aggComm,
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

