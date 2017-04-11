#include "mpi.h"
#include "xolotlMemUsage/standard/MemUsageStats.h"

namespace xolotlMemUsage {

// TODO This is very ugly.  Surely there's a better way.
void
MemUsageStats::Aggregate(std::shared_ptr<MemUsageStats> localValue,
                    int myRank,
                    bool localIsValid,
                    uint32_t processCount) {

    // Find global minimums for each metric we know about.
    MinReduceKnown<uint64_t>(&MemUsageStats::vmSize, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
    MinReduceKnown<uint64_t>(&MemUsageStats::vmRSS, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
    MinReduceKnown<uint64_t>(&MemUsageStats::rss, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
    MinReduceKnown<uint64_t>(&MemUsageStats::text, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
    MinReduceKnown<uint64_t>(&MemUsageStats::dataAndStack, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
    
    // Find global maximums for each metric we know about.
    MaxReduceKnown<uint64_t>(&MemUsageStats::vmSize, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
    MaxReduceKnown<uint64_t>(&MemUsageStats::vmRSS, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
    MaxReduceKnown<uint64_t>(&MemUsageStats::rss, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
    MaxReduceKnown<uint64_t>(&MemUsageStats::text, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
    MaxReduceKnown<uint64_t>(&MemUsageStats::dataAndStack, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);

    // Find global "average of averages" for each metric we know about.
    // Assumes our processCount value has already been set.
    AvgReduceKnown<double>(&MemUsageStats::vmSize, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
    AvgReduceKnown<double>(&MemUsageStats::vmRSS, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
    AvgReduceKnown<double>(&MemUsageStats::rss, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
    AvgReduceKnown<double>(&MemUsageStats::text, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
    AvgReduceKnown<double>(&MemUsageStats::dataAndStack, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);

    // Find global "standard deviation of averages" for each metric we know about.
    StdevReduceKnown<double>(&MemUsageStats::vmSize, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
    StdevReduceKnown<double>(&MemUsageStats::vmRSS, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
    StdevReduceKnown<double>(&MemUsageStats::rss, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
    StdevReduceKnown<double>(&MemUsageStats::text, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
    StdevReduceKnown<double>(&MemUsageStats::dataAndStack, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
}

} // namespace xolotlMemUsage

