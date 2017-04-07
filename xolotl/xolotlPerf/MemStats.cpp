#include "mpi.h"
#include "MemStats.h"

namespace xolotlPerf {

// TODO This is very ugly.  Surely there's a better way.
void
MemStats::Aggregate(const MemStats& localValue,
                    int myRank,
                    bool localIsValid,
                    uint32_t processCount) {

    // Find global minimums for each metric we know about.
    MinReduceKnown<uint64_t>(&MemStats::vmSize, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
    MinReduceKnown<uint64_t>(&MemStats::vmRSS, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
    MinReduceKnown<uint64_t>(&MemStats::rss, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
    MinReduceKnown<uint64_t>(&MemStats::text, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
    MinReduceKnown<uint64_t>(&MemStats::dataAndStack, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
    
    // Find global maximums for each metric we know about.
    MaxReduceKnown<uint64_t>(&MemStats::vmSize, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
    MaxReduceKnown<uint64_t>(&MemStats::vmRSS, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
    MaxReduceKnown<uint64_t>(&MemStats::rss, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
    MaxReduceKnown<uint64_t>(&MemStats::text, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);
    MaxReduceKnown<uint64_t>(&MemStats::dataAndStack, localValue, myRank, MPI_UNSIGNED_LONG, localIsValid);

    // Find global "average of averages" for each metric we know about.
    // Assumes our processCount value has already been set.
    AvgReduceKnown<double>(&MemStats::vmSize, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
    AvgReduceKnown<double>(&MemStats::vmRSS, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
    AvgReduceKnown<double>(&MemStats::rss, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
    AvgReduceKnown<double>(&MemStats::text, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
    AvgReduceKnown<double>(&MemStats::dataAndStack, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);

    // Find global "standard deviation of averages" for each metric we know about.
    StdevReduceKnown<double>(&MemStats::vmSize, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
    StdevReduceKnown<double>(&MemStats::vmRSS, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
    StdevReduceKnown<double>(&MemStats::rss, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
    StdevReduceKnown<double>(&MemStats::text, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
    StdevReduceKnown<double>(&MemStats::dataAndStack, localValue, myRank, MPI_DOUBLE, localIsValid, processCount);
}

} // namespace xolotlPerf

