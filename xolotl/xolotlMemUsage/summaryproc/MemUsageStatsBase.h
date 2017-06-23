#ifndef XMEMUSAGE_MEM_STATS_BASE_H
#define XMEMUSAGE_MEM_STATS_BASE_H

#include "mpi.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <limits>
#include "xolotlMemUsage/IMemUsageSampler.h"
#include "xolotlMemUsage/memUsageConfig.h"

#if defined(HAVE_STATM)
    #include "xolotlMemUsage/summaryproc/Statm/StatmSampler.h"
    namespace PerProcDataSource = xolotlMemUsage::Statm;
#else
    #error "Configuration error: thought we had a per-proc data source, but no recognized source available."
#endif // defined(HAVE_STATM)


namespace xolotlMemUsage {

struct MemUsagePageStats : public PerProcDataSource::PageStats {

    double stdev;

    /**
     * Construct a MemUsagePageStats with useful initial values.
     */
    MemUsagePageStats(void)
      : stdev(std::numeric_limits<double>::quiet_NaN())
    {
        // Nothing else to do.
    }


    /**
     * Construct a MemUsagePageStats from existing valid data.
     */
    MemUsagePageStats(uint64_t _min,
                    uint64_t _max,
                    double _avg,
                    double _stdev = std::numeric_limits<double>::quiet_NaN())
      : PerProcDataSource::PageStats(_min, _max, _avg),
        stdev(_stdev)
    {
        // Nothing else to do.
    }


    /**
     * Construct a MemUsagePageStats from existing valid data.
     */
    MemUsagePageStats(const PerProcDataSource::RunningPageStats& _rps, uint64_t _nSamples)
      : PerProcDataSource::PageStats(_rps, _nSamples)
    {
        // Nothing else to do.
    }


    void OutputTo(std::ostream& os, std::string prefix = "") const {

        os << prefix 
            << "min: " << min 
            << " max: " << max 
            << " avg: " << average
            << " stdev: " << stdev;
    }

};

struct MemUsageStatsBase : public IMemUsageSampler::MemUsageData {


    /**
     * Construct a MemUsageStatsBase with useful initial values.
     */
    MemUsageStatsBase(void) = default;


    template<typename T>
    void MinReduceKnown(const T& myValue,
                    T* suggestedTarget,
                    MPI_Comm aggComm,
                    int myRank,
                    MPI_Datatype mpiDatatype,
                    bool knowObject) {

        T* target = (myRank == 0) ? suggestedTarget : nullptr;
        T myCandidate = knowObject ? myValue : std::numeric_limits<T>::max();
        
        MPI_Reduce(&myCandidate, target, 1, mpiDatatype, MPI_MIN, 0, aggComm);
    }

    template<typename T>
    void MaxReduceKnown(const T& myValue,
                    T* suggestedTarget,
                    MPI_Comm aggComm,
                    int myRank,
                    MPI_Datatype mpiDatatype,
                    bool knowObject) {
        
        T* target = (myRank == 0) ? suggestedTarget : nullptr;
        T myCandidate = knowObject ? myValue : std::numeric_limits<T>::min();
        MPI_Reduce(&myCandidate, target, 1, mpiDatatype, MPI_MAX, 0, aggComm);
    }

    template<typename T>
    void AvgReduceKnown(const T& myValue,
                    T* suggestedTarget,
                    MPI_Comm aggComm,
                    int myRank,
                    MPI_Datatype mpiDatatype,
                    bool knowObject,
                    uint32_t processCount) {

        T sum;
        T myContribution = knowObject ? myValue : 0;
        MPI_Reduce(&myContribution, &sum, 1, mpiDatatype, MPI_SUM, 0, aggComm);

        if(myRank == 0) {
            *suggestedTarget = sum / processCount;
        }
    }

    template<typename T>
    void StdevReduceKnown(const T& myValueAverage,
                            const T& aggValueAverage,
                            T* suggestedTarget,
                            MPI_Comm aggComm,
                            int myRank,
                            MPI_Datatype mpiDatatype,
                            bool knowObject,
                            uint32_t processCount) {

        T valSquaredSum;
        T myContribution = knowObject ? (myValueAverage * myValueAverage) : 0;
        MPI_Reduce(&myContribution, &valSquaredSum, 1, mpiDatatype, MPI_SUM, 0, aggComm);

        if(myRank == 0) {
            *suggestedTarget =
                    std::sqrt((valSquaredSum / processCount)
                        - (aggValueAverage * aggValueAverage));
        }
    }
};

} // end namespace xolotlMemUsage

#endif // XMEMUSAGE_MEM_STATS_BASE_H
