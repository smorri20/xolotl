#ifndef XMEMUSAGE_MEM_STATS_H
#define XMEMUSAGE_MEM_STATS_H

#include "mpi.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <limits>
#include "xolotlMemUsage/IMemUsageSampler.h"
#include "xolotlMemUsage/summaryproc/StatmSampler.h"


namespace xolotlMemUsage {

struct MemUsagePageStats : public Statm::PageStats {

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
      : Statm::PageStats(_min, _max, _avg),
        stdev(_stdev)
    {
        // Nothing else to do.
    }


    /**
     * Construct a MemUsagePageStats from existing valid data.
     */
    MemUsagePageStats(const Statm::RunningPageStats& _rps, uint64_t _nSamples)
      : Statm::PageStats(_rps, _nSamples)
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

struct MemUsageStats : public IMemUsageSampler::MemUsageData {
    MemUsagePageStats vmSize;
    MemUsagePageStats vmRSS;
    MemUsagePageStats rss;
    MemUsagePageStats text;
    MemUsagePageStats dataAndStack;


    /**
     * Construct a MemUsageStats with useful initial values.
     */
    MemUsageStats(void) = default;


    /**
     * Construct a MemUsageStats from existing data.
     */
    MemUsageStats(const MemUsagePageStats& _vmSize,
                const MemUsagePageStats& _vmRSS,
                const MemUsagePageStats& _rss,
                const MemUsagePageStats& _text,
                const MemUsagePageStats& _dataAndStack)
      : vmSize(_vmSize),
        vmRSS(_vmRSS),
        rss(_rss),
        text(_text),
        dataAndStack(_dataAndStack)
    {
        // Nothing else to do.
    }



    void OutputTo(std::ostream& os, std::string prefix = "") const {
        // TODO should we just define an ostream inserter for the MemUsagePageStats?

        os << std::fixed << std::setprecision(2);
        vmSize.OutputTo(os, prefix + "vmSize: ");
        os << '\n';
        vmRSS.OutputTo(os, prefix + "vmRSS: ");
        os << '\n';
        rss.OutputTo(os, prefix + "rss: ");
        os << '\n';
        text.OutputTo(os, prefix + "text: ");
        os << '\n';
        dataAndStack.OutputTo(os, prefix + "dataAndStack: ");
        os << '\n';
        os << std::defaultfloat;
    }

    template<typename T>
    void MinReduceKnown(MemUsagePageStats MemUsageStats::*pMember,
                    std::shared_ptr<MemUsageStats> myValue,
                    int myRank,
                    int mpiDatatype,
                    bool knowObject) {

        T* target = (myRank == 0) ? &((this->*pMember).min) : nullptr;
        T myCandidate = knowObject ? (myValue.get()->*pMember).min : std::numeric_limits<T>::max();
        
        MPI_Reduce(&myCandidate, target, 1, mpiDatatype, MPI_MIN, 0, MPI_COMM_WORLD);
    }

    template<typename T>
    void MaxReduceKnown(MemUsagePageStats MemUsageStats::*pMember,
                    std::shared_ptr<MemUsageStats> myValue,
                    int myRank,
                    int mpiDatatype,
                    bool knowObject) {
        
        T* target = (myRank == 0) ? &((this->*pMember).max) : nullptr;
        T myCandidate = knowObject ? (myValue.get()->*pMember).max : std::numeric_limits<T>::min();
        MPI_Reduce(&myCandidate, target, 1, mpiDatatype, MPI_MAX, 0, MPI_COMM_WORLD);
    }

    template<typename T>
    void AvgReduceKnown(MemUsagePageStats MemUsageStats::*pMember,
                    std::shared_ptr<MemUsageStats> myValue,
                    int myRank,
                    int mpiDatatype,
                    bool knowObject,
                    uint32_t processCount) {

        T sum;
        T myContribution = knowObject ? (myValue.get()->*pMember).average : 0;
        MPI_Reduce(&myContribution, &sum, 1, mpiDatatype, MPI_SUM, 0, MPI_COMM_WORLD);

        if(myRank == 0) {
            (this->*pMember).average = sum / processCount;
        }
    }

    template<typename T>
    void StdevReduceKnown(MemUsagePageStats MemUsageStats::*pMember,
                            std::shared_ptr<MemUsageStats> myValue,
                            int myRank,
                            int mpiDatatype,
                            bool knowObject,
                            uint32_t processCount) {

        T valSquaredSum;
        T myContribution = knowObject ? ((myValue.get()->*pMember).average * (myValue.get()->*pMember).average) : 0;
        MPI_Reduce(&myContribution, &valSquaredSum, 1, mpiDatatype, MPI_SUM, 0, MPI_COMM_WORLD);

        if(myRank == 0) {
            (this->*pMember).stdev = 
                    std::sqrt((valSquaredSum / processCount)
                        - ((this->*pMember).average * (this->*pMember).average));
        }
    }

    void Aggregate(std::shared_ptr<MemUsageStats> localValue,
                    int myRank,
                    bool localIsValid,
                    uint32_t processCount);
};


} // end namespace xolotlMemUsage

#endif // XMEMUSAGE_MEM_STATS_H
