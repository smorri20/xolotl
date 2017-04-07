#ifndef XOLOTLPERF_MEM_STATS_H
#define XOLOTLPERF_MEM_STATS_H

#include "mpi.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <limits>

namespace xolotlPerf {

struct MemPageStats {
    uint64_t min;
    uint64_t max;
    double avg;
    double stdev;

    /**
     * Construct a MemPageStats with useful initial values.
     */
    MemPageStats(void)
      : min(std::numeric_limits<uint64_t>::max()),
        max(std::numeric_limits<uint64_t>::min()),
        avg(std::numeric_limits<double>::quiet_NaN()),
        stdev(std::numeric_limits<double>::quiet_NaN())
    {
        // Nothing else to do.
    }


    /**
     * Construct a MemPageStats from existing valid data.
     */
    MemPageStats(uint64_t _min,
                    uint64_t _max,
                    double _avg,
                    double _stdev = std::numeric_limits<double>::quiet_NaN())
      : min(_min),
        max(_max),
        avg(_avg),
        stdev(_stdev)
    {
        // Nothing else to do.
    }

    void OutputTo(std::ostream& os, std::string prefix = "") const {

        os << prefix 
            << "min: " << min 
            << " max: " << max 
            << " avg: " << avg
            << " stdev: " << stdev;
    }

};

struct MemStats {
    MemPageStats vmSize;
    MemPageStats vmRSS;
    MemPageStats rss;
    MemPageStats text;
    MemPageStats dataAndStack;


    /**
     * Construct a MemStats with useful initial values.
     */
    MemStats(void) = default;


    /**
     * Construct a MemStats from existing data.
     */
    MemStats(const MemPageStats& _vmSize,
                const MemPageStats& _vmRSS,
                const MemPageStats& _rss,
                const MemPageStats& _text,
                const MemPageStats& _dataAndStack)
      : vmSize(_vmSize),
        vmRSS(_vmRSS),
        rss(_rss),
        text(_text),
        dataAndStack(_dataAndStack)
    {
        // Nothing else to do.
    }



    void OutputTo(std::ostream& os, std::string prefix = "") const {
        // TODO should we just define an ostream inserter for the MemPageStats?

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
    void MinReduceKnown(MemPageStats MemStats::*pMember,
                    const MemStats& myValue,
                    int myRank,
                    int mpiDatatype,
                    bool knowObject) {

        T* target = (myRank == 0) ? &((this->*pMember).min) : nullptr;
        T myCandidate = knowObject ? (myValue.*pMember).min : std::numeric_limits<T>::max();
        
        MPI_Reduce(&myCandidate, target, 1, mpiDatatype, MPI_MIN, 0, MPI_COMM_WORLD);
    }

    template<typename T>
    void MaxReduceKnown(MemPageStats MemStats::*pMember,
                    const MemStats& myValue,
                    int myRank,
                    int mpiDatatype,
                    bool knowObject) {
        
        T* target = (myRank == 0) ? &((this->*pMember).max) : nullptr;
        T myCandidate = knowObject ? (myValue.*pMember).max : std::numeric_limits<T>::min();
        MPI_Reduce(&myCandidate, target, 1, mpiDatatype, MPI_MAX, 0, MPI_COMM_WORLD);
    }

    template<typename T>
    void AvgReduceKnown(MemPageStats MemStats::*pMember,
                    const MemStats& myValue,
                    int myRank,
                    int mpiDatatype,
                    bool knowObject,
                    uint32_t processCount) {

        T sum;
        T myContribution = knowObject ? (myValue.*pMember).avg : 0;
        MPI_Reduce(&myContribution, &sum, 1, mpiDatatype, MPI_SUM, 0, MPI_COMM_WORLD);

        if(myRank == 0) {
            (this->*pMember).avg = sum / processCount;
        }
    }

    template<typename T>
    void StdevReduceKnown(MemPageStats MemStats::*pMember,
                            const MemStats& myValue,
                            int myRank,
                            int mpiDatatype,
                            bool knowObject,
                            uint32_t processCount) {

        T valSquaredSum;
        T myContribution = knowObject ? ((myValue.*pMember).avg * (myValue.*pMember).avg) : 0;
        MPI_Reduce(&myContribution, &valSquaredSum, 1, mpiDatatype, MPI_SUM, 0, MPI_COMM_WORLD);

        if(myRank == 0) {
            (this->*pMember).stdev = 
                    std::sqrt((valSquaredSum / processCount)
                        - ((this->*pMember).avg * (this->*pMember).avg));
        }
    }

    void Aggregate(const MemStats& localValue,
                    int myRank,
                    bool localIsValid,
                    uint32_t processCount);
};


} // end namespace xolotlPerf

#endif // XOLOTLPERF_MEM_STATS_H
