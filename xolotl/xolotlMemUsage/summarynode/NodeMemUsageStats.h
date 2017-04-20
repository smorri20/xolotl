#ifndef XMEMUSAGE_NODE_MEM_STATS_H
#define XMEMUSAGE_NODE_MEM_STATS_H

#include "mpi.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <limits>
#include "xolotlMemUsage/IMemUsageSampler.h"
#include "xolotlMemUsage/summarynode/SysInfoSampler.h"


namespace xolotlMemUsage {

struct NodeMemUsageMetricStats : public SysInfo::MetricStats {

    double stdev;

    /**
     * Construct a NodeMemUsageMetricStats with useful initial values.
     */
    NodeMemUsageMetricStats(void)
      : stdev(std::numeric_limits<double>::quiet_NaN())
    {
        // Nothing else to do.
    }


    /**
     * Construct a NodeMemUsageMetricStats from existing valid data.
     */
    NodeMemUsageMetricStats(uint64_t _min,
                    uint64_t _max,
                    double _avg,
                    double _stdev = std::numeric_limits<double>::quiet_NaN())
      : SysInfo::MetricStats(_min, _max, _avg),
        stdev(_stdev)
    {
        // Nothing else to do.
    }


    /**
     * Construct a NodeMemUsageMetricStats from existing valid data.
     */
    NodeMemUsageMetricStats(const SysInfo::RunningMetricStats& _rps, uint64_t _nSamples)
      : SysInfo::MetricStats(_rps, _nSamples)
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

struct NodeMemUsageStats : public IMemUsageSampler::MemUsageData {
    NodeMemUsageMetricStats totalram;
    NodeMemUsageMetricStats freeram;
    NodeMemUsageMetricStats sharedram;
    NodeMemUsageMetricStats bufferram;
    NodeMemUsageMetricStats totalswap;
    NodeMemUsageMetricStats freeswap;
    NodeMemUsageMetricStats totalhigh;
    NodeMemUsageMetricStats freehigh;
    NodeMemUsageMetricStats mem_unit;

    /**
     * Construct a NodeMemUsageStats with useful initial values.
     */
    NodeMemUsageStats(void) = default;


    /**
     * Construct a NodeMemUsageStats from existing data.
     */
    NodeMemUsageStats(
                const NodeMemUsageMetricStats& _totalram,
                const NodeMemUsageMetricStats& _freeram,
                const NodeMemUsageMetricStats& _sharedram,
                const NodeMemUsageMetricStats& _bufferram,
                const NodeMemUsageMetricStats& _totalswap,
                const NodeMemUsageMetricStats& _freeswap,
                const NodeMemUsageMetricStats& _totalhigh,
                const NodeMemUsageMetricStats& _freehigh,
                const NodeMemUsageMetricStats& _mem_unit
                )
      : totalram(_totalram),
        freeram(_freeram),
        sharedram(_sharedram),
        bufferram(_bufferram),
        totalswap(_totalswap),
        freeswap(_freeswap),
        totalhigh(_totalhigh),
        freehigh(_freehigh),
        mem_unit(_mem_unit)
    {
        // Nothing else to do.
    }



    void OutputTo(std::ostream& os, std::string prefix = "") const {
        // TODO should we just define an ostream inserter for the NodeMemUsageMetricStats?

        os << std::fixed << std::setprecision(2);
        totalram.OutputTo(os, prefix + "totalram: ");
        os << '\n';
        freeram.OutputTo(os, prefix + "freeram: ");
        os << '\n';
        sharedram.OutputTo(os, prefix + "sharedram: ");
        os << '\n';
        bufferram.OutputTo(os, prefix + "bufferram: ");
        os << '\n';
        totalswap.OutputTo(os, prefix + "totalswap: ");
        os << '\n';
        freeswap.OutputTo(os, prefix + "freeswap: ");
        os << '\n';
        totalhigh.OutputTo(os, prefix + "totalhigh: ");
        os << '\n';
        freehigh.OutputTo(os, prefix + "freehigh: ");
        os << '\n';
        mem_unit.OutputTo(os, prefix + "mem_unit: ");
        os << '\n';
        os << std::defaultfloat;
    }

    template<typename T>
    void MinReduceKnown(NodeMemUsageMetricStats NodeMemUsageStats::*pMember,
                    std::shared_ptr<NodeMemUsageStats> myValue,
                    MPI_Comm aggComm,
                    int myRank,
                    int mpiDatatype,
                    bool knowObject) {

        T* target = (myRank == 0) ? &((this->*pMember).min) : nullptr;
        T myCandidate = knowObject ? (myValue.get()->*pMember).min : std::numeric_limits<T>::max();
        
        MPI_Reduce(&myCandidate, target, 1, mpiDatatype, MPI_MIN, 0, aggComm);
    }

    template<typename T>
    void MaxReduceKnown(NodeMemUsageMetricStats NodeMemUsageStats::*pMember,
                    std::shared_ptr<NodeMemUsageStats> myValue,
                    MPI_Comm aggComm,
                    int myRank,
                    int mpiDatatype,
                    bool knowObject) {
        
        T* target = (myRank == 0) ? &((this->*pMember).max) : nullptr;
        T myCandidate = knowObject ? (myValue.get()->*pMember).max : std::numeric_limits<T>::min();
        MPI_Reduce(&myCandidate, target, 1, mpiDatatype, MPI_MAX, 0, aggComm);
    }

    template<typename T>
    void AvgReduceKnown(NodeMemUsageMetricStats NodeMemUsageStats::*pMember,
                    std::shared_ptr<NodeMemUsageStats> myValue,
                    MPI_Comm aggComm,
                    int myRank,
                    int mpiDatatype,
                    bool knowObject,
                    uint32_t processCount) {

        T sum;
        T myContribution = knowObject ? (myValue.get()->*pMember).average : 0;
        MPI_Reduce(&myContribution, &sum, 1, mpiDatatype, MPI_SUM, 0, aggComm);

        if(myRank == 0) {
            (this->*pMember).average = sum / processCount;
        }
    }

    template<typename T>
    void StdevReduceKnown(NodeMemUsageMetricStats NodeMemUsageStats::*pMember,
                            std::shared_ptr<NodeMemUsageStats> myValue,
                            MPI_Comm aggComm,
                            int myRank,
                            int mpiDatatype,
                            bool knowObject,
                            uint32_t processCount) {

        T valSquaredSum;
        T myContribution = knowObject ? ((myValue.get()->*pMember).average * (myValue.get()->*pMember).average) : 0;
        MPI_Reduce(&myContribution, &valSquaredSum, 1, mpiDatatype, MPI_SUM, 0, aggComm);

        if(myRank == 0) {
            (this->*pMember).stdev = 
                    std::sqrt((valSquaredSum / processCount)
                        - ((this->*pMember).average * (this->*pMember).average));
        }
    }

    void Aggregate(std::shared_ptr<NodeMemUsageStats> localValue,
                    MPI_Comm aggComm,
                    int myRank,
                    bool localIsValid,
                    uint32_t processCount);
};


} // end namespace xolotlMemUsage

#endif // XMEMUSAGE_NODE_MEM_STATS_H
