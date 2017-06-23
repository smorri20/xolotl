#ifndef XMEMUSAGE_NODE_MEM_STATS_H
#define XMEMUSAGE_NODE_MEM_STATS_H

#include "xolotlMemUsage/memUsageConfig.h"
#include "xolotlMemUsage/summarynode/NodeMemUsageStatsBase.h"


namespace xolotlMemUsage {


struct NodeMemUsageStats : public NodeMemUsageStatsBase {

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

    void Aggregate(std::shared_ptr<NodeMemUsageStats> localValue,
                    MPI_Comm aggComm,
                    int myRank,
                    bool localIsValid,
                    uint32_t processCount);
};


} // end namespace xolotlMemUsage

#endif // XMEMUSAGE_NODE_MEM_STATS_H
