#ifndef XMEMUSAGE_NODE_MEM_STATS_H
#define XMEMUSAGE_NODE_MEM_STATS_H

#include "xolotlMemUsage/memUsageConfig.h"
#include "xolotlMemUsage/summarynode/NodeMemUsageStatsBase.h"


namespace xolotlMemUsage {


struct NodeMemUsageStats : public NodeMemUsageStatsBase {

    NodeMemUsageMetricStats free_count;
    NodeMemUsageMetricStats active_count;
    NodeMemUsageMetricStats inactive_count;
    NodeMemUsageMetricStats wire_count;
    NodeMemUsageMetricStats zero_fill_count;
    NodeMemUsageMetricStats reactivations;
    NodeMemUsageMetricStats pageins;
    NodeMemUsageMetricStats pageouts;
    NodeMemUsageMetricStats faults;
    NodeMemUsageMetricStats cow_faults;
    NodeMemUsageMetricStats lookups;
    NodeMemUsageMetricStats hits;
    NodeMemUsageMetricStats purges;
    NodeMemUsageMetricStats purgeable_count;

    /**
     * Construct a NodeMemUsageStats with useful initial values.
     */
    NodeMemUsageStats(void) = default;


    /**
     * Construct a NodeMemUsageStats from existing data.
     */
    NodeMemUsageStats(
            const NodeMemUsageMetricStats& _free_count,
            const NodeMemUsageMetricStats& _active_count,
            const NodeMemUsageMetricStats& _inactive_count,
            const NodeMemUsageMetricStats& _wire_count,
            const NodeMemUsageMetricStats& _zero_fill_count,
            const NodeMemUsageMetricStats& _reactivations,
            const NodeMemUsageMetricStats& _pageins,
            const NodeMemUsageMetricStats& _pageouts,
            const NodeMemUsageMetricStats& _faults,
            const NodeMemUsageMetricStats& _cow_faults,
            const NodeMemUsageMetricStats& _lookups,
            const NodeMemUsageMetricStats& _hits,
            const NodeMemUsageMetricStats& _purges,
            const NodeMemUsageMetricStats& _purgeable_count
        )
      : free_count(_free_count),
        active_count(_active_count),
        inactive_count(_inactive_count),
        wire_count(_wire_count),
        zero_fill_count(_zero_fill_count),
        reactivations(_reactivations),
        pageins(_pageins),
        pageouts(_pageouts),
        faults(_faults),
        cow_faults(_cow_faults),
        lookups(_lookups),
        hits(_hits),
        purges(_purges),
        purgeable_count(_purgeable_count)
    {
        // Nothing else to do.
    }



    void OutputTo(std::ostream& os, std::string prefix = "") const {
        // TODO should we just define an ostream inserter for the NodeMemUsageMetricStats?

        os << std::fixed << std::setprecision(2);
        free_count.OutputTo(os, prefix + "free_count: ");
        os << '\n';
        active_count.OutputTo(os, prefix + "active_count: ");
        os << '\n';
        inactive_count.OutputTo(os, prefix + "inactive_count: ");
        os << '\n';
        wire_count.OutputTo(os, prefix + "wire_count: ");
        os << '\n';
        zero_fill_count.OutputTo(os, prefix + "zero_fill_count: ");
        os << '\n';
        reactivations.OutputTo(os, prefix + "reactivations: ");
        os << '\n';
        pageins.OutputTo(os, prefix + "pageins: ");
        os << '\n';
        pageouts.OutputTo(os, prefix + "pageouts: ");
        os << '\n';
        faults.OutputTo(os, prefix + "faults: ");
        os << '\n';
        cow_faults.OutputTo(os, prefix + "cow_faults: ");
        os << '\n';
        lookups.OutputTo(os, prefix + "lookups: ");
        os << '\n';
        hits.OutputTo(os, prefix + "hits: ");
        os << '\n';
        purges.OutputTo(os, prefix + "purges: ");
        os << '\n';
        purgeable_count.OutputTo(os, prefix + "purgeable_count: ");
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
