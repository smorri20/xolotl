#ifndef XMEMUSAGE_MEM_STATS_H
#define XMEMUSAGE_MEM_STATS_H

#include "xolotlMemUsage/memUsageConfig.h"
#include "xolotlMemUsage/summaryproc/MemUsageStatsBase.h"


namespace xolotlMemUsage {


struct MemUsageStats : public MemUsageStatsBase {

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

    void Aggregate(std::shared_ptr<MemUsageStats> localValue,
                    MPI_Comm aggComm,
                    int myRank,
                    bool localIsValid,
                    uint32_t processCount);
};


} // end namespace xolotlMemUsage

#endif // XMEMUSAGE_MEM_STATS_H
