#ifndef XMEMUSAGE_IMEMUSAGE_SAMPLER_H
#define XMEMUSAGE_IMEMUSAGE_SAMPLER_H

#include "xolotlCore/IIdentifiable.h"
#include "xolotlMemUsage/MemUsageStats.h"
#include "xolotlMemUsage/MemUsageObjStatistics.h"

namespace xolotlMemUsage {

class IMemUsageSampler : public virtual xolotlCore::IIdentifiable {

public:
    /// Type of the aggregate data we collect
    typedef MemUsageStats ValType;

    /**
     * Type of globally aggregated value statistics.
     * TODO - do we need MemUsageObjStatistics, or is MemUsageStats sufficient?
     */
    typedef MemUsageObjStatistics<MemUsageStats> GlobalStatsType;

    /// Destroy the memory usage sampler.
    virtual ~IMemUsageSampler(void) { }

    /// Start collecting memory usage samples.
    virtual void start(void) = 0;

    /// Stop collecting memory usage samples.
    virtual void stop(void) = 0;

    /// Obtain the sampler's current values.
    virtual ValType getValue(void) const = 0;
};

}//end namespace xolotlMemUsage

#endif // XMEMUSAGE_IMEMUSAGE_SAMPLER_H
