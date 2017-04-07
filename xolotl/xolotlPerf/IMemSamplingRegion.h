#ifndef IMEM_SAMPLING_REGION_H
#define IMEM_SAMPLING_REGION_H

#include "xolotlCore/IIdentifiable.h"
#include "xolotlPerf/MemStats.h"
#include "xolotlPerf/PerfObjStatistics.h"

namespace xolotlPerf {

class IMemSamplingRegion : public virtual xolotlCore::IIdentifiable {

public:
    /// Type of the aggregate data we collect
    typedef MemStats ValType;

    /**
     * Type of globally aggregated value statistics.
     * TODO - do we need PerfObjStatistics, or is MemStats sufficient?
     */
    typedef PerfObjStatistics<MemStats> GlobalStatsType;

    /// Destroy the memory usage sampling region.
    virtual ~IMemSamplingRegion(void) { }

    /// Start collecting memory usage samples for the given region.
    virtual void start(void) = 0;

    /// Stop collecting memory usage samples for the given region.
    virtual void stop(void) = 0;

    /// Obtain the sampling region's current value.
    virtual ValType getValue(void) const = 0;
};

}//end namespace xolotlPerf

#endif // IMEM_SAMPLING_REGION_H
