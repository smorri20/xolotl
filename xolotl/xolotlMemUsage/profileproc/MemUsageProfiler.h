#ifndef XMEMUSAGE_MEM_USAGE_PROFILER_H
#define XMEMUSAGE_MEM_USAGE_PROFILER_H

#include "xolotlMemUsage/memUsageConfig.h"
#include "xolotlMemUsage/common/MemUsageSamplerBase.h"

#if defined(HAVE_STATM)
    #include "xolotlMemUsage/profileproc/Statm/StatmProfiler.h"
    namespace PerProcDataSource = xolotlMemUsage::Statm;
#else
    #error "Configuration error: thought we had a per-proc data source, but no actual data source available."
#endif // defined(HAVE_STATM)



namespace xolotlMemUsage {

class MemUsageProfiler : public MemUsageSamplerBase<PerProcDataSource::Profiler> {

public:

    /// Disallow building an object without a name.
    MemUsageProfiler(void) = delete;


    /**
     * Construct a MemUsageProfiler with given name.
     *
     * @param name The object's name.
     */
    MemUsageProfiler(const std::string& name)
      : MemUsageSamplerBase<PerProcDataSource::Profiler>(name) {

        // Nothing else to do.
    }

    /**
     * Obtain the data we've collected.
     */
    virtual std::shared_ptr<IMemUsageSampler::MemUsageData> getValue(void) const {
        return GetCurrentProfile();
    }
};

} // namespace xolotlMemUsage

#endif // XMEMUSAGE_MEM_USAGE_PROFILER_H
