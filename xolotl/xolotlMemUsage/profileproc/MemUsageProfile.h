#ifndef XMEMUSAGE_MEM_USAGE_PROFILE_H
#define XMEMUSAGE_MEM_USAGE_PROFILE_H

#include "xolotlMemUsage/memUsageConfig.h"

#if defined(HAVE_STATM)
    #include "xolotlMemUsage/profileproc/Statm/StatmProfiler.h"
    namespace PerProcDataSource = xolotlMemUsage::Statm;
#else
    #error "Configuration error: thought we had a per-proc data source, but no actual data source available."
#endif // defined(HAVE_STATM)


namespace xolotlMemUsage {

struct MemUsageProfile : public IMemUsageSampler::MemUsageData
{
    PerProcDataSource::ProfilerTimeHistogram profile;

    MemUsageProfile(const PerProcDataSource::ProfilerTimeHistogram& _profile)
      : profile(_profile)
    {
        // Nothing else to do.
    }

    void outputTo(std::ostream& os) const;
};

} // namespace xolotlMemUsage

#endif // XMEMUSAGE_MEM_USAGE_PROFILE_H
