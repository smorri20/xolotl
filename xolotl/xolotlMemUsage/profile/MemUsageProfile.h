#ifndef XMEMUSAGE_MEM_USAGE_PROFILE_H
#define XMEMUSAGE_MEM_USAGE_PROFILE_H

#include "xolotlMemUsage/profile/StatmProfiler.h"

namespace xolotlMemUsage {

struct MemUsageProfile : public IMemUsageSampler::MemUsageData
{
    Statm::ProfilerTimeHistogram profile;

    MemUsageProfile(const Statm::ProfilerTimeHistogram& _profile)
      : profile(_profile)
    {
        // Nothing else to do.
    }

    void outputTo(std::ostream& os) const;
};

} // namespace xolotlMemUsage

#endif // XMEMUSAGE_MEM_USAGE_PROFILE_H
