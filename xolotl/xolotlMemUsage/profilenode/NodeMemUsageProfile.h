#ifndef XMEMUSAGE_NODE_MEM_USAGE_PROFILE_NODE_H
#define XMEMUSAGE_NODE_MEM_USAGE_PROFILE_NODE_H

#include "xolotlMemUsage/profilenode/SysInfoProfiler.h"

namespace xolotlMemUsage {

struct NodeMemUsageProfile : public IMemUsageSampler::MemUsageData
{
    SysInfo::ProfilerTimeHistogram profile;

    NodeMemUsageProfile(const SysInfo::ProfilerTimeHistogram& _profile)
      : profile(_profile)
    {
        // Nothing else to do.
    }

    void outputTo(std::ostream& os) const;
};

} // namespace xolotlMemUsage

#endif // XMEMUSAGE_NODE_MEM_USAGE_PROFILE_NODE_H
