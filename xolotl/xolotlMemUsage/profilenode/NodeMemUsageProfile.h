#ifndef XMEMUSAGE_NODE_MEM_USAGE_PROFILE_NODE_H
#define XMEMUSAGE_NODE_MEM_USAGE_PROFILE_NODE_H

#include "xolotlMemUsage/memUsageConfig.h"

#if defined(HAVE_SYSINFO)
    #include "xolotlMemUsage/profilenode/SysInfo/SysInfoProfiler.h"
    namespace PerNodeDataSource = xolotlMemUsage::SysInfo;
#elif defined(HAVE_MACH_HOST_STATISTICS)
    #include "xolotlMemUsage/profilenode/OSX/OSXProfiler.h"
    namespace PerNodeDataSource = xolotlMemUsage::OSX;
#else
    #error "Configuration error: thought we had a per-node data source, but no actual data source available."
#endif // defined(HAVE_SYSINFO)


namespace xolotlMemUsage {

struct NodeMemUsageProfile : public IMemUsageSampler::MemUsageData
{
    PerNodeDataSource::ProfilerTimeHistogram profile;

    NodeMemUsageProfile(const PerNodeDataSource::ProfilerTimeHistogram& _profile)
      : profile(_profile)
    {
        // Nothing else to do.
    }

    void outputTo(std::ostream& os) const;
};

} // namespace xolotlMemUsage

#endif // XMEMUSAGE_NODE_MEM_USAGE_PROFILE_NODE_H
