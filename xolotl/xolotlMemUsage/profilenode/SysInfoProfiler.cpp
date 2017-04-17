#include "xolotlMemUsage/profilenode/SysInfoProfiler.h"
#include "xolotlMemUsage/profilenode/NodeMemUsageProfile.h"

namespace xolotlMemUsage {

namespace SysInfo {

const unsigned int Profiler::nTimeHistogramBins = 40;

const AsyncSamplingThreadBase::ClockType::duration Profiler::initialBinWidth = std::chrono::duration<uint64_t, std::milli>(1000);

std::shared_ptr<IMemUsageSampler::MemUsageData>
Profiler::GetCurrentProfile(void) const
{
    return std::make_shared<NodeMemUsageProfile>(GetRunningSampleData().hist);
}

} // end namespace SysInfo

} // namespace xolotlMemUsage

