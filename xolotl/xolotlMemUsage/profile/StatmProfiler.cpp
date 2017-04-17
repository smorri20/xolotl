#include "xolotlMemUsage/profile/StatmProfiler.h"
#include "xolotlMemUsage/profile/MemUsageProfile.h"

namespace xolotlMemUsage {

namespace Statm {

const unsigned int Profiler::nTimeHistogramBins = 40;

const AsyncSamplingThreadBase::ClockType::duration Profiler::initialBinWidth = std::chrono::duration<uint64_t, std::milli>(1000);

std::shared_ptr<IMemUsageSampler::MemUsageData>
Profiler::GetCurrentProfile(void) const
{
    return std::make_shared<MemUsageProfile>(GetRunningSampleData().hist);
}

} // end namespace Statm

} // namespace xolotlMemUsage

