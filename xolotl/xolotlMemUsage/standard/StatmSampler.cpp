#include <fstream>
#include "xolotlMemUsage/standard/StatmSampler.h"
#include "xolotlMemUsage/standard/MemUsageStats.h"

namespace xolotlMemUsage {

namespace Statm {

std::shared_ptr<IMemUsageSampler::MemUsageData>
Sampler::GetCurrentStats(void) const
{
    auto const& runningData = GetRunningSampleData();
    auto nSamples = runningData.nSamples;
    return std::make_shared<MemUsageStats>(MemUsagePageStats(runningData.vmSize, nSamples),
                                MemUsagePageStats(runningData.vmRSS, nSamples),
                                MemUsagePageStats(runningData.rss, nSamples),
                                MemUsagePageStats(runningData.text, nSamples),
                                MemUsagePageStats(runningData.dataAndStack, nSamples));
}


} // end namespace Statm

} // end namespace xolotlMemUsage

