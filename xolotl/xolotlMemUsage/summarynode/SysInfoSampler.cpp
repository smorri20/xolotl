#include "xolotlMemUsage/summarynode/SysInfoSampler.h"
#include "xolotlMemUsage/summarynode/NodeMemUsageStats.h"

namespace xolotlMemUsage {

namespace SysInfo {

std::shared_ptr<IMemUsageSampler::MemUsageData>
Sampler::GetCurrentStats(void) const
{
    auto const& runningData = GetRunningSampleData();
    auto nSamples = runningData.nSamples;
    return std::make_shared<NodeMemUsageStats>(
                NodeMemUsageMetricStats(runningData.totalram, nSamples),
                NodeMemUsageMetricStats(runningData.freeram, nSamples),
                NodeMemUsageMetricStats(runningData.sharedram, nSamples),
                NodeMemUsageMetricStats(runningData.bufferram, nSamples),
                NodeMemUsageMetricStats(runningData.totalswap, nSamples),
                NodeMemUsageMetricStats(runningData.freeswap, nSamples),
                NodeMemUsageMetricStats(runningData.totalhigh, nSamples),
                NodeMemUsageMetricStats(runningData.freehigh, nSamples),
                NodeMemUsageMetricStats(runningData.mem_unit, nSamples));
}

} // end namespace SysInfo

} // end namespace xolotlMemUsage

