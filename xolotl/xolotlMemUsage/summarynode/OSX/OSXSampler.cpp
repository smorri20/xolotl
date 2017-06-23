#include "xolotlMemUsage/summarynode/OSX/OSXSampler.h"
#include "xolotlMemUsage/summarynode/OSX/NodeMemUsageStats.h"

namespace xolotlMemUsage {

namespace OSX {

std::shared_ptr<IMemUsageSampler::MemUsageData>
Sampler::GetCurrentStats(void) const
{
    auto const& runningData = GetRunningSampleData();
    auto nSamples = runningData.nSamples;
    return std::make_shared<NodeMemUsageStats>(
            NodeMemUsageMetricStats(runningData.free_count, nSamples),
            NodeMemUsageMetricStats(runningData.active_count, nSamples),
            NodeMemUsageMetricStats(runningData.inactive_count, nSamples),
            NodeMemUsageMetricStats(runningData.wire_count, nSamples),
            NodeMemUsageMetricStats(runningData.zero_fill_count, nSamples),
            NodeMemUsageMetricStats(runningData.reactivations, nSamples),
            NodeMemUsageMetricStats(runningData.pageins, nSamples),
            NodeMemUsageMetricStats(runningData.pageouts, nSamples),
            NodeMemUsageMetricStats(runningData.faults, nSamples),
            NodeMemUsageMetricStats(runningData.cow_faults, nSamples),
            NodeMemUsageMetricStats(runningData.lookups, nSamples),
            NodeMemUsageMetricStats(runningData.hits, nSamples),
            NodeMemUsageMetricStats(runningData.purges, nSamples),
            NodeMemUsageMetricStats(runningData.purgeable_count, nSamples)
        );
}

} // end namespace OSX

} // end namespace xolotlMemUsage

