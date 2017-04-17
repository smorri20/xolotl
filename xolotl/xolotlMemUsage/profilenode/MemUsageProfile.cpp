#include <iomanip>
#include "xolotlMemUsage/profilenode/NodeMemUsageProfile.h"

namespace xolotlMemUsage {

void
NodeMemUsageProfile::outputTo(std::ostream& os) const
{
    using Seconds = std::chrono::duration<double>;

    auto const& bins = profile.GetBins();

    std::time_t startTimeT = AsyncSamplingThreadBase::ClockType::to_time_t(profile.GetStartTimestamp());

    os << "\tnBins: " << bins.size()
        << "\tstartTimestamp: " << std::put_time(std::localtime(&startTimeT), "%F %T")
        << "\tbinWidth_sec: " << Seconds(profile.GetBinWidth()).count() << '\n'
        << "# bin\tnSamples\ttotalram\tfreeram\tsharedram\tbufferram\ttotalswap\tfreeswap\ttotahigh\tfreehigh\n";

    auto binIdx = 0;
    for(auto currBin : bins)
    {
        auto const& currBinMetric = currBin.GetMetricValue();

        os << binIdx << '\t'
            << currBin.GetNumSamples() << '\t'
            << currBin.GetMetricValue()
            << std::endl;

        ++binIdx;
    }
}

} // namespace xolotlMemUsage

