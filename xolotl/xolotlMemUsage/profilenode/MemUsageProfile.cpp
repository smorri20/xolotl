#include <iomanip>
#include "xolotlMemUsage/profilenode/NodeMemUsageProfile.h"

namespace xolotlMemUsage {

void
NodeMemUsageProfile::outputTo(std::ostream& os) const
{
    using Seconds = std::chrono::duration<double>;

    auto const& bins = profile.GetBins();

    std::time_t startTimeT = AsyncSamplingThreadBase::ClockType::to_time_t(profile.GetStartTimestamp());

    // Output rest of metadata line - we assume someone has
    // started the line with a '#' and the name of the histogram.
    os << "\tnBins: " << bins.size()
        << "\tstartTimestamp: " << std::put_time(std::localtime(&startTimeT), "%F %T")
        << "\tbinWidth_sec: " << Seconds(profile.GetBinWidth()).count() << '\n';


    // Output the CSV header line.
    os << "binNumber,nSamples,totalram_bytes,freeram_bytes,sharedram_bytes,bufferram_bytes,totalswap_bytes,freeswap_bytes,totalhigh_bytes,freehigh_bytes\n";

    auto binIdx = 0;
    for(auto currBin : bins)
    {
        auto const& currBinMetric = currBin.GetMetricValue();

        os << binIdx << ','
            << currBin.GetNumSamples() << ','
            << currBin.GetMetricValue()
            << std::endl;

        ++binIdx;
    }
}

} // namespace xolotlMemUsage

