#include <iomanip>
#include "xolotlMemUsage/profileproc/MemUsageProfile.h"

namespace xolotlMemUsage {

void
MemUsageProfile::outputTo(std::ostream& os) const
{
    using Seconds = std::chrono::duration<double>;

    auto const& bins = profile.GetBins();

    std::time_t startTimeT = AsyncSamplingThreadBase::ClockType::to_time_t(profile.GetStartTimestamp());

    // Output metadata line.  We assume comment char and name have already
    // been output.
    os << "\tnBins: " << bins.size()
        << "\tstartTimestamp: " << std::put_time(std::localtime(&startTimeT), "%F %T")
        << "\tbinWidth_sec: " << Seconds(profile.GetBinWidth()).count() << '\n';

    // Output the CSV header line.
    os << "binNumber,nSamples,vmSize_pages,vmRSS_pages,rss_pages,text_pages,dataAndStack_pages";

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

