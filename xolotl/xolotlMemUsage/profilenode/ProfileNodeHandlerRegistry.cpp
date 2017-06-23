#include "mpi.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <tuple>
#include <cmath>
#include <cstring>
#include "xolotlMemUsage/profilenode/ProfileNodeHandlerRegistry.h"
#include "xolotlMemUsage/profilenode/NodeMemUsageProfiler.h"

// Pull in definitions of TimeHistogram template methods so we can
// instantiate the methods we need.
#include "xolotlMemUsage/common/TimeHistogram.cpp"

namespace xolotlMemUsage {

std::shared_ptr<IHandlerRegistry::GlobalData>
ProfileNodeHandlerRegistry::collectData(void) const {

    auto ret = std::make_shared<GlobalMemUsageData>();

    // Determine my MPI rank.
	int myRank;
    int cwSize;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &cwSize);

#if READY
    // TODO implement aggregation of memory usage profiles.
    // There are lots of questions to answer to do this.
    // Need handle:
    // * bin widths that might not be the same across all processes
    // * start times that almost certainly don't match
    // * metrics that make sense to aggregate.
#else
    // Output per-process memory profile data.
    // Replace '%r' in output file name with our MPI rank.
    std::string actualFilename = oFilename;
    auto pos = actualFilename.find("%r");
    if(pos != std::string::npos) {
        // How many digits wide is our maximum rank?
        uint32_t maxWidth = std::rint(std::ceil(std::log10(cwSize - 1)));

        // Build our output file name.
        std::ostringstream rstr;
        rstr << std::setw(maxWidth) << std::setfill('0') << myRank;
        actualFilename.replace(pos, 2, rstr.str());
    }

    // Output the profile(s) we collected.
    std::ofstream pstr(actualFilename);
    for(auto const& currMapItem : allSamplers) {

        auto currProfiler = std::dynamic_pointer_cast<NodeMemUsageProfiler>(currMapItem.second);
        auto currProfile = std::dynamic_pointer_cast<NodeMemUsageProfile>(currProfiler->GetCurrentProfile());

        pstr << "# name: " << currMapItem.first;
        currProfile->outputTo(pstr);
        pstr << "\n\n";
    }
#endif // READY

    return ret;
}

} // namespace xolotlMemUsage

