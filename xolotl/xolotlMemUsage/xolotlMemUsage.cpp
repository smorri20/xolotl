#include <iostream>
#include <sstream>
#include "xolotlMemUsage/memUsageConfig.h"
#include "xolotlMemUsage/xolotlMemUsage.h"
#include "xolotlMemUsage/dummy/DummyHandlerRegistry.h"
#include "xolotlMemUsage/standard/StdHandlerRegistry.h"
#include "xolotlMemUsage/profileproc/ProfileProcHandlerRegistry.h"
#include "xolotlMemUsage/profilenode/ProfileNodeHandlerRegistry.h"


namespace xolotlMemUsage {

static std::shared_ptr<IHandlerRegistry> theHandlerRegistry;

// Create the desired type of handler registry.
void initialize(IHandlerRegistry::RegistryType rtype,
                IHandlerRegistry::SamplingInterval samplingInterval,
                std::string profileFilename) {

    // Set our sampling interval to the one given.
    AsyncSamplingThreadBase::SetSamplingInterval(samplingInterval);

	switch (rtype) {
	case IHandlerRegistry::dummy:
		theHandlerRegistry = std::make_shared<DummyHandlerRegistry>();
		break;

	case IHandlerRegistry::std:
        theHandlerRegistry = std::make_shared<StdHandlerRegistry>();
		break;

    case IHandlerRegistry::profileproc:
        theHandlerRegistry = std::make_shared<ProfileProcHandlerRegistry>(profileFilename);
        break;

    case IHandlerRegistry::profilenode:
        // The user wants per-node profiling.  We only need one process
        // on the node to collect the profile data.
        // We select the one process that will collect profile data
        // by creating a node-local communicator, and "bless" whatever
        // process is rank 0 in that communicator.
        // Everyone else will build a dummy registry.
        // TODO this requires MPI-3 functionality - what to do if it
        // isn't available?
        // TODO how to handle situations where the "blessed" process 
        // doesn't define some memory profiling regions that others do?
        {
            MPI_Comm nodeLocalComm;
            MPI_Comm_split_type(MPI_COMM_WORLD,
                                MPI_COMM_TYPE_SHARED,
                                0,  // key
                                MPI_INFO_NULL,
                                &nodeLocalComm);

            int myNodeLocalRank;
            MPI_Comm_rank(nodeLocalComm, &myNodeLocalRank);

            if(myNodeLocalRank == 0)
            {
                theHandlerRegistry = std::make_shared<ProfileNodeHandlerRegistry>(profileFilename);
            }
            else
            {
                theHandlerRegistry = std::make_shared<DummyHandlerRegistry>();
            }
        }
        break;

	default:
		throw std::invalid_argument(
				"unrecognized memory usage handler registry type requested");
		break;
	}
}

// Provide access to our handler registry.
std::shared_ptr<IHandlerRegistry> getHandlerRegistry(void) {
	if (!theHandlerRegistry) {
		throw std::runtime_error(
				"Request for xolotlMemUsage handler registry before xolotlMemUsage library has been initialized");
	}
	return theHandlerRegistry;
}

} // end namespace xolotlMemUsage

