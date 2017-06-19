#include "mpi.h"
#include <iostream>
#include <sstream>
#include "xolotlMemUsage/memUsageConfig.h"
#include "xolotlMemUsage/xolotlMemUsage.h"
#include "xolotlMemUsage/dummy/DummyHandlerRegistry.h"

#if defined(HAVE_PER_NODE_DATA_SOURCE)
#include "xolotlMemUsage/summarynode/SummaryNodeHandlerRegistry.h"
#include "xolotlMemUsage/profilenode/ProfileNodeHandlerRegistry.h"
#endif // defined(HAVE_PER_NODE_DATA_SOURCE)

#if defined(HAVE_PER_PROC_DATA_SOURCE)
#include "xolotlMemUsage/summaryproc/SummaryProcHandlerRegistry.h"
#include "xolotlMemUsage/profileproc/ProfileProcHandlerRegistry.h"
#endif // defined(HAVE_PER_PROC_DATA_SOURCE)


namespace xolotlMemUsage {

static std::shared_ptr<IHandlerRegistry> theHandlerRegistry;

// Create the desired type of handler registry.
void initialize(IHandlerRegistry::RegistryType rtype,
                IHandlerRegistry::SamplingInterval samplingInterval,
                std::string profileFilename) {

    // Set our sampling interval to the one given.
    AsyncSamplingThreadBase::SetSamplingInterval(samplingInterval);

    // If user wants node-level information, figure out if we are 
    // going to be participating in the memory usage data collection.
    // TODO this requires MPI-3 functionality.  What to do if that
    // isn't available?
    // TODO how to handle situations where the "blessed" process
    // doesn't define some memory regions that others do?
    bool participating = true;
    MPI_Comm aggComm;
    int aggCommRank = -1;
    if((rtype == IHandlerRegistry::summaryNode) or
            (rtype == IHandlerRegistry::profileNode)) {

        MPI_Comm nodeLocalComm;
        MPI_Comm_split_type(MPI_COMM_WORLD,
                            MPI_COMM_TYPE_SHARED,
                            0,  // key
                            MPI_INFO_NULL,
                            &nodeLocalComm);

        int myNodeLocalRank;
        MPI_Comm_rank(nodeLocalComm, &myNodeLocalRank);
        participating = (myNodeLocalRank == 0);

        // We need a communicator that only the 'blessed' processes
        // belong to, since only those will have data to be aggregated.
        // We pass in a color and a key when creating the new communicator.
        // All communicators with same color are in the same communicator.
        // Any that pass in MPI_UNDEFINED for the color get put into 
        // a communicator of one.
        // The key determines the ordering of the ranks within the
        // new communciator.  We use our rank within MPI_COMM_WORLD,
        // so that (hopefully) rank 0 in MPI_COMM_WORLD will always 
        // be rank 0 in the new communicator.
        int myCWRank = -1;
        MPI_Comm_rank(MPI_COMM_WORLD, &myCWRank);
        MPI_Comm_split(MPI_COMM_WORLD,
                        participating ? 0 : MPI_UNDEFINED,  // color - procs with same color are in same output communicator
                        myCWRank,               // key for ordering ranks in new communicator - should be in same
                        &aggComm);
        if(participating) {
            // Figure out our rank within the aggregating communicator.
            MPI_Comm_rank(aggComm, &aggCommRank);
        }
        else {
            // We won't be providing data, so we can release the
            // communicator we were given.
            // Should be a purely local operation since we 
            // used a color of MPI_UNDEFINED.
            if(aggComm != MPI_COMM_NULL) {
                MPI_Comm_free(&aggComm);
            }
        }
    }
        
    // Build a registry that will produce the desired type of memory
    // usage samplers.
    switch (rtype) {
    case IHandlerRegistry::dummy:
        theHandlerRegistry = std::make_shared<DummyHandlerRegistry>();
        break;

#if defined(HAVE_PER_PROC_DATA_SOURCE)
    case IHandlerRegistry::summaryProc:
        theHandlerRegistry = std::make_shared<SummaryProcHandlerRegistry>();
        break;

    case IHandlerRegistry::profileProc:
        theHandlerRegistry = std::make_shared<ProfileProcHandlerRegistry>(profileFilename);
        break;
#endif // defined(HAVE_PER_PROC_DATA_SOURCE)

#if defined(HAVE_PER_NODE_DATA_SOURCE)
    case IHandlerRegistry::summaryNode:
        if(participating) {
            // Build the registry to use the aggregating communicator.
            assert(aggCommRank != -1);
            theHandlerRegistry = std::make_shared<SummaryNodeHandlerRegistry>(aggComm, aggCommRank);
        }
        else {
            // We won't be providing any data.
            theHandlerRegistry = std::make_shared<DummyHandlerRegistry>();
        }
        break;

    case IHandlerRegistry::profileNode:
        if(participating) {
            // Build the registry to use the aggregating communicator.
            assert(aggCommRank != -1);
            theHandlerRegistry = std::make_shared<ProfileNodeHandlerRegistry>(profileFilename, aggComm, aggCommRank);
        }
        else {
            // We won't be providing any data.
            theHandlerRegistry = std::make_shared<DummyHandlerRegistry>();
        }
        break;
#endif // defined(HAVE_PER_NODE_DATA_SOURCE)

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

