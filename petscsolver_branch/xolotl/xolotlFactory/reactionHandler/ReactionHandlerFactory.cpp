#include "ReactionHandlerFactory.h"
#include <PSIClusterReactionNetwork.h>
#include <fstream>
#include <iostream>
#include <mpi.h>

namespace xolotlFactory {

static std::shared_ptr<xolotlCore::IReactionNetwork> theNetworkHandler;

// Create the desired type of network
bool initializeReactionHandler(std::shared_ptr<xolotlCore::HDF5NetworkLoader> networkLoader,
		xolotlCore::Options &options) {
	// Get the current process ID
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	bool ret = true;

	// Check if we want dummy reactions
	auto map = options.getProcesses();
	if (!map["reaction"]) networkLoader->setDummyReactions();
	// Load the network
	theNetworkHandler = networkLoader->load();

	if (procId == 0) {
		std::cout << "\nFactory Message: " << "Master loaded network of size "
				<< theNetworkHandler->size() << "." << std::endl;
	}

	return ret;
}

// Provide access to the network.
std::shared_ptr<xolotlCore::IReactionNetwork> getNetworkHandler() {
	if (!theNetworkHandler) {
		// We have not yet been initialized.
		throw std::string("\nxolotlFactory network handler requested but "
				"it has not been initialized.");
	}
	return theNetworkHandler;
}

} // end namespace xolotlFactory

