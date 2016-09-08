#include "NEClusterNetworkLoader.h"
#include <NEClusterReactionNetwork.h>
#include <XeCluster.h>
#include <HDF5Utils.h>
#include <xolotlPerf.h>

using namespace xolotlCore;

///**
// * This operation converts a string to a double, taking in to account the fact
// * that the input file may contain keys such as "infinite."
// *
// * @param inString the string to be converted
// * @return the string as a double
// */
//static inline double convertStrToDouble(const std::string& inString) {
//	return (inString.compare("infinite") == 0) ?
//			std::numeric_limits<double>::infinity() :
//			strtod(inString.c_str(), NULL);
//}

std::shared_ptr<NECluster> NEClusterNetworkLoader::createNECluster(int numXe,
		int numV, int numI) {
	// Local Declarations
	std::shared_ptr<NECluster> cluster;

	// Determine the type of the cluster given the number of each species.
	// Create a new cluster by that type and specify the names of the
	// property keys.
	if (numXe > 0) {
		// Create a new XeVCluster
		cluster = std::make_shared<XeCluster>(numXe, handlerRegistry);
	}

	return cluster;
}

NEClusterNetworkLoader::NEClusterNetworkLoader(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	networkStream = nullptr;
	handlerRegistry = registry;
	fileName = "";
	dummyReactions = false;

	return;
}

NEClusterNetworkLoader::NEClusterNetworkLoader(const std::shared_ptr<std::istream> stream,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
	networkStream = stream;
	handlerRegistry = registry;
	fileName = "";
	dummyReactions = false;

	return;
}

std::shared_ptr<IReactionNetwork> NEClusterNetworkLoader::load() {
	// Get the dataset from the HDF5 files
	auto networkVector = xolotlCore::HDF5Utils::readNetwork(fileName);

	// Initialization
	int numXe = 0, numV = 0, numI = 0;
	double formationEnergy = 0.0, migrationEnergy = 0.0;
	double diffusionFactor = 0.0;
	std::vector<std::shared_ptr<Reactant> > reactants;

	// Prepare the network
	std::shared_ptr<NEClusterReactionNetwork> network = std::make_shared
			< NEClusterReactionNetwork > (handlerRegistry);

	// Loop on the networkVector
	for (auto lineIt = networkVector.begin(); lineIt != networkVector.end();
			lineIt++) {

		// Composition of the cluster
		numXe = (int) (*lineIt)[0];
		numV = (int) (*lineIt)[1];
		numI = (int) (*lineIt)[2];
		// Create the cluster
		auto nextCluster = createNECluster(numXe, numV, numI);

		// Energies
		formationEnergy = (*lineIt)[3];
		migrationEnergy = (*lineIt)[4];
		diffusionFactor = (*lineIt)[5];

		// Set the formation energy
		nextCluster->setFormationEnergy(formationEnergy);
		// Set the diffusion factor and migration energy
		nextCluster->setMigrationEnergy(migrationEnergy);
		nextCluster->setDiffusionFactor(diffusionFactor);

		// Check if we want dummy reactions
		if (dummyReactions) {
			// Create a dummy cluster (Reactant) from the existing cluster
			auto dummyCluster = std::static_pointer_cast<Reactant> (nextCluster->Reactant::clone());
			// Add the cluster to the network
			network->add(dummyCluster);
			// Add it to the list so that we can set the network later
			reactants.push_back(dummyCluster);
		}
		else {
			// Add the cluster to the network
			network->add(nextCluster);
			// Add it to the list so that we can set the network later
			reactants.push_back(nextCluster);
		}
	}

	// Set the reaction network for each reactant
	for (auto reactantsIt = reactants.begin(); reactantsIt != reactants.end();
			++reactantsIt) {
		(*reactantsIt)->setReactionNetwork(network);
	}

	return network;
}
