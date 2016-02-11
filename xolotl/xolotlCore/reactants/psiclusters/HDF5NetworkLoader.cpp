#include "HDF5NetworkLoader.h"
#include <stdio.h>
#include <limits>
#include <algorithm>
#include <vector>
#include "HeCluster.h"
#include "VCluster.h"
#include "InterstitialCluster.h"
#include "HeVCluster.h"
// #include "HeInterstitialCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <xolotlPerf.h>
#include <HDF5Utils.h>

using namespace xolotlCore;

std::shared_ptr<PSIClusterReactionNetwork> HDF5NetworkLoader::load() {
	// Get the dataset from the HDF5 files
	auto networkVector = xolotlCore::HDF5Utils::readNetwork(fileName);

	// Initialization
	int numHe = 0, numV = 0, numI = 0;
	double formationEnergy = 0.0, migrationEnergy = 0.0;
	double diffusionFactor = 0.0;
	std::vector<std::shared_ptr<Reactant> > reactants;

	// Prepare the network
	std::shared_ptr<PSIClusterReactionNetwork> network = std::make_shared
			< PSIClusterReactionNetwork > (handlerRegistry);

	// Loop on the networkVector
	for (auto lineIt = networkVector.begin(); lineIt != networkVector.end();
			lineIt++) {
		// Composition of the cluster
		numHe = (int) (*lineIt)[0];
		numV = (int) (*lineIt)[1];
		numI = (int) (*lineIt)[2];
		// Create the cluster
		auto nextCluster = createCluster(numHe, numV, numI);

		// Energies
		formationEnergy = (*lineIt)[3];
		migrationEnergy = (*lineIt)[4];
		diffusionFactor = (*lineIt)[5];

		// Set the formation energy
		nextCluster->setFormationEnergy(formationEnergy);
		// Set the diffusion factor and migration energy
		nextCluster->setMigrationEnergy(migrationEnergy);
		nextCluster->setDiffusionFactor(diffusionFactor);
		// Add the cluster to the network
		network->add(nextCluster);
		// Add it to the list so that we can set the network later
		reactants.push_back(nextCluster);
	}

	// Set the reaction network for each reactant
	for (auto reactantsIt = reactants.begin(); reactantsIt != reactants.end();
			++reactantsIt) {
		(*reactantsIt)->setReactionNetwork(network);
	}

	// Clear the reactants vector
	reactants.clear();

	// Apply sectional grouping
	applySectionalGrouping(network);

	return network;
}

void HDF5NetworkLoader::applySectionalGrouping(std::shared_ptr<PSIClusterReactionNetwork> network) {
	// Get the HeV cluster map
	auto heVMap = network->getAll(heVType);

	// Create a temporary vector for the loop
	std::vector<PSICluster *> tempVector;

	// Initialize variables for the loop
	PSICluster * cluster;
	std::shared_ptr<SuperCluster> superCluster;
	static std::map<std::string, int> composition;
	int count = 0, superCount = 0, heIndex = 0, vIndex = vMin;
	double heSize = 0.0, vSize = 0.0, radius = 0.0;

	// Map to know which cluster is in which group
	std::map< std::vector<int>, std::pair<int, int> > clusterGroupMap;
	// Map to know which super cluster gathers which group
	std::map< std::pair<int, int>, SuperCluster *> superGroupMap;

	// Get the number of groups in the helium and vacancy directions
	int nVGroup = (network->getAll(vType).size() - vMin) / vSectionWidth + 1;
	int nHeGroup = (network->getAll(vType).size() * 4) / heSectionWidth + 1;

	// Loop on the vacancy groups
	for (int k = 0; k < nVGroup; k++) {
		// Loop on the helium groups
		for (int j = 0; j < nHeGroup; j++) {
			// Loop within the group
			for (int n = k * vSectionWidth + vMin; n < (k + 1) * vSectionWidth + vMin; n++) {
				for (int m = j * heSectionWidth + 1; m < (j + 1) * heSectionWidth + 1; m++) {
					// Get the corresponding cluster
					std::vector<int> compositionVector = { m, n, 0 };
					// Get the product of the same type as the second reactant
					cluster = (PSICluster *) network->getCompound(heVType, compositionVector);

					// Verify if the cluster exists
					if (!cluster) continue;

					// Increment the counter
					count++;

					// Add this cluster to the temporary vector
					tempVector.push_back(cluster);
					heSize += (double) m;
					vSize += (double) n;
					radius += cluster->getReactionRadius();
					// Keep the information of the group
					clusterGroupMap[compositionVector] = std::make_pair(j, k);
				}
			}

			// Check if there were clusters in this group
			if (count == 0) continue;

			// Average all values
			heSize = heSize / (double) count;
			vSize = vSize / (double) count;
			radius = radius / (double) count;
			// Create the cluster
			superCluster = std::make_shared<SuperCluster>(heSize, vSize,
					count, heSectionWidth, vSectionWidth, radius, handlerRegistry);
			// Set the HeV vector
			superCluster->setHeVVector(tempVector);
			// Add this cluster to the network and clusters
			network->addSuper(superCluster);
			clusters.push_back(superCluster);
			// Keep the information of the group
			superGroupMap[std::make_pair(j, k)] = superCluster.get();

			std::cout << j << " " << k << " " << count << " " << superCluster->getName() << std::endl;

			// Reinitialize everything
			heSize = 0.0, vSize = 0.0, radius = 0.0;
			count = 0;
			tempVector.clear();
		}
	}

	// Initialize variables for the loop
	SuperCluster * newCluster;
	// Loop on all the reactants to update the pairs vector with super clusters
	auto reactants = network->getAll();
	// Get the super cluster map
	auto superMap = network->getAll(superType);
	for (int i = 0; i < reactants->size(); i++) {
		// Get the cluster
		cluster = (PSICluster *) reactants->at(i);
		// Get their production and dissociation vectors
		auto react = cluster->reactingPairs;
		auto combi = cluster->combiningReactants;
		auto disso = cluster->dissociatingPairs;
		auto emi = cluster->emissionPairs;

		// Loop on its reacting pairs
		for (int l = 0; l < react.size(); l++) {
			// Test the first reactant
			if (react[l].first->getType() == heVType) {
				// Get its composition
				composition = react[l].first->getComposition();
				// Test its size
				if (composition[vType] >= vMin) {
					// It has to be replaced by a super cluster
					std::vector<int> compositionVector = { composition[heType], composition[vType], 0 };
					newCluster = superGroupMap[clusterGroupMap[compositionVector]];
					react[l].first = newCluster;
					react[l].firstHeDistance = newCluster->getHeDistance(composition[heType]);
					react[l].firstVDistance = newCluster->getVDistance(composition[vType]);
				}
			}

			// Test the second reactant
			if (react[l].second->getType() == heVType) {
				// Get its composition
				composition = react[l].second->getComposition();
				// Test its size
				if (composition[vType] >= vMin) {
					// It has to be replaced by a super cluster
					std::vector<int> compositionVector = { composition[heType], composition[vType], 0 };
					newCluster = superGroupMap[clusterGroupMap[compositionVector]];
					react[l].second = newCluster;
					react[l].secondHeDistance = newCluster->getHeDistance(composition[heType]);
					react[l].secondVDistance = newCluster->getVDistance(composition[vType]);
				}
			}
		}

		// Loop on its combining reactants
		for (int l = 0; l < combi.size(); l++) {
			// Test the combining reactant
			if (combi[l].combining->getType() == heVType) {
				// Get its composition
				composition = combi[l].combining->getComposition();
				// Test its size
				if (composition[vType] >= vMin) {
					// It has to be replaced by a super cluster
					std::vector<int> compositionVector = { composition[heType], composition[vType], 0 };
					newCluster = superGroupMap[clusterGroupMap[compositionVector]];
					combi[l].combining = newCluster;
					combi[l].heDistance = newCluster->getHeDistance(composition[heType]);
					combi[l].vDistance = newCluster->getVDistance(composition[vType]);
				}
			}
		}

		// Loop on its dissociating pairs
		for (int l = 0; l < disso.size(); l++) {
			// Test the first reactant
			if (disso[l].first->getType() == heVType) {
				// Get its composition
				composition = disso[l].first->getComposition();
				// Test its size
				if (composition[vType] >= vMin) {
					// It has to be replaced by a super cluster
					std::vector<int> compositionVector = { composition[heType], composition[vType], 0 };
					newCluster = superGroupMap[clusterGroupMap[compositionVector]];
					disso[l].first = newCluster;
					disso[l].firstHeDistance = newCluster->getHeDistance(composition[heType]);
					disso[l].firstVDistance = newCluster->getVDistance(composition[vType]);
				}
			}

			// Test the second reactant
			if (disso[l].second->getType() == heVType) {
				// Get its composition
				composition = disso[l].second->getComposition();
				// Test its size
				if (composition[vType] >= vMin) {
					// It has to be replaced by a super cluster
					std::vector<int> compositionVector = { composition[heType], composition[vType], 0 };
					newCluster = superGroupMap[clusterGroupMap[compositionVector]];
					disso[l].second = newCluster;
					disso[l].secondHeDistance = newCluster->getHeDistance(composition[heType]);
					disso[l].secondVDistance = newCluster->getVDistance(composition[vType]);
				}
			}
		}

		// Loop on its emission pairs
		for (int l = 0; l < emi.size(); l++) {
			// Test the first reactant
			if (emi[l].first->getType() == heVType) {
				// Get its composition
				composition = emi[l].first->getComposition();
				// Test its size
				if (composition[vType] >= vMin) {
					// It has to be replaced by a super cluster
					std::vector<int> compositionVector = { composition[heType], composition[vType], 0 };
					newCluster = superGroupMap[clusterGroupMap[compositionVector]];
					emi[l].first = newCluster;
					emi[l].firstHeDistance = newCluster->getHeDistance(composition[heType]);
					emi[l].firstVDistance = newCluster->getVDistance(composition[vType]);
				}
			}

			// Test the second reactant
			if (emi[l].second->getType() == heVType) {
				// Get its composition
				composition = emi[l].second->getComposition();
				// Test its size
				if (composition[vType] >= vMin) {
					// It has to be replaced by a super cluster
					std::vector<int> compositionVector = { composition[heType], composition[vType], 0 };
					newCluster = superGroupMap[clusterGroupMap[compositionVector]];
					emi[l].second = newCluster;
					emi[l].secondHeDistance = newCluster->getHeDistance(composition[heType]);
					emi[l].secondVDistance = newCluster->getVDistance(composition[vType]);
				}
			}
		}

		// Set their production and dissociation vectors
		cluster->reactingPairs = react;
		cluster->combiningReactants = combi;
		cluster->dissociatingPairs = disso;
		cluster->emissionPairs = emi;
	}

	// Set the reaction network for each super reactant
	for (auto reactantsIt = superMap.begin(); reactantsIt != superMap.end();
			++reactantsIt) {
		(*reactantsIt)->setReactionNetwork(network);
	}

	// Remove HeV clusters bigger than vMin from the network
	// Loop on the HeV clusters
	for (int i = 0; i < heVMap.size(); i++) {
		// Get the cluster and its composition
		cluster = (PSICluster *) heVMap.at(i);
		composition = cluster->getComposition();

		// Skip the clusters that are too small
		if (composition[vType] < vMin) continue;

		// Remove the reactant from the network
		network->removeReactant(heVMap.at(i));
	}

	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeNetwork();

	return;
}
