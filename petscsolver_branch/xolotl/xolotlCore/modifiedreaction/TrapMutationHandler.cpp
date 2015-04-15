// Includes
#include <TrapMutationHandler.h>
#include <MathUtils.h>

namespace xolotlCore {

void TrapMutationHandler::initialize(PSIClusterReactionNetwork *network,
		std::vector<double> grid) {
	// Add the needed reaction connectivity
	// Each (He_i)(V) cluster and I clusters are connected to He_i

	// Get all the He clusters from the network
	auto heClusters = network->getAll(heType);
	// Get all the HeV bubbles from the network
	auto bubbles = network->getAll(heVType);
	// Get the single interstitial cluster
	auto singleInterstitial = (PSICluster *) network->get(iType, 1);

	// If the single I cluster is not in the network,
	// there is no trap-mutation
	if (!singleInterstitial) {
		// Loop on the grid points
		for (int i = 0; i < grid.size(); i++) {
			// Create the list (vector) of indices at this grid point
			std::vector<int> indices;
			// And give it empty to the index vector
			indexVector.push_back(indices);
		}

		return;
	}

	// Loop on the He clusters
	for (int i = 0; i < heClusters.size(); i++) {
		// Get the cluster and its size
		auto cluster = (PSICluster *) heClusters[i];
		int heSize = cluster->getSize();

		// The helium cluster is connected to itself
		cluster->setDissociationConnectivity(cluster->getId());

		// The single interstitial cluster is connected to He
		singleInterstitial->setDissociationConnectivity(cluster->getId());

		// Loop on the bubbles
		for (int j = 0; j < bubbles.size(); j++) {
			// Get the bubble and its composition
			auto bubble =  (PSICluster *) bubbles[j];
			auto comp = bubble->getComposition();

			// We are only interested in bubbles with one vacancy
			if (comp[vType] > 1) continue;

			// Connect with He if the number of helium in the bubble is the same
			if (comp[heType] == heSize) {
				bubble->setDissociationConnectivity(cluster->getId());
			}
		}
	}

	// Two vectors to define the modified trap-mutation: the first one is the
	// depth up to which the modified trap-mutation is valid, and the second
	// one is the minimum size of helium cluster undergoing trap-mutation at
	// that depth.
	std::vector<double> depthVec = {0.0, 0.2, 0.5, 0.7, 1.0, 1.2, 1.5, 1.7, 2.0};
	std::vector<int> sizeVec = {0, 1, 2, 3, 4, 5, 6, 7, 8};
//	std::vector<double> depthVec = {-0.1, 0.6, 0.7, 0.8, 0.9};
//	std::vector<int> sizeVec = {1, 2, 3, 4, 8};

	// Clear the vector of HeV indices created by He undergoing trap-mutation
	// at each grid point
	indexVector.clear();
	// Loop on the grid points
	for (int i = 0; i < grid.size(); i++) {
		// Get the depth
		double depth = grid[i];
		// Create the list (vector) of indices at this grid point
		std::vector<int> indices;

		// If the depth is greater than the last depth stored in depthVec
		// there is no modified trap-mutation
		if (depth > depthVec[depthVec.size() - 1]) {
			indexVector.push_back(indices);
			continue;
		}

		// Loop to determine at which depth interval we are
		int l = 0;
		for (l = 0; l < depthVec.size() - 1; l++) {
			if (depth > depthVec[l] && depth <= depthVec[l+1]) break;
		}

		// Loop on the helium clusters
		for (int j = 0; j < heClusters.size(); j++) {
			// Get the cluster and its size
			auto cluster = (PSICluster *) heClusters[j];
			int heSize = cluster->getSize();

			// Skip if the helium cluster is too small for this depth
			if (heSize <= sizeVec[l]) continue;

			// Loop on the bubbles
			for (int k = 0; k < bubbles.size(); k++) {
				// Get the bubble and its composition
				auto bubble =  (PSICluster *) bubbles[k];
				auto comp = bubble->getComposition();
				if (comp[vType] > 1) continue;
				if (comp[heType] == heSize) {
					// Add this bubble to the indices
					indices.push_back(k);
				}
			}
		}

		// Add indices to the index vector
		indexVector.push_back(indices);
	}

	// Update the bubble bursting rate
	updateTrapMutationRate(network);

	return;
}

void TrapMutationHandler::updateTrapMutationRate(PSIClusterReactionNetwork *network) {
	// Get all the He clusters from the network
	auto heClusters = network->getAll(heType);

	// Update the rate by finding the biggest rate in all the He clusters
	kMutation = 0.0;
	// Loop on the helium
	for (int j = 0; j < heClusters.size(); j++) {
		// Get the bubble and its rate
		auto cluster =  (PSICluster *) heClusters[j];
		double rate = cluster->getBiggestRate();

		// Compare it to the bursting rate
		if (kMutation < rate) kMutation = rate;
	}
	// Multiply it by 1000.0
	kMutation = 1000.0 * kMutation;

	return;
}

void TrapMutationHandler::computeTrapMutation(PSIClusterReactionNetwork *network,
		int xi, double *concOffset, double *updatedConcOffset) {
	// Get all the HeV bubbles
	auto bubbles = network->getAll(heVType);

	// Get the single interstitial cluster
	auto singleInterstitial = (PSICluster *) network->get(iType, 1);

	// Don't do anything if I is not in the network
	if (!singleInterstitial) return;

	// Get its ID
	int interstitialIndex = singleInterstitial->getId() - 1;

	// Get the pointer to list of indices at this grid point
	std::vector<int> * indices = &indexVector[xi];
	// Loop on the list
	for (int i = 0; i < indices->size(); i++) {
		// Get the stored bubble and its ID
		auto bubble = (PSICluster *) bubbles[indices->at(i)];
		int bubbleIndex = bubble->getId() - 1;

		// Get the helium cluster with the same number of He and its ID
		auto comp = bubble->getComposition();
		auto heCluster = (PSICluster *) network->get(heType, comp[heType]);
		int heIndex = heCluster->getId() - 1;

		// Get the initial concentration of helium
		double oldConc = concOffset[heIndex];

		// Update the concentrations (the bubble loses its concentration)
		updatedConcOffset[heIndex] -= kMutation * oldConc;
		updatedConcOffset[bubbleIndex] += kMutation * oldConc;
		updatedConcOffset[interstitialIndex] += kMutation * oldConc;
	}

	return;
}

int TrapMutationHandler::computePartialsForTrapMutation(
		PSIClusterReactionNetwork *network, double *val,
		int *indices, int xi) {
	// Get all the HeV bubbles
	auto bubbles = network->getAll(heVType);

	// Get the single interstitial cluster
	auto singleInterstitial = (PSICluster *) network->get(iType, 1);

	// Don't do anything if I is not in the network
	if (!singleInterstitial) return 0;

	// Get its ID
	int interstitialIndex = singleInterstitial->getId() - 1;

	// Get the pointer to list of indices at this grid point
	std::vector<int> * clusterIndices = &indexVector[xi];
	// Loop on the list
	for (int i = 0; i < clusterIndices->size(); i++) {
		// Get the stored bubble and its ID
		auto bubble = (PSICluster *) bubbles[clusterIndices->at(i)];
		int bubbleIndex = bubble->getId() - 1;

		// Get the helium cluster with the same number of He and its ID
		auto comp = bubble->getComposition();
		auto heCluster = (PSICluster *) network->get(heType, comp[heType]);
		int heIndex = heCluster->getId() - 1;

		// Set the helium cluster partial derivative
		indices[i * 3] = heIndex;
		val[i * 3] = -kMutation;

		// Set the bubble cluster partial derivative
		indices[(i * 3) + 1] = bubbleIndex;
		val[(i * 3) + 1] = kMutation;

		// Set the interstitial cluster partial derivative
		indices[(i * 3) + 2] = interstitialIndex;
		val[(i * 3) + 2] = kMutation;
	}

	return clusterIndices->size();
}

}/* end namespace xolotlCore */

