// Includes
#include <TrapMutationHandler.h>
#include <MathUtils.h>

namespace xolotlCore {

void TrapMutationHandler::initialize(int surfacePos, PSIClusterReactionNetwork *network,
		std::vector<double> grid) {
	// Add the needed reaction (dissociation) connectivity
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

	// This method fills one vector to define the modified trap-mutation: the first value
	// correspond to the depth at which the He1 cluster undergo trap-mutation (if the value
	// is negative it means that it doesn't TM), the second value correspond to He2, etc.
	initializeDepthSize();

	// Method that will fill the index vector that is actually used during the solving steps
	initializeIndex(surfacePos, network, grid);

	// Update the bubble bursting rate
	updateTrapMutationRate(network);

	return;
}

void TrapMutationHandler::initializeIndex(int surfacePos, PSIClusterReactionNetwork *network,
		std::vector<double> grid) {
	// Clear the vector of HeV indices created by He undergoing trap-mutation
	// at each grid point
	indexVector.clear();

	// Get all the He clusters from the network
	auto heClusters = network->getAll(heType);
	// Get all the HeV bubbles from the network
	auto bubbles = network->getAll(heVType);

	// Loop on the grid points
	for (int i = 0; i < grid.size(); i++) {
		// Create the list (vector) of indices at this grid point
		std::vector<int> indices;

		// If we are on the left side of the surface there is no
		// modified trap-mutation
		if (i < surfacePos) {
			indexVector.push_back(indices);
			continue;
		}

		// Get the depth
		double depth = grid[i] - grid[surfacePos];

		// Loop on the depth vector
		for (int j = 0; j < depthVec.size(); j++) {
			// Check if a helium cluster undergo TM at this depth
			if (std::fabs(depth - depthVec[j])  < 0.01) {
				// Add the bubble of size j+1 to the indices
				// Loop on the bubbles
				for (int k = 0; k < bubbles.size(); k++) {
					// Get the bubble and its composition
					auto bubble =  (PSICluster *) bubbles[k];
					auto comp = bubble->getComposition();
					// We are interested in bubbles with only one vacancy
					if (comp[vType] > 1) continue;
					if (comp[heType] == j+1) {
						// Add this bubble to the indices
						indices.push_back(k);
					}
				}
			}
		}

//		std::cout << "Depth " << depth << " at i " << i << ": ";
//		for (int j = 0; j < indices.size(); j++) {
//			std::cout << indices[j] << " ";
//		}
//		std::cout << std::endl;

		// Add indices to the index vector
		indexVector.push_back(indices);
	}

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
	// Multiply it by 1000.0 so that trap-mutation overcomes any other reaction
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

	// Initialize the rate of the reaction
	double rate = 0.0;
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

		// Check the desorption
		if (comp[heType] == desorp.size) {
			// Get the left side rate (combination + emission)
			double totalRate = heCluster->getLeftSideRate();
			// Define the trap-mutation rate taking into account the desorption
			rate = totalRate * (1.0 - desorp.portion) / desorp.portion;
		}
		else {
			rate = kMutation;
		}

		// Update the concentrations (the helium cluster loses its concentration)
		updatedConcOffset[heIndex] -= rate * oldConc;
		updatedConcOffset[bubbleIndex] += rate * oldConc;
		updatedConcOffset[interstitialIndex] += rate * oldConc;
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
	// Initialize the rate of the reaction
	double rate = 0.0;
	// Loop on the list
	for (int i = 0; i < clusterIndices->size(); i++) {
		// Get the stored bubble and its ID
		auto bubble = (PSICluster *) bubbles[clusterIndices->at(i)];
		int bubbleIndex = bubble->getId() - 1;

		// Get the helium cluster with the same number of He and its ID
		auto comp = bubble->getComposition();
		auto heCluster = (PSICluster *) network->get(heType, comp[heType]);
		int heIndex = heCluster->getId() - 1;

		// Check the desorption
		if (comp[heType] == desorp.size) {
			// Get the left side rate (combination + emission)
			double totalRate = heCluster->getLeftSideRate();
			// Define the trap-mutation rate taking into account the desorption
			rate = totalRate * (1.0 - desorp.portion) / desorp.portion;
		}
		else {
			rate = kMutation;
		}

		// Set the helium cluster partial derivative
		indices[i * 3] = heIndex;
		val[i * 3] = -rate;

		// Set the bubble cluster partial derivative
		indices[(i * 3) + 1] = bubbleIndex;
		val[(i * 3) + 1] = rate;

		// Set the interstitial cluster partial derivative
		indices[(i * 3) + 2] = interstitialIndex;
		val[(i * 3) + 2] = rate;
	}

	return clusterIndices->size();
}

}/* end namespace xolotlCore */

