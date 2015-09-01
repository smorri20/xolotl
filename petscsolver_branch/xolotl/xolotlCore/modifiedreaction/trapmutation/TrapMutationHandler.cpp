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
	// Get the double interstitial cluster
	auto doubleInterstitial = (PSICluster *) network->get(iType, 2);
	// Get the triple interstitial cluster
	auto tripleInterstitial = (PSICluster *) network->get(iType, 3);

	// If the I clusters are not in the network,
	// there is no trap-mutation
	if (!singleInterstitial || !doubleInterstitial || !tripleInterstitial) {
		// Clear the vector of HeV indices created by He undergoing trap-mutation
		// at each grid point
		indexVector.clear();

		// Loop on the grid points
		for (int i = 0; i < grid.size(); i++) {
			// Create the list (vector) of indices at this grid point
			std::vector<int> indices;
			// And give it empty to the index vector
			indexVector.push_back(indices);
		}

		// Inform the user
		std::cout << "The modified trap-mutation won't happen because "
				"the interstitial clusters are missing." << std::endl;

		return;
	}

	// Loop on the He clusters
	for (int i = 0; i < heClusters.size(); i++) {
		// Get the cluster and its size
		auto cluster = (PSICluster *) heClusters[i];
		int heSize = cluster->getSize();

		// The helium cluster is connected to itself
		cluster->setDissociationConnectivity(cluster->getId());

		// The single and double interstitial clusters are connected to He
		singleInterstitial->setDissociationConnectivity(cluster->getId());
		doubleInterstitial->setDissociationConnectivity(cluster->getId());
		tripleInterstitial->setDissociationConnectivity(cluster->getId());

		// Loop on the bubbles
		for (int j = 0; j < bubbles.size(); j++) {
			// Get the bubble and its composition
			auto bubble =  (PSICluster *) bubbles[j];
			auto comp = bubble->getComposition();

			// We are only interested in bubbles with one, two, or three vacancies
			if (comp[vType] > 3) continue;

			// Connect with He if the number of helium in the bubble is the same
			if (comp[heType] == heSize) {
				bubble->setDissociationConnectivity(cluster->getId());
			}
		}
	}

	// This method fills two vectors to define the modified trap-mutation: for the first one,
	// the first value corresponds to the depth at which the He1 cluster undergo trap-mutation
	// (if the value is negative it means that it doesn't TM), the second value correspond
	// to He2, etc.; the second vector gives the size of the vacancies into which He
	// trap-mutates. Information about desorption is also initialized here.
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
		if (i <= surfacePos) {
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
					// Get the correct bubble
					if (comp[heType] == j+1 && comp[vType] == sizeVec[j]) {
						// Add this bubble to the indices
						indices.push_back(k);
					}
				}
			}
		}

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
	// Initialyze the pointers to interstitial and helium clusters and their ID
	PSICluster * iCluster, * heCluster, * bubble;
	int iIndex, heIndex, bubbleIndex;

	// Initialize the rate of the reaction
	double rate = 0.0;

	// Get the pointer to list of indices at this grid point
	std::vector<int> * indices = &indexVector[xi];
	// Loop on the list
	for (int i = 0; i < indices->size(); i++) {
		// Get the stored bubble and its ID
		bubble = (PSICluster *) bubbles[indices->at(i)];
		bubbleIndex = bubble->getId() - 1;

		// Get the helium cluster with the same number of He and its ID
		auto comp = bubble->getComposition();
		heCluster = (PSICluster *) network->get(heType, comp[heType]);
		heIndex = heCluster->getId() - 1;

		// Get the interstitial cluster with the same number of I as the number
		// of vacancies in the bubble and its ID
		iCluster = (PSICluster *) network->get(iType, comp[vType]);
		iIndex = iCluster->getId() - 1;

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
		updatedConcOffset[iIndex] += rate * oldConc;
	}

	return;
}

int TrapMutationHandler::computePartialsForTrapMutation(
		PSIClusterReactionNetwork *network, double *val,
		int *indices, int xi) {
	// Get all the HeV bubbles
	auto bubbles = network->getAll(heVType);
	// Initialyze the pointers to interstitial and helium clusters and their ID
	PSICluster * iCluster, * heCluster, * bubble;
	int iIndex, heIndex, bubbleIndex;

	// Initialize the rate of the reaction
	double rate = 0.0;

	// Get the pointer to list of indices at this grid point
	std::vector<int> * clusterIndices = &indexVector[xi];
	// Loop on the list
	for (int i = 0; i < clusterIndices->size(); i++) {
		// Get the stored bubble and its ID
		bubble = (PSICluster *) bubbles[clusterIndices->at(i)];
		bubbleIndex = bubble->getId() - 1;

		// Get the helium cluster with the same number of He and its ID
		auto comp = bubble->getComposition();
		heCluster = (PSICluster *) network->get(heType, comp[heType]);
		heIndex = heCluster->getId() - 1;

		// Get the interstitial cluster with the same number of I as the number
		// of vacancies in the bubble and its ID
		iCluster = (PSICluster *) network->get(iType, comp[vType]);
		iIndex = iCluster->getId() - 1;

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
		indices[(i * 3) + 2] = iIndex;
		val[(i * 3) + 2] = rate;
	}

	return clusterIndices->size();
}

}/* end namespace xolotlCore */

