// Includes
#include "BubbleBurstingHandler.h"

namespace xolotlCore {

void BubbleBurstingHandler::initialize(std::shared_ptr<PSIClusterReactionNetwork> network,
		double hx, int nGrid) {
	// Add the needed reaction connectivity
	// Each V cluster connects to every HeV clusters with the same number of V

	// Get all the V clusters from the network
	auto vClusters = network->getAll(vType);
	// Get all the HeV bubbles from the network
	auto bubbles = network->getAll(heVType);
	// Loop on the V clusters
	for (int i = 0; i < vClusters.size(); i++) {
		// Get the cluster and its size
		auto cluster = (PSICluster *) vClusters[i];
		int vSize = cluster->getSize();

		// Loop on the bubbles
		for (int j = 0; j < bubbles.size(); j++) {
			// Get the bubble and its composition
			auto bubble =  (PSICluster *) bubbles[j];
			auto comp = bubble->getComposition();

			// Connect if their vSize is the same
			if (comp[vType] == vSize) {
				cluster->setReactionConnectivity(bubble->getId());
			}
		}
	}

	// Clear the vector of HeV bubble bursting at each grid point
	indexVector.clear();
	// Loop on the grid points
	for (int i = 0; i < nGrid; i++) {
		// Compute the distance between the surface and this grid point
		double surfaceDistance = (double) i * hx;
		// Create the list (vector) of indices at this grid point
		std::vector<int> indices;

		// Loop on all the bubbles
		for (int j = 0; j < bubbles.size(); j++) {
			// Get the bubble and its radius
			auto bubble = (PSICluster *) bubbles[j];
			double radius = bubble->getReactionRadius();

			// If the radius is bigger than tha distance to the surface there is bursting
			if (radius > surfaceDistance) {
				// Add the bubble index to the list of indices
				indices.push_back(j);
			}
		}

		// Add indices to the index vector
		indexVector.push_back(indices);
	}

	return;
}

void BubbleBurstingHandler::computeBursting(std::shared_ptr<PSIClusterReactionNetwork> network,
		int xi, double *concOffset, double *updatedConcOffset) {
	// Get all the HeV bubbles
	auto bubbles = network->getAll(heVType);

	// Get the pointer to list of indices at this grid point
	std::vector<int> * indices = &indexVector[xi];
	// Loop on the list
	for (int i = 0; i < indices->size(); i++) {
		// Get the stored bubble and its ID
		auto bubble = (PSICluster *) bubbles[indices->at(i)];
		int bubbleIndex = bubble->getId() - 1;

		// Get the initial concentration
		double oldConc = concOffset[bubbleIndex];

		// Go to the next bubble if this concentration is 0.0
		if (oldConc == 0.0) continue;

		// Get the V cluster with the same number of V
		auto comp = bubble->getComposition();
		auto vCluster = (PSICluster *) network->get(vType, comp[vType]);
		// And its ID
		int vIndex = vCluster->getId() - 1;

		// Update the concentrations (the bubble loses its concentration)
		updatedConcOffset[bubbleIndex] -= kBursting * oldConc;
		updatedConcOffset[vIndex] += kBursting * oldConc;
	}

	return;
}

int BubbleBurstingHandler::computePartialsForBursting(std::shared_ptr<PSIClusterReactionNetwork> network,
		double *val, int *row, int *col, int xi, int xs) {
	// Get all the HeV bubbles
	auto bubbles = network->getAll(heVType);
	// Get the size of the network
	int size = network->getAll()->size();

	// Get the pointer to list of indices at this grid point
	std::vector<int> * indices = &indexVector[xi];
	// Loop on the list
	for (int i = 0; i < indices->size(); i++) {
		// Get the stored bubble and its ID
		auto bubble = (PSICluster *) bubbles[indices->at(i)];
		int bubbleIndex = bubble->getId() - 1;

		// Set the row and column indices. These indices are computed
		// by using xi and xi-1, and the arrays are shifted to
		// (xs+1)*size to properly account for the neighboring ghost
		// cells.
		row[i * 2] = (xi - xs + 1) * size + bubbleIndex;
		col[i] = (xi - xs + 1) * size + bubbleIndex;
		// Set its partial derivative
		val[i * 2] = -kBursting;

		// Get the V cluster with the same number of V
		auto comp = bubble->getComposition();
		auto vCluster = (PSICluster *) network->get(vType, comp[vType]);
		// And its ID
		int vIndex = vCluster->getId() - 1;
		// Set its partial derivative (the column is the same as previously)
		row[(i * 2) + 1] = (xi - xs + 1) * size + vIndex;
		val[(i * 2) + 1] = kBursting;
	}

	return indices->size();
}

}/* end namespace xolotlCore */

