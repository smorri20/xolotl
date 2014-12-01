// Includes
#include "DiffusionHandler.h"

namespace xolotlCore {

void DiffusionHandler::initializeOFill(std::shared_ptr<PSIClusterReactionNetwork> network,
		int *ofill) {
	// Get all the reactant
	auto reactants = network->getAll();
	int size = reactants->size();

	// Clear the index vector
	indexVector.clear();

	// Loop on them
	for (int i = 0; i < size; i++) {
		// Get the i-th cluster
		auto cluster = (PSICluster *) reactants->at(i);
		// Get its diffusion coefficient
		double diffFactor = cluster->getDiffusionFactor();

		// Don't do anything if the diffusion factor is 0.0
		if (diffFactor == 0.0) continue;

		// Add it's index (i) to the vector of indices
		indexVector.push_back(i);

		// Get its id
		int index = cluster->getId() - 1;
		// Set the ofill value to 1 for this cluster
		ofill[index * size + index] = 1;
	}

	return;
}

void DiffusionHandler::computeDiffusion(std::shared_ptr<PSIClusterReactionNetwork> network,
		double sx, double *concOffset, double *leftConcOffset,
		double *rightConcOffset, double *updatedConcOffset) {
	// Get all the reactant
	auto reactants = network->getAll();
	// Get the number of diffusing cluster
	int nDiff = indexVector.size();

	// Loop on them
	for (int i = 0; i < nDiff; i++) {
		// Get the diffusing cluster
		auto cluster = (PSICluster *) reactants->at(indexVector[i]);
		// Get the index of the cluster
		int index = cluster->getId() - 1;

		// Get the initial concentrations
		double oldConc = concOffset[index];
		double oldLeftConc = leftConcOffset[index];
		double oldRightConc = rightConcOffset[index];

		// Use a simple midpoint stencil to compute the concentration
		double conc = cluster->getDiffusionCoefficient()
				* (-2.0 * oldConc + oldLeftConc + oldRightConc) * sx;

		// Update the concentration of the cluster
		updatedConcOffset[index] += conc;
	}

	return;
}

void DiffusionHandler::computePartialsForDiffusion(
		std::shared_ptr<PSIClusterReactionNetwork> network,
		double sx, double *val, int *indices) {
	// Get all the reactant
	auto reactants = network->getAll();
	// And the size of the network
	int size = reactants->size();
	// Get the number of diffusing cluster
	int nDiff = indexVector.size();

	// Loop on them
	for (int i = 0; i < nDiff; i++) {
		// Get the diffusing cluster
		auto cluster = (PSICluster *) reactants->at(indexVector[i]);
		// Get the index of the cluster
		int index = cluster->getId() - 1;
		// Get the diffusion coefficient of the cluster
		double diffCoeff = cluster->getDiffusionCoefficient();

		// Set the cluster index, the PetscSolver will use it to compute
		// the row and column indices for the Jacobian
		indices[i] = index;

		// Compute the partial derivatives for diffusion of this cluster
		// for the left, middle, and right grid point
		val[i * 3] = diffCoeff * sx;
		val[(i * 3) + 1] = -2.0 * diffCoeff * sx;
		val[(i * 3) + 2] = diffCoeff * sx;
	}

	return;
}

}/* end namespace xolotlCore */
