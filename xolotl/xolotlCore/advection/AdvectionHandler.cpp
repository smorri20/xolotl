// Includes
#include "AdvectionHandler.h"

namespace xolotlCore {

void AdvectionHandler::initialize(std::shared_ptr<PSIClusterReactionNetwork> network) {
	// Get all the reactant
	auto reactants = network->getAll();
	int size = reactants->size();

	// Clear the index and sink strength vectors
	indexVector.clear();
	sinkStrengthVector.clear();

	// Loop on them
	for (int i = 0; i < size; i++) {
		// Get the i-th cluster
		auto cluster = (PSICluster *) reactants->at(i);
		// Get its diffusion coefficient
		double diffFactor = cluster->getDiffusionFactor();

		// Don't do anything if the diffusion factor is 0.0
		if (diffFactor == 0.0) continue;

		// Only helium clusters
		if (cluster->getType() != heType) continue;

		// Get its size
		int size = cluster->getSize();

		// Switch on it to get the sink strength (in eV.nm3)
		double sinkStrength = 0.0;
		switch (size) {
			case 1:
				sinkStrength = 2.28e-3;
				break;
			case 2:
				sinkStrength = 5.06e-3;
				break;
			case 3:
				sinkStrength = 7.26e-3;
				break;
			case 4:
				sinkStrength = 15.87e-3;
				break;
			case 5:
				sinkStrength = 16.95e-3;
				break;
			case 6:
				sinkStrength = 27.16e-3;
				break;
			case 7:
				sinkStrength = 35.56e-3;
				break;
		}

		// If the sink strength is still 0.0, this cluster is not advecting
		if (sinkStrength == 0.0) continue;

		// Add it's index (i) to the vector of indices
		indexVector.push_back(i);

		// Add the sink strength to the vector
		sinkStrengthVector.push_back(sinkStrength);
	}

	return;
}

void AdvectionHandler::computeAdvection(std::shared_ptr<PSIClusterReactionNetwork> network,
		double hx, int xi, PetscScalar *concOffset, PetscScalar *rightConcOffset,
		PetscScalar *updatedConcOffset) {
	// Get all the reactant
	auto reactants = network->getAll();
	// Get the number of advecting cluster
	int nAdvec = indexVector.size();

	// Loop on them
	for (int i = 0; i < nAdvec; i++) {
		// Get the advecting cluster
		auto cluster = (PSICluster *) reactants->at(indexVector[i]);
		// Get the index of the cluster
		int index = cluster->getId() - 1;

		// Get the initial concentrations
		double oldConc = concOffset[index];
		double oldRightConc = rightConcOffset[index];

		// Compute the concentration as explained in the description of the method
		double conc = (3.0 * sinkStrengthVector[i] * cluster->getDiffusionCoefficient())
				* ((oldRightConc / pow((xi + 1) * hx, 4)) - (oldConc / pow(xi * hx, 4)))
				/ (xolotlCore::kBoltzmann * cluster->getTemperature() * hx);

//		double conc = sinkStrengthVector[i] * (oldRightConc - oldConc) / hx;

		// Update the concentration of the cluster
		updatedConcOffset[index] += conc;
	}

	return;
}

void AdvectionHandler::computePartialsForAdvection(
		std::shared_ptr<PSIClusterReactionNetwork> network,
		double hx, PetscReal *val, PetscInt *row, PetscInt *col, PetscInt xi,
		PetscInt xs) {
	// Get all the reactant
	auto reactants = network->getAll();
	// And the size of the network
	int size = reactants->size();
	// Get the number of diffusing cluster
	int nAdvec = indexVector.size();

	// Loop on them
	for (int i = 0; i < nAdvec; i++) {
		// Get the diffusing cluster
		auto cluster = (PSICluster *) reactants->at(indexVector[i]);
		// Get the index of the cluster
		int index = cluster->getId() - 1;
		// Get the diffusion coefficient of the cluster
		double diffCoeff = cluster->getDiffusionCoefficient();
		// Get the sink strenght value
		double sinkStrength = sinkStrengthVector[i];

		// Set the row and column indices. These indices are computed
		// by using xi and xi-1, and the arrays are shifted to
		// (xs+1)*size to properly account for the neighboring ghost
		// cells.

		// Set the row index
		row[i] = (xi - xs + 1) * size + index;

		// Set the columns indices
		col[i * 2] = ((xi - 1) - xs + 1) * size + index;
		col[(i * 2) + 1] = (xi - xs + 1) * size + index;

		// Compute the partial derivatives for advection of this cluster as
		// explained in the description of this method
		val[i * 2] = (3.0 * sinkStrength * diffCoeff)
						/ (xolotlCore::kBoltzmann * cluster->getTemperature()
								* hx * pow(xi * hx, 4));
		val[(i * 2) + 1] = -(3.0 * sinkStrength * diffCoeff)
								/ (xolotlCore::kBoltzmann * cluster->getTemperature()
										* hx * pow(xi * hx, 4));

//		val[i * 2] = sinkStrength / hx;
//		val[(i * 2) + 1] = - sinkStrength / hx;
	}

	return;
}

}/* end namespace xolotlCore */
