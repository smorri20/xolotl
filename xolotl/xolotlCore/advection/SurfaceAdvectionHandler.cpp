// Includes
#include "SurfaceAdvectionHandler.h"

namespace xolotlCore {

void SurfaceAdvectionHandler::computeAdvection(
		PSIClusterReactionNetwork *network, double *h, std::vector<double> &pos,
		double **concVector, double *updatedConcOffset) {
	// Get all the reactant
	auto reactants = network->getAll();
	// Get the number of advecting cluster
	int nAdvec = indexVector.size();

	// Loop on the advecting clusters
	for (int i = 0; i < nAdvec; i++) {
		// Get a specific one and its index
		auto cluster = (PSICluster *) reactants->at(indexVector[i]);
		int index = cluster->getId() - 1;

		// Get the initial concentrations
		double oldConc = concVector[0][index]; // middle
		double oldRightConc = concVector[2][index]; // right

		// Compute the concentration as explained in the description of the method
		double conc = (3.0 * sinkStrengthVector[i]
				* cluster->getDiffusionCoefficient())
				* ((oldRightConc / pow(pos[0] + h[0], 4))
						- (oldConc / pow(pos[0], 4)))
				/ (xolotlCore::kBoltzmann * cluster->getTemperature() * h[0]);

		// Update the concentration of the cluster
		updatedConcOffset[index] += conc;
	}

	return;
}

void SurfaceAdvectionHandler::computePartialsForAdvection(
		PSIClusterReactionNetwork *network, double *h, double *val,
		int *indices, std::vector<double> &pos) {
	// Get all the reactant
	auto reactants = network->getAll();
	// Get the number of advecting cluster
	int nAdvec = indexVector.size();

	// Loop on the advecting clusters
	for (int i = 0; i < nAdvec; i++) {
		// Get a specific one and its index
		auto cluster = (PSICluster *) reactants->at(indexVector[i]);
		int index = cluster->getId() - 1;
		// Get the diffusion coefficient of the cluster
		double diffCoeff = cluster->getDiffusionCoefficient();
		// Get the sink strength value
		double sinkStrength = sinkStrengthVector[i];

		// Set the cluster index that will be used by PetscSolver
		// to compute the row and column indices for the Jacobian
		indices[i] = index;

		// Compute the partial derivatives for advection of this cluster as
		// explained in the description of this method
		val[i * 2] = -(3.0 * sinkStrength * diffCoeff)
				/ (xolotlCore::kBoltzmann * cluster->getTemperature() * h[0]
						* pow(pos[0], 4)); // middle
		val[(i * 2) + 1] = (3.0 * sinkStrength * diffCoeff)
				/ (xolotlCore::kBoltzmann * cluster->getTemperature() * h[0]
						* pow(pos[0] + h[0], 4)); // right
	}

	return;
}

std::vector<int> SurfaceAdvectionHandler::getStencilForAdvection(
		std::vector<double> &pos) {
	// Always return (1, 0, 0)
	return {1, 0, 0};
}

}/* end namespace xolotlCore */
