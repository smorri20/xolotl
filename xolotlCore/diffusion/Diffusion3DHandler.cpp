// Includes
#include "Diffusion3DHandler.h"

namespace xolotlCore {

void Diffusion3DHandler::computeDiffusion(PSIClusterReactionNetwork *network,
		double s, double **concVector, double *updatedConcOffset) {
	// Get all the reactant
	auto reactants = network->getAll();
	// Get the number of diffusing cluster
	int nDiff = indexVector.size();

	// Get the number of degrees of freedom which is the size of the network
	int dof = reactants->size();

	// Loop on them
	for (int i = 0; i < nDiff; i++) {
		// Get the diffusing cluster
		auto cluster = (PSICluster *) reactants->at(indexVector[i]);
		// Get the index of the cluster
		int index = cluster->getId() - 1;

		// Get the initial concentrations
		double oldConc = concVector[0][index];
		double oldLeftConc = concVector[1][index];
		double oldRightConc = concVector[2][index];
		double oldTopConc = concVector[3][index];
		double oldBottomConc = concVector[4][index];
		double oldFrontConc = concVector[5][index];
		double oldBackConc = concVector[6][index];

		// Use a simple midpoint stencil to compute the concentration
		double conc = cluster->getDiffusionCoefficient()
						* (-6.0 * oldConc + oldLeftConc + oldRightConc + oldTopConc
								+ oldBottomConc + oldFrontConc + oldBackConc) * s;

		// Update the concentration of the cluster
		updatedConcOffset[index] += conc;
	}

	return;
}

void Diffusion3DHandler::computePartialsForDiffusion(
		PSIClusterReactionNetwork *network,
		double s, double *val, int *indices) {
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
		// for the middle, left, right, top, bottom, front, and back grid point
		val[i * 7] = -6.0 * diffCoeff * s;
		val[(i * 7) + 1] = diffCoeff * s;
		val[(i * 7) + 2] = diffCoeff * s;
		val[(i * 7) + 3] = diffCoeff * s;
		val[(i * 7) + 4] = diffCoeff * s;
		val[(i * 7) + 5] = diffCoeff * s;
		val[(i * 7) + 6] = diffCoeff * s;
	}

	return;
}

}/* end namespace xolotlCore */
