// Includes
#include "DiffusionHandler.h"

namespace xolotlCore {

void DiffusionHandler::initializeOFill(PSICluster * cluster, int size,
		PetscInt *ofill) {
	// Get the index of the cluster
	int reactantIndex = cluster->getId() - 1;
	// Set the ofill value to 1 for this cluster
	ofill[reactantIndex * size + reactantIndex] = 1;

	return;
}

void DiffusionHandler::computeDiffusion(PSICluster * cluster, double sx,
		PetscScalar *concOffset, PetscScalar *leftConcOffset,
		PetscScalar *rightConcOffset, PetscScalar *updatedConcOffset) {
	// Get the index of the cluster
	int reactantIndex = cluster->getId() - 1;
	// Get the concentrations from the previous solution
	double oldConc = concOffset[reactantIndex];
	double oldLeftConc = leftConcOffset[reactantIndex];
	double oldRightConc = rightConcOffset[reactantIndex];

	// Use a simple midpoint stencil to compute the concentration
	double conc = cluster->getDiffusionCoefficient()
			* (-2.0 * oldConc + oldLeftConc + oldRightConc) * sx;

	// Update the concentration of the cluster
	updatedConcOffset[reactantIndex] += conc;

	return;
}

void DiffusionHandler::computePartialsForDiffusion(PSICluster * cluster, PetscReal sx,
			PetscReal val[3], PetscInt row[1], PetscInt col[3], PetscInt xi,
			PetscInt xs, int size) {
		// Get the diffusion coefficient of the cluster
		double diffCoeff = cluster->getDiffusionCoefficient();
		// Compute the partial derivatives for diffusion of this cluster
		val[0] = diffCoeff * sx;
		val[1] = -2.0 * diffCoeff * sx;
		val[2] = diffCoeff * sx;

		// Get the reactant index
		int reactantIndex = cluster->getId() - 1;
		// Set the row and column indices. These indices are computed
		// by using xi, xi-1 and xi+1 and the arrays are shifted to
		// (xs+1)*size to properly account for the neighboring ghost
		// cells.
		row[0] = (xi - xs + 1) * size + reactantIndex;
		col[0] = ((xi - 1) - xs + 1) * size + reactantIndex;
		col[1] = (xi - xs + 1) * size + reactantIndex;
		col[2] = ((xi + 1 + 1) - xs) * size + reactantIndex;

		// Boundary conditions
		if (xi == 0) {
			val[0] = 0.0;
			val[1] = 0.0;
			val[2] = 0.0;
		}

	return;
}

}/* end namespace xolotlCore */
