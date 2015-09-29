// Includes
#include "ZGBAdvectionHandler.h"

namespace xolotlCore {

void ZGBAdvectionHandler::initialize(PSIClusterReactionNetwork *network,
		int *ofill) {
	// Get all the reactants and their number
	auto reactants = network->getAll();
	int size = reactants->size();

	// Clear the index and sink strength vectors
	indexVector.clear();
	sinkStrengthVector.clear();

	// Loop on all the reactants
	for (int i = 0; i < size; i++) {
		// Get the i-th cluster
		auto cluster = (PSICluster *) reactants->at(i);
		// Get its diffusion coefficient
		double diffFactor = cluster->getDiffusionFactor();

		// Don't do anything if the diffusion factor is 0.0
		if (xolotlCore::equal(diffFactor, 0.0))
			continue;

		// Keep only the helium clusters
		if (cluster->getType() != heType)
			continue;

		// Get its size
		int heSize = cluster->getSize();

		// Switch on the size to get the sink strength (in eV.nm3)
		double sinkStrength = 0.0;
		switch (heSize) {
		case 1:
			sinkStrength = 0.54e-3;
			break;
		case 2:
			sinkStrength = 1.01e-3;
			break;
		case 3:
			sinkStrength = 3.03e-3;
			break;
		case 4:
			sinkStrength = 3.93e-3;
			break;
		case 5:
			sinkStrength = 7.24e-3;
			break;
		case 6:
			sinkStrength = 10.82e-3;
			break;
		case 7:
			sinkStrength = 19.26e-3;
			break;
		}

		// If the sink strength is still 0.0, this cluster is not advecting
		if (xolotlCore::equal(sinkStrength, 0.0))
			continue;

		// Add its index (i) to the vector of indices
		indexVector.push_back(i);

		// Add the sink strength to the vector
		sinkStrengthVector.push_back(sinkStrength);

		// Set the off-diagonal part for the Jacobian to 1
		// Get its id
		int index = cluster->getId() - 1;
		// Set the ofill value to 1 for this cluster
		ofill[index * size + index] = 1;
	}

	return;
}

void ZGBAdvectionHandler::computeAdvection(PSIClusterReactionNetwork *network,
		std::vector<double> &pos, double **concVector,
		double *updatedConcOffset, double hxLeft, double hxRight, double hy,
		double hz) {
	// Get all the reactant
	auto reactants = network->getAll();
	// Get the number of advecting cluster
	int nAdvec = indexVector.size();

	// Loop on the advecting clusters
	for (int i = 0; i < nAdvec; i++) {
		// Get a specific one and its index
		auto cluster = (PSICluster *) reactants->at(indexVector[i]);
		int index = cluster->getId() - 1;

		// If we are on the sink, the behavior is not the same
		// Both sides are giving their concentrations to the center
		if (isPointOnSink(pos)) {
			double oldFrontConc = concVector[5][index]; // front
			double oldBackConc = concVector[6][index]; // back

			double conc = (3.0 * sinkStrengthVector[i]
					* cluster->getDiffusionCoefficient())
					* ((oldFrontConc / pow(hz, 5)) + (oldBackConc / pow(hz, 5)))
					/ (xolotlCore::kBoltzmann * cluster->getTemperature());

			// Update the concentration of the cluster
			updatedConcOffset[index] += conc;

			// Removing the diffusion of this cluster
			double oldConc = concVector[0][index]; // middle
			double oldLeftConc = concVector[1][index]; // left
			double oldRightConc = concVector[2][index]; // right
			double oldBottomConc = concVector[3][index]; // bottom
			double oldTopConc = concVector[4][index]; // top
			double sy = 1.0 / (hy * hy);
			double sz = 1.0 / (hz * hz);
			conc = cluster->getDiffusionCoefficient()
					* (2.0
							* (oldLeftConc + (hxLeft / hxRight) * oldRightConc
									- (1.0 + (hxLeft / hxRight)) * oldConc)
							/ (hxLeft * (hxLeft + hxRight))
							+ (oldBottomConc + oldTopConc) * sy
							- 2.0 * oldConc * (sy + sz));
			// Update the concentration of this cluster
			updatedConcOffset[index] -= conc;
		}
		// Here we are NOT on the GB sink
		else {
			// Get the initial concentrations
			double oldConc = concVector[0][index]; // middle
			double oldRightConc = concVector[6 * (pos[2] > location)
					+ 5 * (pos[2] < location)][index]; // back or front

			// Get the a=d and b=d+h positions
			double a = fabs(location - pos[2]);
			double b = fabs(location - pos[2]) + hz;

			// Compute the concentration as explained in the description of the method
			double conc = (3.0 * sinkStrengthVector[i]
					* cluster->getDiffusionCoefficient())
					* ((oldRightConc / pow(b, 4)) - (oldConc / pow(a, 4)))
					/ (xolotlCore::kBoltzmann * cluster->getTemperature() * hz);

			// Update the concentration of the cluster
			updatedConcOffset[index] += conc;

			// If the position is next to the advection sink location
			// we must remove the diffusion of this cluster
			std::vector<double> newPos = { 0.0, 0.0, 0.0 };
			newPos[2] = pos[2] + hz;
			if (isPointOnSink(newPos)) {
				// We are on the front side of the sink location
				// So we won't receive the diffusion from the back side
				oldConc = concVector[6][index]; // back
				double sz = 1.0 / (hz * hz);
				conc = cluster->getDiffusionCoefficient() * oldConc * sz;
				// Update the concentration of this cluster
				updatedConcOffset[index] -= conc;
			}
			newPos[2] = pos[2] - hz;
			if (isPointOnSink(newPos)) {
				// We are on the bakc side of the sink location
				// So we won't receive the diffusion from the front side
				oldConc = concVector[5][index]; // front
				double sz = 1.0 / (hz * hz);
				conc = cluster->getDiffusionCoefficient() * oldConc * sz;
				// Update the concentration of this cluster
				updatedConcOffset[index] -= conc;
			}
		}
	}

	return;
}

void ZGBAdvectionHandler::computePartialsForAdvection(
		PSIClusterReactionNetwork *network, double *val, int *indices,
		std::vector<double> &pos, double hxLeft, double hxRight, double hy,
		double hz) {
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

		// If we are on the sink, the partial derivatives are not the same
		// Both sides are giving their concentrations to the center
		if (isPointOnSink(pos)) {
			val[i * 7] = (3.0 * sinkStrength * diffCoeff)
					/ (xolotlCore::kBoltzmann * cluster->getTemperature()
							* pow(hz, 5)); // back or front
			val[(i * 7) + 1] = val[i * 7]; // back or front

			// Removing the diffusion on the middle point
			double sy = 1.0 / (hy * hy);
			double sz = 1.0 / (hz * hz);
			val[(i * 7) + 2] = diffCoeff * 2.0
					* ((1.0 / (hxLeft * hxRight)) + sy + sz); // middle
			val[(i * 7) + 3] = -diffCoeff * 2.0 / (hxLeft * (hxLeft + hxRight)); // left
			val[(i * 7) + 4] = -diffCoeff * 2.0
					/ (hxRight * (hxLeft + hxRight)); // right
			val[(i * 7) + 5] = -diffCoeff * sy; // bottom
			val[(i * 7) + 6] = val[(i * 7) + 5]; // top
		}
		// Here we are NOT on the GB sink
		else {
			// Get the a=d and b=d+h positions
			double a = fabs(location - pos[2]);
			double b = fabs(location - pos[2]) + hz;

			// If we are on a grid point just next to the sink location
			// We have to remove the diffusion received from the sink location
			std::vector<double> newPosA = { 0.0, 0.0, 0.0 };
			newPosA[2] = pos[2] - hz;
			std::vector<double> newPosB = { 0.0, 0.0, 0.0 };
			newPosB[2] = pos[2] + hz;
			if (isPointOnSink(newPosA) || isPointOnSink(newPosB)) {
				// Compute the partial derivatives for advection of this cluster as
				// explained in the description of this method
				val[i * 3] = -(3.0 * sinkStrength * diffCoeff)
						/ (xolotlCore::kBoltzmann * cluster->getTemperature()
								* hz * pow(a, 4)); // middle
				val[(i * 3) + 1] = (3.0 * sinkStrength * diffCoeff)
						/ (xolotlCore::kBoltzmann * cluster->getTemperature()
								* hz * pow(b, 4)); // back or front

				// Remove the diffusion
				double sz = 1.0 / (hz * hz);
				val[(i * 3) + 2] = -diffCoeff * sz; // opposite of back or front
			} else {
				// Compute the partial derivatives for advection of this cluster as
				// explained in the description of this method
				val[i * 2] = -(3.0 * sinkStrength * diffCoeff)
						/ (xolotlCore::kBoltzmann * cluster->getTemperature()
								* hz * pow(a, 4)); // middle
				val[(i * 2) + 1] = (3.0 * sinkStrength * diffCoeff)
						/ (xolotlCore::kBoltzmann * cluster->getTemperature()
								* hz * pow(b, 4)); // back or front
			}
		}
	}

	return;
}

std::vector<int> ZGBAdvectionHandler::getStencilForAdvection(
		std::vector<double> &pos) {
	// The third index is positive by convention if we are on the sink
	if (isPointOnSink(pos))
		return {0, 0, 1};
	// The third index is positive if pos[2] > location
	// negative if pos[2] < location
	return {0, 0, (pos[2] > location) - (pos[2] < location)};
}

bool ZGBAdvectionHandler::isPointOnSink(std::vector<double> &pos) {
	// Return true if pos[2] is equal to location
	return fabs(location - pos[2]) < 0.001;
}

}/* end namespace xolotlCore */
