// Includes
#include "XGBAdvectionHandler.h"

namespace xolotlCore {

void XGBAdvectionHandler::initialize(PSIClusterReactionNetwork *network,
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

void XGBAdvectionHandler::computeAdvection(PSIClusterReactionNetwork *network,
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
			double oldLeftConc = concVector[1][index]; // left
			double oldRightConc = concVector[2][index]; // right

			double conc = (3.0 * sinkStrengthVector[i]
					* cluster->getDiffusionCoefficient())
					* ((oldLeftConc / pow(hxLeft, 5))
							+ (oldRightConc / pow(hxRight, 5)))
					/ (xolotlCore::kBoltzmann * cluster->getTemperature());

			// Update the concentration of the cluster
			updatedConcOffset[index] += conc;

			// Removing the diffusion of this cluster
			double oldConc = concVector[0][index]; // middle
			conc = -cluster->getDiffusionCoefficient() * 2.0 * oldConc
					/ (hxLeft * hxRight);
			// Update the concentration of this cluster
			updatedConcOffset[index] -= conc;

			// In 2D/3D there are more things to remove for the diffusion
			if (dimension > 1) {
				double oldBottomConc = concVector[3][index]; // bottom
				double oldTopConc = concVector[4][index]; // top
				double sy = 1.0 / (hy * hy);
				conc = cluster->getDiffusionCoefficient() * sy
						* (oldBottomConc + oldTopConc - 2.0 * oldConc);
				// Update the concentration of this cluster
				updatedConcOffset[index] -= conc;

				// And more things in 3D
				if (dimension == 3) {
					double oldFrontConc = concVector[5][index]; // front
					double oldBackConc = concVector[6][index]; // back
					double sz = 1.0 / (hz * hz);
					conc = cluster->getDiffusionCoefficient() * sz
							* (oldFrontConc + oldBackConc - 2.0 * oldConc);
					// Update the concentration of this cluster
					updatedConcOffset[index] -= conc;
				}
			}

		}
		// Here we are NOT on the GB sink
		else {
			// Get the initial concentrations
			double oldConc = concVector[0][index]; // middle
			double oldRightConc = concVector[2 * (pos[0] > location)
					+ 1 * (pos[0] < location)][index]; // left or right

			// Get the a=d and b=d+h positions
			double a = fabs(location - pos[0]);
			double b = fabs(location - pos[0]) + hxRight * (pos[0] > location)
					+ hxLeft * (pos[0] < location);

			// Compute the concentration as explained in the description of the method
			double conc = (3.0 * sinkStrengthVector[i]
					* cluster->getDiffusionCoefficient())
					* ((oldRightConc / pow(b, 4)) - (oldConc / pow(a, 4)))
					/ (xolotlCore::kBoltzmann * cluster->getTemperature()
							* (hxRight * (pos[0] > location)
									+ hxLeft * (pos[0] < location)));

			// Update the concentration of the cluster
			updatedConcOffset[index] += conc;

			// If the position is next to the advection sink location
			// we must remove the diffusion of this cluster
			std::vector<double> newPos = { 0.0, 0.0, 0.0 };
			newPos[0] = pos[0] + hxRight;
			if (isPointOnSink(newPos)) {
				// We are on the left side of the sink location
				// So we won't receive the diffusion from the right side
				oldConc = concVector[2][index]; // right
				conc = cluster->getDiffusionCoefficient() * oldConc * 2.0
						/ (hxRight * (hxLeft + hxRight));
				// Update the concentration of this cluster
				updatedConcOffset[index] -= conc;
			}
			newPos[0] = pos[0] - hxLeft;
			if (isPointOnSink(newPos)) {
				// We are on the right side of the sink location
				// So we won't receive the diffusion from the left side
				oldConc = concVector[1][index]; // left
				conc = cluster->getDiffusionCoefficient() * oldConc * 2.0
						/ (hxLeft * (hxLeft + hxRight));
				// Update the concentration of this cluster
				updatedConcOffset[index] -= conc;
			}
		}
	}

	return;
}

void XGBAdvectionHandler::computePartialsForAdvection(
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
			// 1D case
			if (dimension == 1) {
				val[i * 3] = (3.0 * sinkStrength * diffCoeff)
						/ (xolotlCore::kBoltzmann * cluster->getTemperature()
								* pow(hxLeft, 5)); // left
				val[(i * 3) + 1] = (3.0 * sinkStrength * diffCoeff)
						/ (xolotlCore::kBoltzmann * cluster->getTemperature()
								* pow(hxRight, 5)); // right

				// Removing the diffusion on the middle point
				val[(i * 3) + 2] = diffCoeff * 2.0 / (hxLeft * hxRight); // middle
			}
			// 2D case
			else if (dimension == 2) {
				val[i * 5] = (3.0 * sinkStrength * diffCoeff)
						/ (xolotlCore::kBoltzmann * cluster->getTemperature()
								* pow(hxLeft, 5)); // left
				val[(i * 5) + 1] = (3.0 * sinkStrength * diffCoeff)
						/ (xolotlCore::kBoltzmann * cluster->getTemperature()
								* pow(hxRight, 5)); // right

				// Removing the diffusion on the middle point
				double sy = 1.0 / (hy * hy);
				val[(i * 5) + 2] = diffCoeff * 2.0
						* (1.0 / (hxLeft * hxRight) + sy); // middle
				val[(i * 5) + 3] = -diffCoeff * sy; // bottom
				val[(i * 5) + 4] = val[(i * 5) + 3]; // top
			}
			// 3D case
			else if (dimension == 3) {
				val[i * 7] = (3.0 * sinkStrength * diffCoeff)
						/ (xolotlCore::kBoltzmann * cluster->getTemperature()
								* pow(hxLeft, 5)); // left
				val[(i * 7) + 1] = (3.0 * sinkStrength * diffCoeff)
						/ (xolotlCore::kBoltzmann * cluster->getTemperature()
								* pow(hxRight, 5)); // right

				// Removing the diffusion on the middle point
				double sy = 1.0 / (hy * hy);
				double sz = 1.0 / (hz * hz);
				val[(i * 7) + 2] = diffCoeff * 2.0
						* (1.0 / (hxLeft * hxRight) + sy + sz); // middle
				val[(i * 7) + 3] = -diffCoeff * sy; // bottom
				val[(i * 7) + 4] = val[(i * 7) + 3]; // top
				val[(i * 7) + 5] = -diffCoeff * sz; // front
				val[(i * 7) + 6] = val[(i * 7) + 5]; // back
			}
		}
		// Here we are NOT on the GB sink
		else {
			// Get the a=d and b=d+h positions
			double a = fabs(location - pos[0]);
			double b = fabs(location - pos[0]) + hxRight * (pos[0] > location)
					+ hxLeft * (pos[0] < location);

			// If the position is next to the advection sink location
			// we must remove the diffusion of this cluster
			std::vector<double> newPosA = { 0.0, 0.0, 0.0 };
			std::vector<double> newPosB = { 0.0, 0.0, 0.0 };
			newPosA[0] = pos[0] + hxRight;
			newPosB[0] = pos[0] - hxLeft;
			if (isPointOnSink(newPosA) || isPointOnSink(newPosB)) {
				// Compute the partial derivatives for advection of this cluster as
				// explained in the description of this method
				val[i * 3] = -(3.0 * sinkStrength * diffCoeff)
						/ (xolotlCore::kBoltzmann * cluster->getTemperature()
								* pow(a, 4)
								* (hxRight * (pos[0] > location)
										+ hxLeft * (pos[0] < location))); // middle
				val[(i * 3) + 1] = (3.0 * sinkStrength * diffCoeff)
						/ (xolotlCore::kBoltzmann * cluster->getTemperature()
								* pow(b, 4)
								* (hxRight * (pos[0] > location)
										+ hxLeft * (pos[0] < location))); // left or right

				// Remove the diffusion
				if (isPointOnSink(newPosA)) {
					val[(i * 3) + 2] = -diffCoeff * 2.0
							/ (hxRight * (hxLeft + hxRight));
				} else {
					val[(i * 3) + 2] = -diffCoeff * 2.0
							/ (hxLeft * (hxLeft + hxRight));
				}
			} else {
				// Compute the partial derivatives for advection of this cluster as
				// explained in the description of this method
				val[i * 2] = -(3.0 * sinkStrength * diffCoeff)
						/ (xolotlCore::kBoltzmann * cluster->getTemperature()
								* pow(a, 4)
								* (hxRight * (pos[0] > location)
										+ hxLeft * (pos[0] < location))); // middle
				val[(i * 2) + 1] = (3.0 * sinkStrength * diffCoeff)
						/ (xolotlCore::kBoltzmann * cluster->getTemperature()
								* pow(b, 4)
								* (hxRight * (pos[0] > location)
										+ hxLeft * (pos[0] < location))); // left or right
			}
		}
	}

	return;
}

std::vector<int> XGBAdvectionHandler::getStencilForAdvection(
		std::vector<double> &pos) {
	// The first index is positive by convention if we are on the sink
	if (isPointOnSink(pos))
		return {1, 0, 0};
	// The first index is positive if pos[0] > location
	// negative if pos[0] < location
	return {(pos[0] > location) - (pos[0] < location), 0, 0};
}

bool XGBAdvectionHandler::isPointOnSink(std::vector<double> &pos) {
	// Return true if pos[0] is equal to location
	return fabs(location - pos[0]) < 0.001;
}

}/* end namespace xolotlCore */
