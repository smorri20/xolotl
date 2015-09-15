#ifndef YGBADVECTIONHANDLER_H
#define YGBADVECTIONHANDLER_H

// Includes
#include "AdvectionHandler.h"
#include <MathUtils.h>

namespace xolotlCore {

/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile helium cluster.
 */
class YGBAdvectionHandler: public AdvectionHandler {
private:

	//! The location of the GB along the Y axis
	double location;

public:

	//! The Constructor
	YGBAdvectionHandler() : location(0.0) {}

	//! The Destructor
	~YGBAdvectionHandler() {}

	/**
	 * This function initialize the list of clusters that will move through advection for
	 * grain boundaries.
	 *
	 * @param network The network
	 */
	void initialize(PSIClusterReactionNetwork *network) {
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
			if (xolotlCore::equal(diffFactor, 0.0)) continue;

			// Keep only the helium clusters
			if (cluster->getType() != heType) continue;

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
			if (xolotlCore::equal(sinkStrength, 0.0)) continue;

			// Add its index (i) to the vector of indices
			indexVector.push_back(i);

			// Add the sink strength to the vector
			sinkStrengthVector.push_back(sinkStrength);
		}

		return;
	}

	/**
	 * Set the position of the sink (location).
	 *
	 * @param pos The position of the sink
	 */
	void setPosition(double pos) {
		location = pos;

		return;
	}

	/**
	 * Compute the flux due to the advection for all the helium clusters,
	 * given the space parameters and the position.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * @param network The network
	 * @param h The space parameters in the three directions
	 * @param pos The position on the grid
	 * @param concVector The pointer to the pointer of arrays of concentration at middle,
	 * left, and right grid points
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the advection is computed used to find the next solution
	 */
	void computeAdvection(PSIClusterReactionNetwork *network, double *h,
			std::vector<double> &pos, double **concVector, double *updatedConcOffset) {
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
				double oldBottomConc = concVector[3][index]; // bottom
				double oldTopConc = concVector[4][index]; // top

				double conc = (3.0 * sinkStrengthVector[i] * cluster->getDiffusionCoefficient())
							* ((oldBottomConc / pow(h[1], 4)) + (oldTopConc / pow(h[1], 4)))
							/ (xolotlCore::kBoltzmann * cluster->getTemperature() * h[1]);

				// Update the concentration of the cluster
				updatedConcOffset[index] += conc;

				// Removing the diffusion of this cluster
				double oldConc = concVector[0][index]; // middle
				double oldLeftConc = concVector[1][index]; // left
				double oldRightConc = concVector[2][index]; // right
				double sx = 1.0 / (h[0] * h[0]);
				double sy = 1.0 / (h[1] * h[1]);
				conc = cluster->getDiffusionCoefficient()
						* ((oldLeftConc + oldRightConc) * sx
								- 2.0 * oldConc * (sx + sy));
				// Update the concentration of this cluster
				updatedConcOffset[index] -= conc;

				// In 3D there are more things to remove for diffusion
				if (dimension == 3) {
					double oldFrontConc = concVector[5][index]; // front
					double oldBackConc = concVector[6][index]; // back
					double sz = 1.0 / (h[2] * h[2]);
					conc = cluster->getDiffusionCoefficient()
							* sz * (oldFrontConc + oldBackConc - 2.0 * oldConc);
					// Update the concentration of this cluster
					updatedConcOffset[index] -= conc;
				}
			}
			else {
				// Get the initial concentrations
				double oldConc = concVector[0][index]; // middle
				double oldRightConc = concVector[4*(pos[1] > location) + 3*(pos[1] < location)][index]; // top or bottom

				// Get the a=y and b=y+h positions
				double a = abs(location - pos[1]);
				double b = abs(location - pos[1]) + h[1];

				// Compute the concentration as explained in the description of the method
				double conc = (3.0 * sinkStrengthVector[i] * cluster->getDiffusionCoefficient())
							* ((oldRightConc / pow(b, 4)) - (oldConc / pow(a, 4)))
							/ (xolotlCore::kBoltzmann * cluster->getTemperature() * h[1]);

				// Update the concentration of the cluster
				updatedConcOffset[index] += conc;

				// If the position is next to the advection sink location
				// we must remove the diffusion of this cluster
				std::vector<double> newPos = { 0.0, 0.0, 0.0 };
				newPos[1] = pos[1] + h[1];
				if (isPointOnSink(newPos)) {
					// We are on the bottom side of the sink location
					// So we won't receive the diffusion from the top side
					oldConc = concVector[4][index]; // top
					double sy = 1.0 / (h[1] * h[1]);
					conc = cluster->getDiffusionCoefficient() * oldConc * sy;
					// Update the concentration of this cluster
					updatedConcOffset[index] -= conc;
				}
				newPos[1] = pos[1] - h[1];
				if (isPointOnSink(newPos)) {
					// We are on the top side of the sink location
					// So we won't receive the diffusion from the bottom side
					oldConc = concVector[3][index]; // bottom
					double sy = 1.0 / (h[1] * h[1]);
					conc = cluster->getDiffusionCoefficient() * oldConc * sy;
					// Update the concentration of this cluster
					updatedConcOffset[index] -= conc;
				}
			}
		}

		return;
	}

	/**
	 * Compute the partials due to the advection of all the helium clusters given
	 * the space parameters and the position.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * @param network The network
	 * @param h The space parameters in the three directions
	 * @param val The pointer to the array that will contain the values of partials
	 * for the advection
	 * @param indices The pointer to the array that will contain the indices of the
	 * advecting cluster in the network
	 * @param pos The position on the grid
	 */
	void computePartialsForAdvection(PSIClusterReactionNetwork *network,
			double *h, double *val, int *indices, std::vector<double> &pos){
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
				// 2D case (there is no 1D case)
				if (dimension == 2) {
					val[i * 5] = (3.0 * sinkStrength * diffCoeff)
								/ (xolotlCore::kBoltzmann * cluster->getTemperature()
										* h[1] * pow(h[1], 4)); // top or bottom
					val[(i * 5) + 1] = val[i * 5]; // top or bottom

					// Removing the diffusion on the middle point
					double sx = 1.0 / (h[0] * h[0]);
					double sy = 1.0 / (h[1] * h[1]);
					val[(i * 5) + 2] = diffCoeff * 2.0 * (sx + sy); // middle
					val[(i * 5) + 3] = - diffCoeff * sx; // left
					val[(i * 5) + 4] = val[(i * 5) + 3]; // right
				}
				// 3D case
				else if (dimension == 3) {
					val[i * 7] = (3.0 * sinkStrength * diffCoeff)
								/ (xolotlCore::kBoltzmann * cluster->getTemperature()
										* h[1] * pow(h[1], 4)); // top or bottom
					val[(i * 7) + 1] = val[i * 7]; // top or bottom

					// Removing the diffusion on the middle point
					double sx = 1.0 / (h[0] * h[0]);
					double sy = 1.0 / (h[1] * h[1]);
					double sz = 1.0 / (h[2] * h[2]);
					val[(i * 7) + 2] = diffCoeff * 2.0 * (sx + sy + sz); // middle
					val[(i * 7) + 3] = - diffCoeff * sz; // front
					val[(i * 7) + 4] = val[(i * 7) + 3]; // back
					val[(i * 7) + 5] = - diffCoeff * sx; // left
					val[(i * 7) + 6] = val[(i * 7) + 5]; // right
				}
			}
			else {
				// Get the a=y and b=y+h positions
				double a = abs(location - pos[1]);
				double b = abs(location - pos[1]) + h[1];

				// If we are on a grid point just next to the sink location
				// We have to remove the diffusion received from the sink location
				std::vector<double> newPosA = { 0.0, 0.0, 0.0 };
				newPosA[1] = pos[1] - h[1];
				std::vector<double> newPosB = { 0.0, 0.0, 0.0 };
				newPosB[1] = pos[1] + h[1];
				if (isPointOnSink(newPosA) || isPointOnSink(newPosB)) {
					// Compute the partial derivatives for advection of this cluster as
					// explained in the description of this method
					val[i * 3] = -(3.0 * sinkStrength * diffCoeff)
									/ (xolotlCore::kBoltzmann * cluster->getTemperature()
											* h[1] * pow(a, 4)); // middle
					val[(i * 3) + 1] = (3.0 * sinkStrength * diffCoeff)
									/ (xolotlCore::kBoltzmann * cluster->getTemperature()
											* h[1] * pow(b, 4)); // top or bottom

					// Remove the diffusion
					double sy = 1.0 / (h[1] * h[1]);
					val[(i * 3) + 2] = - diffCoeff * sy; // opposite of top or bottom
				}
				else {
					// Compute the partial derivatives for advection of this cluster as
					// explained in the description of this method
					val[i * 2] = -(3.0 * sinkStrength * diffCoeff)
									/ (xolotlCore::kBoltzmann * cluster->getTemperature()
											* h[1] * pow(a, 4)); // middle
					val[(i * 2) + 1] = (3.0 * sinkStrength * diffCoeff)
									/ (xolotlCore::kBoltzmann * cluster->getTemperature()
											* h[1] * pow(b, 4)); // top or bottom
				}
			}
		}

		return;
	}

	/**
	 * Compute the indices that will determine where the partial derivatives will
	 * be put in the Jacobian.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * Here we consider GB in the Y direction.
	 *
	 * @param pos The position on the grid
	 * @return The indices for the position in the Jacobian
	 */
	std::vector<int> getStencilForAdvection(std::vector<double> &pos) {
		// The second index is positive if pos[1] > location
		// negative if pos[1] < location
		return {0, (pos[1] > location) - (pos[1] < location) + isPointOnSink(pos), 0};
	}

	/**
	 * Check whether the grid point is located on the sink surface or not.
	 *
	 * @param pos The position on the grid
	 * @return True if the point is on the sink
	 */
	bool isPointOnSink(std::vector<double> &pos) {
		// Return true if pos[1] is equal to location
		return abs(location - pos[1]) < 0.001;
	}

};
//end class YGBAdvectionHandler

} /* end namespace xolotlCore */
#endif
