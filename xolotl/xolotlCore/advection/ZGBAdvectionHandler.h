#ifndef ZGBADVECTIONHANDLER_H
#define ZGBADVECTIONHANDLER_H

// Includes
#include "AdvectionHandler.h"
#include <MathUtils.h>

namespace xolotlCore {

/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile helium cluster.
 */
class ZGBAdvectionHandler: public AdvectionHandler {
private:

	//! The location of the GB along the Z axis
	double location;

public:

	//! The Constructor
	ZGBAdvectionHandler() : location(0.0) {}

	//! The Destructor
	~ZGBAdvectionHandler() {}

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
				double oldBottomConc = concVector[5][index]; // front
				double oldTopConc = concVector[6][index]; // back

				double conc = (3.0 * sinkStrengthVector[i] * cluster->getDiffusionCoefficient())
							* ((oldBottomConc / pow(h[2], 4)) + (oldTopConc / pow(h[2], 4)))
							/ (xolotlCore::kBoltzmann * cluster->getTemperature() * h[2]);

				// Update the concentration of the cluster
				updatedConcOffset[index] += conc;
			}
			else {
				// Get the initial concentrations
				double oldConc = concVector[0][index]; // middle
				double oldRightConc = concVector[6*(pos[2] > location) + 5*(pos[2] < location)][index]; // back or front

				// Get the a=y and b=y+h positions
				double a = abs(location - pos[2]);
				double b = abs(location - pos[2]) + h[2];

				// Compute the concentration as explained in the description of the method
				double conc = (3.0 * sinkStrengthVector[i] * cluster->getDiffusionCoefficient())
							* ((oldRightConc / pow(b, 4)) - (oldConc / pow(a, 4)))
							/ (xolotlCore::kBoltzmann * cluster->getTemperature() * h[2]);

				// Update the concentration of the cluster
				updatedConcOffset[index] += conc;
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
				val[i * 2] = (3.0 * sinkStrength * diffCoeff)
							/ (xolotlCore::kBoltzmann * cluster->getTemperature()
									* h[2] * pow(h[2], 4)); // back or front
				val[(i * 2) + 1] = val[i * 2]; // back or front
			}
			else {
				// Get the a=y and b=y+h positions
				double a = abs(location - pos[2]);
				double b = abs(location - pos[2]) + h[2];

				// Compute the partial derivatives for advection of this cluster as
				// explained in the description of this method
				val[i * 2] = -(3.0 * sinkStrength * diffCoeff)
							/ (xolotlCore::kBoltzmann * cluster->getTemperature()
									* h[2] * pow(a, 4)); // middle
				val[(i * 2) + 1] = (3.0 * sinkStrength * diffCoeff)
							/ (xolotlCore::kBoltzmann * cluster->getTemperature()
									* h[2] * pow(b, 4)); // back or front
			}
		}

		return;
	}

	/**
	 * Compute the indices that will determine where the partial derivatives will
	 * be put in the Jacobian.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * Here we consider GB in the Z direction.
	 *
	 * @param pos The position on the grid
	 * @return The indices for the position in the Jacobian
	 */
	std::vector<int> getStencilForAdvection(std::vector<double> &pos) {
		// The third index is positive if pos[2] > location
		// negative if pos[2] < location
		return {0, 0, (pos[2] > location) - (pos[2] < location) + isPointOnSink(pos)};
	}

	/**
	 * Check whether the grid point is located on the sink surface or not.
	 *
	 * @param pos The position on the grid
	 * @return True if the point is on the sink
	 */
	bool isPointOnSink(std::vector<double> &pos) {
		// Return true if pos[2] is equal to location
		return abs(location - pos[2]) < 0.001;
	}

};
//end class ZGBAdvectionHandler

} /* end namespace xolotlCore */
#endif
