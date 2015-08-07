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
public:

	//! The Constructor
	YGBAdvectionHandler() {}

	//! The Destructor
	~YGBAdvectionHandler() {}

	/**
	 * This function initialize the list of clusters that will move through advection for
	 * grain boundaries.
	 *
	 * @param network The network
	 * @param ofill The pointer to the array that will contain the value 1 at the indices
	 * of the advecting clusters
	 */
	void initialize(PSIClusterReactionNetwork *network, int *ofill) {
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

			// Set the off-diagonal part for the Jacobian to 1
			// Get its id
			int index = cluster->getId() - 1;
			// Set the ofill value to 1 for this cluster
			ofill[index * size + index] = 1;
		}

		return;
	}

	/**
	 * Compute the flux due to the advection for all the helium clusters,
	 * given the space parameters and the position.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * @param network The network
	 * @param pos The position on the grid
	 * @param concVector The pointer to the pointer of arrays of concentration at middle,
	 * left, and right grid points
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the advection is computed used to find the next solution
	 * @param hxLeft The step size on the left side of the point in the x direction
	 * @param hxRight The step size on the right side of the point in the x direction
	 * @param hy The step size in the y direction
	 * @param hz The step size in the z direction
	 */
	void computeAdvection(PSIClusterReactionNetwork *network,
			std::vector<double> &pos, double **concVector, double *updatedConcOffset,
			double hxLeft, double hxRight, double hy = 0.0, double hz = 0.0) {
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
							* ((oldBottomConc / pow(hy, 5)) + (oldTopConc / pow(hy, 5)))
							/ (xolotlCore::kBoltzmann * cluster->getTemperature());

				// Update the concentration of the cluster
				updatedConcOffset[index] += conc;
			}
			else {
				// Get the initial concentrations
				double oldConc = concVector[0][index]; // middle
				double oldRightConc = concVector[4*(pos[1] > location) + 3*(pos[1] < location)][index]; // top or bottom

				// Get the a=y and b=y+h positions
				double a = abs(location - pos[1]);
				double b = abs(location - pos[1]) + hy;

				// Compute the concentration as explained in the description of the method
				double conc = (3.0 * sinkStrengthVector[i] * cluster->getDiffusionCoefficient())
							* ((oldRightConc / pow(b, 4)) - (oldConc / pow(a, 4)))
							/ (xolotlCore::kBoltzmann * cluster->getTemperature() * hy);

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
	 * @param val The pointer to the array that will contain the values of partials
	 * for the advection
	 * @param indices The pointer to the array that will contain the indices of the
	 * advecting cluster in the network
	 * @param pos The position on the grid
	 * @param hxLeft The step size on the left side of the point in the x direction
	 * @param hxRight The step size on the right side of the point in the x direction
	 * @param hy The step size in the y direction
	 * @param hz The step size in the z direction
	 */
	void computePartialsForAdvection(PSIClusterReactionNetwork *network,
			double *val, int *indices, std::vector<double> &pos,
			double hxLeft, double hxRight, double hy = 0.0, double hz = 0.0){
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
									* pow(hy, 5)); // top or bottom
				val[(i * 2) + 1] = val[i * 2]; // top or bottom
			}
			else {
				// Get the a=y and b=y+h positions
				double a = abs(location - pos[1]);
				double b = abs(location - pos[1]) + hy;

				// Compute the partial derivatives for advection of this cluster as
				// explained in the description of this method
				val[i * 2] = -(3.0 * sinkStrength * diffCoeff)
							/ (xolotlCore::kBoltzmann * cluster->getTemperature()
									* hy * pow(a, 4)); // middle
				val[(i * 2) + 1] = (3.0 * sinkStrength * diffCoeff)
							/ (xolotlCore::kBoltzmann * cluster->getTemperature()
									* hy * pow(b, 4)); // top or bottom
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
