#ifndef ADVECTIONHANDLER_H
#define ADVECTIONHANDLER_H

// Includes
#include "IAdvectionHandler.h"

namespace xolotlCore {

/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile cluster.
 */
class AdvectionHandler: public IAdvectionHandler {
protected:

	//! The vector containing the indices of the advecting clusters
	std::vector<int> indexVector;

	//! The vector containing the value of the sink strength (called A) of the advecting clusters
	std::vector<double> sinkStrengthVector;

public:

	//! The Constructor
	AdvectionHandler() {}

	//! The Destructor
	~AdvectionHandler() {}

	/**
	 * The off-diagonal part of the Jacobian is already initialized by the diffusion handler.
	 * This function initialize the list of clusters that will move through advection.
	 * It must be overridden by the subclasses for specific surface orientations.
	 *
	 * @param network The network
	 */
	virtual void initialize(PSIClusterReactionNetwork *network) {return;}

	/**
	 * Compute the flux due to the advection for all the cluster,
	 * given the space parameter h and the position.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * If D is the diffusion coefficient, and C_r, C_m the right and middle concentration
	 * of this cluster, A the sink strength, K the Boltzmann constant, T the temperature,
	 * the value to add to the updated concentration is:
	 *
	 * [(3 * A * D) / (K * T * h)] * [(C_r / [pos_x + h]^4) - (C_m / (pos_x)^4)]
	 *
	 * @param network The network
	 * @param h The space parameter, here the grid step size
	 * @param pos The position on the grid
	 * @param concVector The pointer to the pointer of arrays of concentration at middle,
	 * left, and right grid points
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the advection is computed used to find the next solution
	 */
	void computeAdvection(PSIClusterReactionNetwork *network, double h,
			std::vector<double> &pos, double **concVector, double *updatedConcOffset);

	/**
	 * Compute the partials due to the advection of all the clusters given
	 * the space parameter h and the position.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * The partial derivative on the left grid point is the opposite of the one
	 * on this grid point which is (same notation as for the computeAdvection method)
	 *
	 * - (3 * A * D) / (K * T * h * (pos_x)^4)
	 *
	 * @param network The network
	 * @param h The space parameter, here the grid step size
	 * @param val The pointer to the array that will contain the values of partials
	 * for the advection
	 * @param indices The pointer to the array that will contain the indices of the
	 * advecting cluster in the network
	 * @param pos The position on the grid
	 */
	void computePartialsForAdvection(PSIClusterReactionNetwork *network,
			double h, double *val, int *indices, std::vector<double> &pos);

	/**
	 * Get the total number of advecting clusters in the network.
	 *
	 * @return The number of advecting clusters
	 */
	int getNumberOfAdvecting() {return indexVector.size();}

};
//end class AdvectionHandler

} /* end namespace xolotlCore */
#endif
