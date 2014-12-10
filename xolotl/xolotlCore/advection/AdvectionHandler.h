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
	 * Initialize the off-diagonal part of the Jacobian for advection.
	 * This function initialize the list of clusters that will move through advection.
	 *
	 * @param network The network
	 * @param ofill The pointer to the array that will contain the value 1 at the indices
	 * of the advecting clusters
	 */
	virtual void initializeOFill(PSIClusterReactionNetwork *network, int *ofill) {return;}

	/**
	 * Compute the flux due to the advection for all the cluster,
	 * given the space parameter hx and the position index xi.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * If D is the diffusion coefficient, and C_r, C_m the right and middle concentration
	 * of this cluster, A the sink strength, K the Boltzmann constant, T the temperature,
	 * the value to add to the updated concentration is:
	 *
	 * [(3 * A * D) / (K * T * hx)] * [(C_r / [(xi + 1) * hx]^4) - (C_m / (xi * hx)^4)]
	 *
	 * @param network The network
	 * @param hx The space parameter, here the grid step size
	 * @param xi The index of the position on the grid
	 * @param surfacePos The index of the position on the surface
	 * @param concOffset The pointer to the array of concentration at the grid
	 * point where the advection is computed
	 * @param rightConcOffset The pointer to the array of concentration at the grid
	 * point to the right of where the advection is computed
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the advection is computed used to find the next solution
	 */
	void computeAdvection(PSIClusterReactionNetwork *network, double hx,
			int xi, int surfacePos, double *concOffset, double *rightConcOffset,
			double *updatedConcOffset);

	/**
	 * Compute the partials due to the advection of all the clusters given
	 * the space parameter hx and the position index xi.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * The partial derivative on the left grid point is the opposite of the one
	 * on this grid point which is (same notation as for the computeAdvection method)
	 *
	 * - (3 * A * D) / (K * T * hx * (xi * hx)^4)
	 *
	 * @param network The network
	 * @param hx The space parameter, here the grid step size
	 * @param val The pointer to the array that will contain the values of partials
	 * for the advection
	 * @param indices The pointer to the array that will contain the indices of the
	 * advecting cluster in the network
	 * @param xi The index of the grid point
	 * @param surfacePos The index of the position on the surface
	 */
	void computePartialsForAdvection(PSIClusterReactionNetwork *network,
			double hx, double *val, int *indices, int xi, int surfacePos);

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
