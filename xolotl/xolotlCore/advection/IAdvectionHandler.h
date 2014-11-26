#ifndef IADVECTIONHANDLER_H
#define IADVECTIONHANDLER_H

// Includes
#include <PSICluster.h>
#include <PSIClusterReactionNetwork.h>
#include <memory>

namespace xolotlCore {

/**
 * Realizations of this interface are responsible for all the physical parts
 * for the advection of mobile helium cluster. The solver call these methods to handle
 * the advection.
 */
class IAdvectionHandler {

public:

	/**
	 * The destructor
	 */
	virtual ~IAdvectionHandler() {}

	/**
	 * Initialize the off-diagonal part of the Jacobian for advection.
	 * This function initialize the list of clusters that will move through advection.
	 *
	 * @param network The network
	 * @param ofill The pointer to the array that will contain the value 1 at the indices
	 * of the advecting clusters
	 */
	virtual void initialize(std::shared_ptr<PSIClusterReactionNetwork> network,
			int *ofill) = 0;

	/**
	 * Compute the flux due to the advection for all the cluster,
	 * given the space parameter hx and the position index xi.
	 * This method is called by the RHSFunction from the PetscSolver.
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
	virtual void computeAdvection(std::shared_ptr<PSIClusterReactionNetwork> network, double hx,
			int xi, int surfacePos, double *concOffset, double *rightConcOffset,
			double *updatedConcOffset) = 0;

	/**
	 * Compute the partials due to the advection of all the clusters given
	 * the space parameter hx and the position index xi.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * @param network The network
	 * @param hx The space parameter, here the grid step size
	 * @param val The pointer to the array that will contain the values of partials
	 * for the advection
	 * @param row The pointer to the array that will contain the indices of the row
	 * for the Jacobian
	 * @param col The pointer to the array that will contain the indices of the columns
	 * for the Jacobian
	 * @param xi The index of the grip point
	 * @param xs The index of the first grid point on the locally owned grid
	 * @param surfacePos The index of the position on the surface
	 */
	virtual void computePartialsForAdvection(std::shared_ptr<PSIClusterReactionNetwork> network,
			double hx, double *val, int *row, int *col, int xi,
			int xs, int surfacePos) = 0;

	/**
	 * Get the total number of advecting clusters in the network.
	 *
	 * @return The number of advecting clusters
	 */
	virtual int getNumberOfAdvecting() = 0;

};
//end class IAdvectionHandler

} /* namespace xolotlCore */
#endif
