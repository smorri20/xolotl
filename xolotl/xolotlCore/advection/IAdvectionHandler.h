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
	 * The off-diagonal part of the Jacobian is already initialized by the diffusion handler.
	 * This function initialize the list of clusters that will move through advection.
	 *
	 * @param network The network
	 */
	virtual void initialize(PSIClusterReactionNetwork *network) = 0;

	/**
	 * Compute the flux due to the advection for all the cluster,
	 * given the space parameter h and the position.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * @param network The network
	 * @param h The space parameter, here the grid step size
	 * @param pos The position on the grid
	 * @param surfacePos The index of the position on the surface
	 * @param concVector The pointer to the pointer of arrays of concentration at middle,
	 * left, and right grid points
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the advection is computed used to find the next solution
	 */
	virtual void computeAdvection(PSIClusterReactionNetwork *network, double h,
			std::vector<double> &pos, int surfacePos, double **concVector,
			double *updatedConcOffset) = 0;

	/**
	 * Compute the partials due to the advection of all the clusters given
	 * the space parameter h and the position.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * @param network The network
	 * @param h The space parameter, here the grid step size
	 * @param val The pointer to the array that will contain the values of partials
	 * for the advection
	 * @param indices The pointer to the array that will contain the indices of the
	 * advecting cluster in the network
	 * @param pos The position on the grid
	 * @param surfacePos The index of the position on the surface
	 */
	virtual void computePartialsForAdvection(PSIClusterReactionNetwork *network,
			double h, double *val, int *indices, std::vector<double> &pos,
			int surfacePos) = 0;

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
