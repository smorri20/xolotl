#ifndef IBUBBLEBURSTINGHANDLER_H
#define IBUBBLEBURSTINGHANDLER_H

// Includes
#include <PSICluster.h>
#include <IReactionNetwork.h>
#include <memory>

namespace xolotlCore {

/**
 * Realizations of this interface are responsible for all the physical parts
 * for the bursting of HeV bubbles. A HeV bubble bursts when it is close to the surface,
 * and looses all its helium atoms. The solver call these methods to handle
 * the bubble bursting.
 */
class IBubbleBurstingHandler {

public:

	/**
	 * The destructor
	 */
	virtual ~IBubbleBurstingHandler() {}

	/**
	 * The initialize method has to add connectivity between the V clusters and HeV clusters
	 * of same number of V. It must also initialize the rates of the reactions and define
	 * which bubbles can burst at each grid point by calling initializeIndex.
	 *
	 * @param surfacePos The index of the position of the surface
	 * @param network The network
	 * @param grid The grid on the x axis
	 */
	virtual void initialize(int surfacePos, IReactionNetwork *network,
			std::vector<double> grid) = 0;

	/**
	 * This method defines which bursting is allowed at each grid point.
	 *
	 * @param surfacePos The index of the position of the surface
	 * @param network The network
	 * @param grid The grid on the x axis
	 */
	virtual void initializeIndex(int surfacePos, IReactionNetwork *network,
			std::vector<double> grid) = 0;

	/**
	 * This method updates the rate for the bubble bursting if rates changed in the network,
	 * it should be called when the temperature changes for instance.
	 *
	 * @param network The network
	 */
	virtual void updateBurstingRate(IReactionNetwork *network) = 0;

	/**
	 * Compute the flux due to the bubble bursting for all the cluster,
	 * given the position index xi.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * @param network The network
	 * @param xi The index of the position on the grid
	 * @param concOffset The pointer to the array of concentration at the grid
	 * point where the bursting is computed
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the bursting is computed used to find the next solution
	 */
	virtual void computeBursting(IReactionNetwork *network,
			int xi, double *concOffset, double *updatedConcOffset) = 0;

	/**
	 * Compute the partials due to the bubble bursting for all the clusters given
	 * the position index xi. Returns the number of bubbles that can possibly burst at
	 * this grid point.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * @param network The network
	 * @param val The pointer to the array that will contain the values of partials
	 * for the bursting
	 * @param indices The pointer to the array that will contain the indices of the clusters
	 * @param xi The index of the grip point
	 *
	 * @return The number of bubbles that can burst at this grid point
	 */
	virtual int computePartialsForBursting(IReactionNetwork *network,
			double *val, int *indices, int xi) = 0;

};
//end class IBubbleBurstingHandler

} /* namespace xolotlCore */
#endif
