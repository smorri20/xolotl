#ifndef DUMMYBUBBLEBURSTINGHANDLER_H
#define DUMMYBUBBLEBURSTINGHANDLER_H

// Includes
#include <BubbleBurstingHandler.h>

namespace xolotlCore {

/**
 * This class realizes the IBubbleBursting interface, responsible for all the
 * physical parts for the bursting of HeV bubbles. Here it is a dummy class,
 * it won't do anything.
 */
class DummyBubbleBurstingHandler: public BubbleBurstingHandler {

public:

	/**
	 * The constructor
	 */
	DummyBubbleBurstingHandler() : BubbleBurstingHandler() {}

	/**
	 * The destructor
	 */
	~DummyBubbleBurstingHandler() {}

	/**
	 * The initialize method only clears the indexVector to be sure that no bubble
	 * is bursting.
	 *
	 * @param network The network
	 * @param hx The grid step size
	 * @param nGrid The number of points on the grid
	 * @param surfacePos The index of the position on the surface
	 */
	void initialize(std::shared_ptr<PSIClusterReactionNetwork> network, double hx,
			int nGrid, int surfacePos) {
		// Clear the vector of HeV bubble bursting at each grid point
		indexVector.clear();

		// And don't do anything else
		return;
	}

	/**
	 * This method update the rate for the bubble bursting if the rates changed in the network,
	 * it should be called when temperature changes for instance.
	 * Reset the rate to 0.0 just in case for the dummy class
	 *
	 * @param network The network
	 */
	void updateBurstingRate(std::shared_ptr<PSIClusterReactionNetwork> network) {
		kBursting = 0.0;

		return;
	}

	/**
	 * Compute the flux due to the bubble bursting for all the cluster,
	 * given the position index xi.
	 * This method is called by the RHSFunction from the PetscSolver.
	 * Here it won't do anything to the vector of concentration because it is a dummy class.
	 *
	 * @param network The network
	 * @param xi The index of the position on the grid
	 * @param surfacePos The index of the position on the surface
	 * @param concOffset The pointer to the array of concentration at the grid
	 * point where the advection is computed
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the advection is computed used to find the next solution
	 */
	void computeBursting(std::shared_ptr<PSIClusterReactionNetwork> network,
			int xi, int surfacePos, double *concOffset, double *updatedConcOffset) {
		return;
	}

	/**
	 * Compute the partials due to the bubble bursting for all the clusters given
	 * the position index xi. Returns the number of bubbles that can possibly burst.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 * Here it won't do anything to the vector of concentration because it is a dummy class.
	 *
	 * @param network The network
	 * @param val The pointer to the array that will contain the values of partials
	 * for the advection
	 * @param row The pointer to the array that will contain the indices of the row
	 * for the Jacobian
	 * @param col The pointer to the array that will contain the indices of the columns
	 * for the Jacobian
	 * @param xi The index of the grip point
	 * @param xs The index of the first grid point on the locally owned grid
	 * @param surfacePos The index of the position on the surface
	 *
	 * @return The number of bubbles that can burst at this grid point
	 */
	int computePartialsForBursting(std::shared_ptr<PSIClusterReactionNetwork> network,
			double *val, int *row, int *col, int xi, int xs, int surfacePos) {
		return 0;
	}

};
//end class DummyBubbleBurstingHandler

} /* namespace xolotlCore */
#endif
