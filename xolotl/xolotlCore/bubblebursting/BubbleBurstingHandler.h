#ifndef BUBBLEBURSTINGHANDLER_H
#define BUBBLEBURSTINGHANDLER_H

// Includes
#include <PSICluster.h>
#include <PSIClusterReactionNetwork.h>
#include <IBubbleBurstingHandler.h>
#include <memory>

namespace xolotlCore {

/**
 * This class realizes the IBubbleBursting interface, responsible for all the
 * physical parts for the bursting of HeV bubbles.
 */
class BubbleBurstingHandler: public IBubbleBurstingHandler {
protected:

	//! The bursting rate
	double kBursting;

	/**
	 * The vector containing the indices of the clusters taking part in the
	 * bursting process for each grid point
	 */
	std::vector<std::vector<int> > indexVector;


public:

	/**
	 * The constructor
	 */
	BubbleBurstingHandler() : kBursting(0.0){}

	/**
	 * The destructor
	 */
	~BubbleBurstingHandler() {}

	/**
	 * The initialize method has to add connectivity between the V clusters and HeV clusters
	 * of same number of V. It must also initialize the rates of the reactions and define
	 * which bubbles can burst at each grid point.
	 *
	 * @param network The network
	 * @param hx The grid step size
	 * @param nGrid The number of points on the grid
	 * @param surfacePos The index of the position on the surface
	 */
	virtual void initialize(PSIClusterReactionNetwork *network, double hx,
			int nGrid, int surfacePos);

	/**
	 * This method update the rate for the bubble bursting if the rates changed in the network,
	 * it should be called when temperature changes for instance.
	 *
	 * @param network The network
	 */
	virtual void updateBurstingRate(PSIClusterReactionNetwork *network);

	/**
	 * Compute the flux due to the bubble bursting for all the cluster,
	 * given the position index xi.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * A --> B
	 *
	 * F(B) = -F(A) = kBursting * C_A
	 *
	 * @param network The network
	 * @param xi The index of the position on the grid
	 * @param surfacePos The index of the position on the surface
	 * @param concOffset The pointer to the array of concentration at the grid
	 * point where the bursting is computed
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the bursting is computed used to find the next solution
	 */
	virtual void computeBursting(PSIClusterReactionNetwork *network,
			int xi, int surfacePos, double *concOffset, double *updatedConcOffset);

	/**
	 * Compute the partials due to the bubble bursting for all the clusters given
	 * the position index xi. Returns the number of bubbles that can possibly burst.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * A --> B
	 *
	 * dF(B)/dC_A = -dF(A)/dC_A = kBursting
	 *
	 * @param network The network
	 * @param val The pointer to the array that will contain the values of partials
	 * for the bursting
	 * @param indices The pointer to the array that will contain the indices of the clusters
	 * @param xi The index of the grip point
	 * @param surfacePos The index of the position on the surface
	 *
	 * @return The number of bubbles that can burst at this grid point
	 */
	virtual int computePartialsForBursting(PSIClusterReactionNetwork *network,
			double *val, int *indices, int xi, int surfacePos);

};
//end class BubbleBurstingHandler

} /* namespace xolotlCore */
#endif
