#ifndef ITRAPMUTATIONHANDLER_H
#define ITRAPMUTATIONHANDLER_H

// Includes
#include <PSICluster.h>
#include <PSIClusterReactionNetwork.h>
#include <memory>

namespace xolotlCore {

/**
 * Realizations of this interface are responsible for all the physical parts
 * for the modified trap-mutation of small helium clusters close to the surface.
 * The solver call these methods to handle the modified trap-mutation.
 */
class ITrapMutationHandler {

public:

	/**
	 * The destructor
	 */
	virtual ~ITrapMutationHandler() {}

	/**
	 * The initialize method has to add connectivity between the He clusters and
	 * HeV clusters of same number of He, and I. It must also initialize the
	 * rates of the reactions and define which trap-mutation is allowed at
	 * each grid point.
	 *
	 * @param network The network
	 * @param grid The grid on the x axis
	 */
	virtual void initialize(PSIClusterReactionNetwork *network,
			std::vector<double> grid) = 0;

	/**
	 * This method update the rate for the modified trap-mutation if the rates
	 * changed in the network, it should be called when temperature changes
	 * for instance.
	 *
	 * @param network The network
	 */
	virtual void updateTrapMutationRate(PSIClusterReactionNetwork *network) = 0;

	/**
	 * Compute the flux due to the modified trap-mutation for all the cluster,
	 * given the position index xi.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * @param network The network
	 * @param xi The index of the position on the grid
	 * @param concOffset The pointer to the array of concentration at the grid
	 * point where the trap-mutation is computed
	 * @param updatedConcOffset The pointer to the array of the concentration
	 * at the grid point where the trap-mutation is computed used to find the
	 * next solution
	 */
	virtual void computeTrapMutation(PSIClusterReactionNetwork *network,
			int xi, double *concOffset, double *updatedConcOffset) = 0;

	/**
	 * Compute the partials due to the modified trap-mutation for all the
	 * clusters given the position index xi. Returns the number of helium
	 * clusters that are undergoing trap-mutation at this grid point.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * @param network The network
	 * @param val The pointer to the array that will contain the values of
	 * partials for the trap-mutation
	 * @param indices The pointer to the array that will contain the indices
	 * of the clusters
	 * @param xi The index of the grip point
	 *
	 * @return The number of helium clusters that go through modified trap-mutation
	 * at this grid point
	 */
	virtual int computePartialsForTrapMutation(PSIClusterReactionNetwork *network,
			double *val, int *indices, int xi) = 0;

};
//end class ITrapMutationHandler

} /* namespace xolotlCore */
#endif
