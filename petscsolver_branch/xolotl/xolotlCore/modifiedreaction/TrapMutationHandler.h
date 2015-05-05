#ifndef TRAPMUTATIONHANDLER_H
#define TRAPMUTATIONHANDLER_H

// Includes
#include <ITrapMutationHandler.h>

namespace xolotlCore {

/**
 * This class realizes ITrapMutationHandler interface, responsible for the modified
 * trap-mutation of small helium clusters close to the surface.
 */
class TrapMutationHandler: public ITrapMutationHandler {
protected:

	//! The vector containing the different depths for the modified trap mutation
	std::vector<double> depthVec;

	/**
	 * The vector containing the size of the smallest cluster undergoing
	 * trap-mutation for each given depth
	 */
	std::vector<int> sizeVec;

	//! The trap-mutation rate
	double kMutation;

	/**
	 * The vector containing the indices of the clusters undergoing modified
	 * trap-mutation for each grid point
	 */
	std::vector<std::vector<int> > indexVector;

	/**
	 * Method initializing the depth and size vectors.
	 * It needs to be implemented by the daughter classes.
	 */
	virtual void initializeDepthSize() {return;}

public:

	/**
	 * The constructor
	 */
	TrapMutationHandler() : kMutation(0.0){}

	/**
	 * The destructor
	 */
	~TrapMutationHandler() {}

	/**
	 * The initialize method has to add connectivity between the He clusters and
	 * HeV clusters of same number of He, and I. It must also initialize the
	 * rates of the reactions and define which trap-mutation is allowed at
	 * each grid point.
	 *
	 * @param surfacePos The index of the position of the surface
	 * @param network The network
	 * @param grid The grid on the x axis
	 */
	void initialize(int surfacePos, PSIClusterReactionNetwork *network,
			std::vector<double> grid);

	/**
	 * This method defines which trap-mutation is allowed at each grid point.
	 *
	 * @param surfacePos The index of the position of the surface
	 * @param network The network
	 * @param grid The grid on the x axis
	 */
	void initializeIndex(int surfacePos, PSIClusterReactionNetwork *network,
			std::vector<double> grid);

	/**
	 * This method update the rate for the modified trap-mutation if the rates
	 * changed in the network, it should be called when temperature changes
	 * for instance.
	 *
	 * @param network The network
	 */
	void updateTrapMutationRate(PSIClusterReactionNetwork *network);

	/**
	 * Compute the flux due to the modified trap-mutation for all the cluster,
	 * given the position index xi.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * He_i --> (He_i)(V) + I
	 *
	 * F(He_i) = -F[(He_i)(V)] = -F(I) = -kMutation * C_(He_i)
	 *
	 * @param network The network
	 * @param xi The index of the position on the grid
	 * @param concOffset The pointer to the array of concentration at the grid
	 * point where the trap-mutation is computed
	 * @param updatedConcOffset The pointer to the array of the concentration
	 * at the grid point where the trap-mutation is computed used to find the
	 * next solution
	 */
	void computeTrapMutation(PSIClusterReactionNetwork *network,
			int xi, double *concOffset, double *updatedConcOffset);

	/**
	 * Compute the partials due to the modified trap-mutation for all the
	 * clusters given the position index xi. Returns the number of helium
	 * clusters that are undergoing trap-mutation at this grid point.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * He_i --> (He_i)(V) + I
	 *
	 * dF(He_i)/dC_(He_i) = -dF[(He_i)(V)]/dC_(He_i) = -dF(I)/dC_(He_i)
	 * 		= -kMutation
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
	int computePartialsForTrapMutation(PSIClusterReactionNetwork *network,
			double *val, int *indices, int xi);

};
//end class TrapMutationHandler

} /* namespace xolotlCore */
#endif
