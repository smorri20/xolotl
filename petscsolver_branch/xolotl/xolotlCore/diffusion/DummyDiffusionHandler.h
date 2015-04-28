#ifndef DUMMYDIFFUSIONHANDLER_H
#define DUMMYDIFFUSIONHANDLER_H

// Includes
#include "DiffusionHandler.h"

namespace xolotlCore {

/**
 * This class realizes the IDiffusionHandler interface responsible for all
 * the physical parts for the diffusion of mobile cluster. Here it is a dummy class,
 * meaning that it should not do anything.
 */
class DummyDiffusionHandler: public DiffusionHandler {
public:

	//! The Constructor
	DummyDiffusionHandler() {}

	//! The Destructor
	~DummyDiffusionHandler() {}

	/**
	 * Initialize the off-diagonal part of the Jacobian. If this step is skipped it
	 * won't be possible to set the partials for the diffusion.
	 *
	 * We don't want any cluster to diffuse, so nothing is set to 1 in ofill, and no index
	 * is added to the vector.
	 *
	 * @param network The network
	 * @param ofill The pointer to the array that will contain the value 1 at the indices
	 * of the diffusing clusters
	 */
	void initializeOFill(PSIClusterReactionNetwork *network, int *ofill) {
		// Clear the index vector
		indexVector.clear();

		// And don't do anything else
		return;
	}

	/**
	 * Compute the flux due to the diffusion for all the cluster that are diffusing,
	 * given the space parameters.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * Here it won't do anything because it is a dummy class.
	 *
	 * @param network The network
	 * @param concVector The pointer to the pointer of arrays of concentration at middle,
	 * left, and right grid points
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the diffusion is computed used to find the next solution
	 * @param hxLeft The step size on the left side of the point in the x direction (a)
	 * @param hxRight The step size on the right side of the point in the x direction (b)
	 * @param sy The space parameter, depending on the grid step size in the y direction
	 * @param sz The space parameter, depending on the grid step size in the z direction
	 */
	void computeDiffusion(PSIClusterReactionNetwork *network,
			double **concVector, double *updatedConcOffset,
			double hxLeft, double hxRight, double sy = 0.0, double sz = 0.0) {return;}

	/**
	 * Compute the partials due to the diffusion of all the diffusing clusters given
	 * the space parameters.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * Here it won't do anything because it is a dummy class.
	 *
	 * @param network The network
	 * @param val The pointer to the array that will contain the values of partials
	 * for the diffusion
	 * @param indices The pointer to the array that will contain the indices of the
	 * diffusing clusters in the network
	 * @param hxLeft The step size on the left side of the point in the x direction (a)
	 * @param hxRight The step size on the right side of the point in the x direction (b)
	 * @param sy The space parameter, depending on the grid step size in the y direction
	 * @param sz The space parameter, depending on the grid step size in the z direction
	 */
	void computePartialsForDiffusion(PSIClusterReactionNetwork *network,
			double *val, int *indices, double hxLeft, double hxRight,
			double sy = 0.0, double sz = 0.0) {return;}

};
//end class DummyDiffusionHandler

} /* end namespace xolotlCore */
#endif
