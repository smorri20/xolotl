#ifndef DIFFUSION1DHANDLER_H
#define DIFFUSION1DHANDLER_H

// Includes
#include "DiffusionHandler.h"

namespace xolotlCore {

/**
 * This class is a subclass of DiffusionHandler for the isotropic diffusion in 1D.
 */
class Diffusion1DHandler: public DiffusionHandler {
public:

	//! The Constructor
	Diffusion1DHandler() {}

	//! The Destructor
	~Diffusion1DHandler() {}

	/**
	 * Compute the flux due to the diffusion for all the cluster that are diffusing,
	 * given the different space parameters.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * If D is the diffusion coefficient, C_l, C_r, C_m the left, right, and
	 * middle concentration of this cluster, and a and b are the step sizes on the
	 * left and right side on the grid point, the value to add to the updated
	 * concentration is:
	 *
	 * D * (2.0 / [a * (a + b)]) * (C_l + [a/b] * C_r - [1.0 + (a/b)] * C_m)
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
			double hxLeft, double hxRight, double sy = 0.0, double sz = 0.0);

	/**
	 * Compute the partials due to the diffusion of all the diffusing clusters given the
	 * different space parameters.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * Using the same notation as for computeDiffusion, the partial derivative on the
	 * left grid point should be:
	 *
	 * D * (2.0 / [a * (a + b)])
	 *
	 * on the right grid point:
	 *
	 * D * (2.0 / [b * (a + b)])
	 *
	 * and on the middle grid point:
	 *
	 * - D * (2.0 / [a*b])
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
			double sy = 0.0, double sz = 0.0);

};
//end class Diffusion1DHandler

} /* end namespace xolotlCore */
#endif
