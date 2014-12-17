#ifndef DIFFUSION2DHANDLER_H
#define DIFFUSION2DHANDLER_H

// Includes
#include "DiffusionHandler.h"

namespace xolotlCore {

/**
 * This class is a subclass of DiffusionHandler for the isotropic diffusion in 2D.
 */
class Diffusion2DHandler: public DiffusionHandler {
public:

	//! The Constructor
	Diffusion2DHandler() {}

	//! The Destructor
	~Diffusion2DHandler() {}

	/**
	 * Compute the flux due to the diffusion for all the cluster that are diffusing,
	 * given the space parameter s.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * If D is the diffusion coefficient, and C_l, C_r, C_m, C_t, C_b the left, right,
	 * middle, top, and bottom concentration of this cluster, the value to add to the
	 * updated concentration is:
	 *
	 * D * s * (C_l + C_r + C_t + C_b - 4 * C_m)
	 *
	 * @param network The network
	 * @param s The space parameter, depending on the grid step size
	 * @param concVector The pointer to the pointer of arrays of concentration at middle, left,
	 * right, top, and bottom grid points
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the diffusion is computed used to find the next solution
	 */
	void computeDiffusion(PSIClusterReactionNetwork *network,
			double s, double **concVector, double *updatedConcOffset);

	/**
	 * Compute the partials due to the diffusion of all the diffusing clusters given
	 * the space parameter s.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * On the directly neighboring grid points the partial will be D * s with D the diffusion
	 * coefficient, and on this grid point the value of the partial will be -4 * D * s.
	 *
	 * @param network The network
	 * @param s The space parameter, depending on the grid step size
	 * @param val The pointer to the array that will contain the values of partials
	 * for the diffusion
	 * @param indices The pointer to the array that will contain the indices of the
	 * diffusing clusters in the network
	 */
	void computePartialsForDiffusion(PSIClusterReactionNetwork *network,
			double s, double *val, int *indices);

};
//end class Diffusion2DHandler

} /* end namespace xolotlCore */
#endif
