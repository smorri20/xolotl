#ifndef DIFFUSION3DHANDLER_H
#define DIFFUSION3DHANDLER_H

// Includes
#include "DiffusionHandler.h"

namespace xolotlCore {

/**
 * This class is a subclass of DiffusionHandler for the isotropic diffusion in 3D.
 */
class Diffusion3DHandler: public DiffusionHandler {
public:

	//! The Constructor
	Diffusion3DHandler() {}

	//! The Destructor
	~Diffusion3DHandler() {}

	/**
	 * Compute the flux due to the diffusion for all the cluster that are diffusing,
	 * given the space parameter s.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * If D is the diffusion coefficient, and C_l, C_r, C_m, C_t, C_b, C_f, C_ba the
	 * left, right, middle, top, bottom, front, and back concentration of this cluster,
	 * the value to add to the updated concentration is:
	 *
	 * D * s * (C_l + C_r + C_t + C_b + C_f + C_ba - 6 * C_m)
	 *
	 * @param network The network
	 * @param s The space parameter, depending on the grid step size
	 * @param concVector The pointer to the pointer of arrays of concentration at left,
	 * middle, right grid points
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
	 * coefficient, and on this grid point the value of the partial will be -6 * D * s.
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
//end class Diffusion3DHandler

} /* end namespace xolotlCore */
#endif
