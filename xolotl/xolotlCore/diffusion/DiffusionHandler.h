#ifndef DIFFUSIONHANDLER_H
#define DIFFUSIONHANDLER_H

// Includes
#include "IDiffusionHandler.h"

namespace xolotlCore {

/**
 * This class realizes the IDiffusionHandler interface responsible for all
 * the physical parts for the diffusion of mobile cluster.
 */
class DiffusionHandler: public IDiffusionHandler {
public:

	//! The Constructor
	DiffusionHandler() {}

	//! The Destructor
	~DiffusionHandler() {}

	/**
	 * Initialize the off-diagonal part of the Jacobian. If this step is skipped it
	 * won't be possible to set the partials for the diffusion.
	 *
	 * @param cluster The diffusing cluster
	 * @param size The size of the network
	 * @param ofill The pointer to the array that will contain the value 1 at the indices
	 * of the diffusing clusters, 0 if they are not diffusing
	 */
	void initializeOFill(PSICluster * cluster, int size, PetscInt *ofill);

	/**
	 * Compute the change of concentration due to the diffusion of the cluster
	 * given the space parameter sx. This method is called by the RHSFunction
	 * from the PetscSolver.
	 *
	 * If D is the diffusion coefficient, and C_l, C_r, C_m the left, right,
	 * middle concentration of this cluster, the value to add to the updated
	 * concentration is:
	 *
	 * D * sx * (C_l + C_r - 2 * C_m)
	 *
	 * @param cluster The diffusing cluster
	 * @param sx The space parameter, depending on the grid step size
	 * @param concOffset The pointer to the array of concentration at the grid
	 * point where the diffusion is computed
	 * @param leftConcOffset The pointer to the array of concentration at the grid
	 * point to the left of where the diffusion is computed
	 * @param rightConcOffset The pointer to the array of concentration at the grid
	 * point to the right of where the diffusion is computed
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the diffusion is computed used to find the next solution
	 */
	void computeDiffusion(PSICluster * cluster, double sx,
			PetscScalar *concOffset, PetscScalar *leftConcOffset,
			PetscScalar *rightConcOffset, PetscScalar *updatedConcOffset);

	/**
	 * Compute the partials due to the diffusion of the cluster given the space parameter sx.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * On the left and right grid point the partial will be D * sx with D the diffusion
	 * coefficient, and on this grid point the value of the partial will be -2 * D * sx.
	 *
	 * @param cluster The diffusing cluster
	 * @param sx The space parameter, depending on the grid step size
	 * @param val The array that will contain the values of partials for the diffusion
	 * @param row The array that will contain the indices of the row for the Jacobian
	 * @param col The array that will contain the indices of the columns for the Jacobian
	 * @param xi The index of the grip point
	 * @param xs The index of the first grid point on the locally owned grid
	 * @param size The size of the network
	 */
	void computePartialsForDiffusion(PSICluster * cluster, PetscReal sx,
			PetscReal val[3], PetscInt row[1], PetscInt col[3], PetscInt xi,
			PetscInt xs, int size);

};
//end class DiffusionHandler

} /* end namespace xolotlCore */
#endif
