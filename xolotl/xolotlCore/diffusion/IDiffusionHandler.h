#ifndef IDIFFUSIONHANDLER_H
#define IDIFFUSIONHANDLER_H

// Includes
#include <petscsys.h>
#include <PSICluster.h>

namespace xolotlCore {

/**
 * Realizations of this interface are responsible for all the physical parts
 * for the diffusion of mobile cluster. The solver call these methods to handle
 * the diffusion.
 */
class IDiffusionHandler {

public:

	/**
	 * The destructor
	 */
	virtual ~IDiffusionHandler(){}

	/**
	 * Initialize the off-diagonal part of the Jacobian. If this step is skipped it
	 * won't be possible to set the partials for the diffusion.
	 *
	 * @param cluster The diffusing cluster
	 * @param size The size of the network
	 * @param ofill The pointer to the array that will contain the value 1 at the indices
	 * of the diffusing clusters, 0 if they are not diffusing
	 */
	virtual void initializeOFill(PSICluster * cluster, int size, PetscInt *ofill) = 0;

	/**
	 * Compute the flux due to the diffusion of the cluster given the space parameter sx.
	 * This method is called by the RHSFunction from the PetscSolver.
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
	virtual void computeDiffusion(PSICluster * cluster, double sx,
			PetscScalar *concOffset, PetscScalar *leftConcOffset,
			PetscScalar *rightConcOffset, PetscScalar *updatedConcOffset) = 0;

	/**
	 * Compute the partials due to the diffusion of the cluster given the space parameter sx.
	 * This method is called by the RHSJacobian from the PetscSolver.
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
	virtual void computePartialsForDiffusion(PSICluster * cluster,
			PetscReal sx, PetscReal val[3], PetscInt row[1], PetscInt col[3],
			PetscInt xi, PetscInt xs, int size) = 0;

}; //end class IDiffusionHandler

} /* namespace xolotlCore */
#endif
