#ifndef PETSCSOLVERHANDLER_H
#define PETSCSOLVERHANDLER_H

// Includes
#include "SolverHandler.h"

namespace xolotlSolver {

/**
 * The last temperature on the grid. In the future this will have to be an
 * array or map, but for now the temperature is isotropic.
 */
extern double lastTemperature;

/**
 * A map for storing the dfill configuration and accelerating the formation of
 * the Jacobian. Its keys are reactant/cluster ids and its values are integer
 * vectors of the column ids that are marked as connected for that cluster in
 * the dfill array.
 */
extern std::unordered_map<int, std::vector<int> > dFillMap;

/**
 * A pointer to all of the reactants in the network. It is retrieved from the
 * network after it is set.
 */
extern std::shared_ptr<std::vector<xolotlCore::Reactant *>> allReactants;

/**
 * A vector for holding the partial derivatives of one cluster. It is sized in
 * the createSolverContext() operation.
 *
 * The vector is used for every cluster and immediately reset to zero before
 * being used for the next. This allows the acquisition of the partial
 * derivatives to take up minimal memory and require no additional dynamic
 * allocations.
 */
extern std::vector<double> clusterPartials;

/**
 * A vector for holding the partial derivatives for one cluster in the order
 * that PETSc expects. It is sized in the createSolverContext() operation.
 *
 * The vector is used for every cluster and immediately reset to zero before
 * being used for the next. This allows the acquisition of the partial
 * derivatives to take up minimal memory and require no additional dynamic
 * allocations.
 */
extern std::vector<double> reactingPartialsForCluster;

/**
 * This class and its subclasses realize the ISolverHandler interface to solve the
 * advection-diffusion-reaction problem with the PETSc solvers from Argonne
 * National Laboratory.
 *
 * This class does NOT implement most of the methods that are needed by the
 * PetscSolver. Only subclasses of this class must be used by the PetscSolver.
 */
class PetscSolverHandler: public SolverHandler {
public:

	//! The Constructor
	PetscSolverHandler() {}

	//! The Destructor
	~PetscSolverHandler() {}

	/**
	 * Get the diagonal fill for the Jacobian, corresponding to the reactions.
	 * \see ISolverHandler.h
	 */
	void getDiagonalFill(PetscInt *diagFill, int diagFillSize) const;

}; //end class PetscSolverHandler

} /* end namespace xolotlSolver */
#endif
