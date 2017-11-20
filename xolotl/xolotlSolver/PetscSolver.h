#ifndef PETSCSOLVER_H
#define PETSCSOLVER_H

// Includes
#include "Solver.h"

namespace xolotlSolver {

#ifndef CHECK_PETSC_ERROR
#define CHECK_PETSC_ERROR
/**
 * This operation checks a PETSc error code and throws an exception with given error message.
 *
 * @param errorCode The PETSc error code.
 * @param errMsg The error message in the thrown exception.
 */
inline void checkPetscError(PetscErrorCode errorCode, const char* errorMsg) {
	if (PetscUnlikely(errorCode))
		throw std::string(errorMsg);
}
#endif

/**
 * This class realizes the ISolver interface to solve the
 * advection-diffusion-reaction problem with the PETSc solvers from Argonne
 * National Laboratory.
 */
class PetscSolver: public Solver {
private:

	/**
	 * This operation configures the initial conditions of the grid in Xolotl.
	 *
	 * @param data The DM (data manager) created by PETSc
	 * @param solutionVector The solution vector that contains the PDE
	 * solution and which needs to be initialized.
	 * @param oldC The solution from the previous loop
	 * @param idMap The map of ids from the previous loop
	 */
	void setupInitialConditions(DM &data, Vec &solutionVector, std::vector<double> &oldC,
			std::map<std::string, int> &idMap);

	/**
	 * This operation configures the bounds to use for the grouping scheme
	 *
	 * @param peak The peak position
	 * @param min The minimum value
	 * @param max The maximum value
	 * @return The vector of bounds
	 */
	std::vector<IReactant::SizeType> generateBounds(int peak, int min, int max);

	/**
	 * This operation creates concentration vector for the next loop
	 * from the current network state.
	 *
	 * @param network The current network
	 * @param bounds1 The bounds on the next loop's network
	 * @param bounds2 The bounds on the next loop's network
	 * @param idMap The map of ids for the concentrations
	 * @param concVec The new vector of concentrations
	 */
	void transformConcentrationVector(IReactionNetwork & network,
			const std::vector<IReactant::SizeType> & bounds1,
			const std::vector<IReactant::SizeType> & bounds2,
			std::map<std::string, int> &idMap,
			std::vector<double> & concVec);

public:

	/**
	 * Default constructor, deleted because we must construct using arguments.
	 */
	PetscSolver() = delete;

	//! The Constructor
	PetscSolver(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	//! The Destructor
	~PetscSolver();

	/**
	 * This operation sets up the mesh that will be used by the solver and
	 * initializes the data on that mesh. This operation will throw an exception
	 * of type std::string if the mesh can not be setup.
	 */
	void setupMesh() override;

	/**
	 * This operation performs all necessary initialization for the solver
	 * possibly including but not limited to setting up MPI and loading initial
	 * conditions. If the solver can not be initialized, this operation will
	 * throw an exception of type std::string.
	 */
	void initialize() override;

	/**
	 * This operation directs the Solver to perform the solve. If the solve
	 * fails, it will throw an exception of type std::string.
	 */
	void solve(Options &options) override;

	/**
	 * This operation performs all necessary finalization for the solver
	 * including but not limited to cleaning up memory, finalizing MPI and
	 * printing diagnostic information. If the solver can not be finalized,
	 * this operation will throw an exception of type std::string.
	 */
	void finalize() override;

};
//end class PetscSolver

} /* end namespace xolotlSolver */

// Some compilers (e.g., recent versions of Intel) define __func__ 
// to include the namespace or class scope when compiled with the C++11
// support enabled.  Others don't.  Because PETSc's PetscFunctionBeginUser
// does a straight string comparison between what we call the function name
// and what it determines from the compiler, we need a way to provide
// either the scoped name or the unscoped name.
#if defined(__ICC) || defined(__INTEL_COMPILER)
#  define Actual__FUNCT__(sname,fname)  sname "::" fname
#else
#  define Actual__FUNCT__(sname,fname)  fname
#endif /* if it is the Intel compiler */

#endif
