// Includes
#include <cassert>
#include <PetscSolver.h>
#include <HDF5Utils.h>
#include <MathUtils.h>
#include <Constants.h>
#include <SolverHandlerFactory.h>
#include <IReactionHandlerFactory.h>
#include <TemperatureHandlerFactory.h>

using namespace xolotlCore;

/*
 C_t =  -D*C_xx + A*C_x + F(C) + R(C) + D(C) from Brian Wirth's SciDAC project.

 D*C_xx  - diffusion of He and V and I
 A*C_x   - advection of He
 F(C)    - forcing function; He being created.
 R(C)    - reaction terms   (clusters combining)
 D(C)    - dissociation terms (cluster breaking up)

 Sample Options:
 -da_grid_x <nx>						 -- number of grid points in the x direction
 -ts_max_steps <maxsteps>                -- maximum number of time-steps to take
 -ts_final_time <time>                   -- maximum time to compute to
 -ts_dt <size>							 -- initial size of the time step

 */

namespace xolotlSolver {

//Timer for RHSFunction()
std::shared_ptr<xolotlPerf::ITimer> RHSFunctionTimer;

////Timer for RHSJacobian()
std::shared_ptr<xolotlPerf::ITimer> RHSJacobianTimer;

//! Help message
static char help[] =
		"Solves C_t =  -D*C_xx + A*C_x + F(C) + R(C) + D(C) from Brian Wirth's SciDAC project.\n";

// ----- GLOBAL VARIABLES ----- //
extern PetscErrorCode setupPetsc0DMonitor(TS);
extern PetscErrorCode setupPetsc1DMonitor(TS);
extern PetscErrorCode setupPetsc2DMonitor(TS);
extern PetscErrorCode setupPetsc3DMonitor(TS);

void PetscSolver::setupInitialConditions(DM &da, Vec &C,
		std::vector<double> & oldC, std::map<std::string, int> &idMap) {
	// Initialize the concentrations in the solution vector
	auto& solverHandler = Solver::getSolverHandler();
	solverHandler.initializeConcentration(da, C, oldC, idMap);

	return;
}

std::vector<IReactant::SizeType> PetscSolver::generateBounds(int peak, int min,
		int max) {
	// Initial declaration
	int value = min;
	double scale = max / 25.0, tau = (double) (max - peak) / 7.0;
	std::vector<IReactant::SizeType> bounds;
	bounds.push_back(1);

	// Generate the bounds
	while (value < max) {
//		std::cout << value << std::endl;

		bounds.push_back(value);

//		// Update the value
//		int x = 0;
//		if (value < peak) x = scale * (1.0 - exp((double)(value - peak)/tau));
//		else x = scale * (1.0 - exp(-(double)(value - peak)/tau));
		value += scale;
	}

	// Add the last value
	bounds.push_back(max + 1);

	return bounds;
}

void PetscSolver::transformConcentrationVector(IReactionNetwork & network,
		const std::vector<IReactant::SizeType> & bounds1,
		const std::vector<IReactant::SizeType> & bounds2,
		std::map<std::string, int> &idMap, std::vector<double> & concVec,
		std::vector<std::vector<double> > & padeVec) {
	// Initial declarations
	int i = 0;
	std::vector<double> coordinates = { -1.0, -0.5, 0.0, 0.5, 1.0 };
	// Reset the concentration vector and id map
	concVec.clear();
	padeVec.clear();
	idMap.clear();

	// Take care of the normal clusters
	// Make a vector of types for the non super clusters
	std::vector<ReactantType> typeVec { ReactantType::He, ReactantType::V,
			ReactantType::I, ReactantType::HeV };
	// Loop on it
	for (auto tvIter = typeVec.begin(); tvIter != typeVec.end(); ++tvIter) {
		// Consider all reactants of the current type.
		auto const& currTypeReactantMap = network.getAll(*tvIter);

		// Update the column in the Jacobian that represents each normal reactant
		for (auto const& currMapItem : currTypeReactantMap) {
			// Save the current cluster concentration
			concVec.push_back(currMapItem.second->getConcentration());
			// Save an empty vector for pade
			padeVec.push_back(std::vector<double>());
			// Save it in the map
			idMap[currMapItem.second->getName()] = i;
			i++;
		}
	}

	// Loop on the bounds
	for (int k = 0; k < bounds2.size() - 1; k++) {
		for (int j = 0; j < bounds1.size() - 1; j++) {
			// Non super cluster case
			if (j == 0 && k == 0)
				continue;
			// Super cluster case
			else {
				// Initialize data for the loop
				double mean1 = (double) (bounds1[j] + bounds1[j + 1] - 1) / 2.0;
				double mean2 = (double) (bounds2[k] + bounds2[k + 1] - 1) / 2.0;
				double conc = 0.0, mom1 = 0.0, mom2 = 0.0;

				// Loop on all the values
				for (int numV = bounds2[k]; numV < bounds2[k + 1]; numV++) {
					for (int numHe = bounds1[j]; numHe < bounds1[j + 1];
							numHe++) {
						// Get the corresponding super cluster
						auto cluster = network.getSuperFromComp(numHe, numV);
						// If it exists
						if (cluster) {
							// Compute the distances
							double heDistance = cluster->getHeDistance(numHe);
							double vDistance = cluster->getVDistance(numV);
							// Get its concentration
							double concentration = cluster->getConcentration(
									heDistance, vDistance);
							// Add it to the running quantities
							conc += concentration;
							if (!xolotlCore::equal(numHe, mean1))
								mom1 += concentration
										/ ((double) numHe - mean1);
							if (!xolotlCore::equal(numV, mean2))
								mom2 += concentration / ((double) numV - mean2);
						}
					}
				}

				// Rescale the running quantities
				double nTot = (double) ((bounds1[j + 1] - bounds1[j])
						* (bounds2[k + 1] - bounds2[k]));
				conc = conc / nTot;
				mom1 = (double) (bounds1[j + 1] - bounds1[j] - 1) * mom1
						/ (nTot * 2.0);
				mom2 = (double) (bounds2[k + 1] - bounds2[k] - 1) * mom2
						/ (nTot * 2.0);

				// Save them
				concVec.push_back(conc);
				std::stringstream nameStream;
				nameStream << "He_" << mean1 << "V_" << mean2;
				std::string clusterName = nameStream.str();
				idMap[clusterName] = i;
				i++;
				concVec.push_back(mom1);
				std::string name = clusterName + "He";
				idMap[name] = i;
				i++;
				concVec.push_back(mom2);
				name = clusterName + "V";
				idMap[name] = i;
				i++;

				// Check if we need the pade approximation
				if ((bounds1[j + 1] - bounds1[j]) > 4
						&& (bounds2[k + 1] - bounds2[k]) > 4) {

//					std::cout << bounds1[j] << " " << bounds2[k] << std::endl;

					// Compute the approximation for 24 data points
					PetscInt dataSize = 24, col[24];
					PetscScalar value[24], rhs[24];
					for (int n = 0; n < dataSize; n++) {
						col[n] = n;
					}

					// Create a PETSc vector for the solution and other side
					// Ax = b
					Vec x, b;
					VecCreate(PETSC_COMM_WORLD, &x);
					PetscObjectSetName((PetscObject) x, "Solution");
					VecSetSizes(x, PETSC_DECIDE, dataSize);
					VecSetFromOptions(x);
					VecSet(x, 0.0);
					VecDuplicate(x, &b);
					VecSet(b, 0.0);

					// Create a matrix for the equations to solve
					Mat A;
					MatCreate(PETSC_COMM_WORLD, &A);
					MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, dataSize,
							dataSize);
					MatSetFromOptions(A);
					MatSetUp(A);

					// Generate a random number to skip one of the grid points
					int toSkip = (int) (((double) rand() / (double) RAND_MAX)
							* 25.0);
					bool alreadySkipped = false, dontPade = false;
					// Loop on the coordinates
					int count = 0;
					for (int b = 0; b < coordinates.size(); b++) {
						// Compute the associated V number
						int numV = (int) (mean2
								+ coordinates[b]
										* (double) (bounds2[k + 1] - bounds2[k])
										/ 2.0);

						for (int a = 0; a < coordinates.size(); a++) {
							// Compute the associated V number
							int numHe = (int) (mean1
									+ coordinates[a]
											* (double) (bounds1[j + 1]
													- bounds1[j]) / 2.0);

							// Check if we should skip this point
							if (count == toSkip && !alreadySkipped) {
								alreadySkipped = true;
								continue;
							}

							// Get the concentration at this point
							// Get the corresponding super cluster
							auto cluster = network.getSuperFromComp(numHe,
									numV);
							// If it exists
							if (cluster) {
								// Compute the distances in the old cluster
								double heDistance = cluster->getHeDistance(
										numHe);
								double vDistance = cluster->getVDistance(numV);
								// Get its concentration
								double concentration =
										cluster->getConcentration(heDistance,
												vDistance);

								// Compute the distances in the new cluster
								heDistance = 2.0 * ((double) numHe - mean1)
										/ (double) (bounds1[j + 1] - bounds1[j]
												- 1);
								vDistance = 2.0 * ((double) numV - mean2)
										/ (double) (bounds2[k + 1] - bounds2[k]
												- 1);

								// Get the difference between the concentration and the moment approximation
								double model = conc + (heDistance * mom1)
										+ (vDistance * mom2);
								concentration = concentration - model;

								// Set all this in the solver
								value[0] = 1.0;
								value[1] = heDistance;
								value[2] = vDistance;
								value[3] = vDistance * vDistance;
								value[4] = heDistance * vDistance;
								value[5] = heDistance * heDistance;
								value[6] = heDistance * heDistance * heDistance;
								value[7] = vDistance * vDistance * vDistance;
								value[8] = heDistance * heDistance * vDistance;
								value[9] = heDistance * vDistance * vDistance;
								value[10] = -concentration * heDistance;
								value[11] = -concentration * vDistance;
								value[12] = -concentration * heDistance
										* heDistance;
								value[13] = -concentration * vDistance
										* vDistance;
								value[14] = -concentration * heDistance
										* vDistance;
								value[15] = -concentration * heDistance
										* heDistance * heDistance;
								value[16] = -concentration * vDistance
										* vDistance * vDistance;
								value[17] = -concentration * heDistance
										* heDistance * vDistance;
								value[18] = -concentration * heDistance
										* vDistance * vDistance;
								value[19] = -concentration * heDistance
										* heDistance * heDistance * heDistance;
								value[20] = -concentration * vDistance
										* vDistance * vDistance * vDistance;
								value[21] = -concentration * heDistance
										* heDistance * heDistance * vDistance;
								value[22] = -concentration * heDistance
										* heDistance * vDistance * vDistance;
								value[23] = -concentration * heDistance
										* vDistance * vDistance * vDistance;
								rhs[count] = concentration;
								MatSetValues(A, 1, &count, dataSize, col, value,
										INSERT_VALUES);

//								std::cout << mean1 << " " << mean2 << " "
//										<< heDistance << " " << vDistance << " "
//										<< concentration << std::endl;

								count++;
							}

							else {
								dontPade = true;
								break;
							}
						}
					}

					if (dontPade) {
						// Save an empty vector for pade 3 times to keep up with the indices
						padeVec.push_back(std::vector<double>());
						padeVec.push_back(std::vector<double>());
						padeVec.push_back(std::vector<double>());

						continue;
					}

					// Do PETSc stuff
					MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
					MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
					VecSetValues(b, dataSize, col, rhs, INSERT_VALUES);
					VecAssemblyBegin(b);
					VecAssemblyEnd(b);

//					MatView(A, PETSC_VIEWER_STDOUT_SELF);
//					VecView(b, PETSC_VIEWER_STDOUT_SELF);

					// Compute the norm of the vector
					double val = 0.0;
					VecNorm(b, NORM_1, &val);

//					// Check if the norm is large enough
//					if (val < 1.0e-20) {
//						// Save an empty vector for pade 3 times to keep up with the indices
//						padeVec.push_back(std::vector<double>());
//						padeVec.push_back(std::vector<double>());
//						padeVec.push_back(std::vector<double>());
//
//						continue;
//					}

//					std::cout << mean1 << " " << mean2 << " " << val
//							<< std::endl;

					// Solve the equation
					KSP ksp;
					PC pc;
					KSPCreate(PETSC_COMM_WORLD, &ksp);
					KSPSetOperators(ksp, A, A);
					KSPGetPC(ksp, &pc);
					PCSetType(pc, PCJACOBI);
					KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT,
					PETSC_DEFAULT,
					PETSC_DEFAULT);
					KSPSetUp(ksp);
					KSPSolve(ksp, b, x);

					// Put the results in a vector
					std::vector<double> result;
					double *r = nullptr;
					VecGetArray(x, &r);
					for (int n = 0; n < dataSize; n++) {
						result.push_back(r[n]);
//						std::cout << r[n] << ", ";
					}
//					std::cout << std::endl;

					VecRestoreArray(x, &r);
					// Save the results to give to the clusters
					padeVec.push_back(result);

					// Save an empty vector for pade 2 times to keep up with the indices for moments
					padeVec.push_back(std::vector<double>());
					padeVec.push_back(std::vector<double>());

					double error = 0.0;
					// Loop on all the values
					for (int numV = bounds2[k]; numV < bounds2[k + 1]; numV++) {
						for (int numHe = bounds1[j]; numHe < bounds1[j + 1];
								numHe++) {
							// Get the corresponding super cluster
							auto cluster = network.getSuperFromComp(numHe,
									numV);
							// If it exists
							if (cluster) {
								// Compute the distances
								double heDistance = cluster->getHeDistance(
										numHe);
								double vDistance = cluster->getVDistance(numV);
								// Get its concentration
								double concentration =
										cluster->getConcentration(heDistance,
												vDistance);
								// Get the current distances
								heDistance = 2.0 * ((double) numHe - mean1)
										/ (double) (bounds1[j + 1] - bounds1[j]
												- 1);
								vDistance = 2.0 * ((double) numV - mean2)
										/ (double) (bounds2[k + 1] - bounds2[k]
												- 1);

								double numerator = result[0]
										+ result[1] * heDistance
										+ result[2] * vDistance
										+ result[3] * heDistance * heDistance
										+ result[4] * vDistance * vDistance
										+ result[5] * heDistance * vDistance
										+ result[6] * heDistance * heDistance
												* heDistance
										+ result[7] * vDistance * vDistance
												* vDistance
										+ result[8] * heDistance * heDistance
												* vDistance
										+ result[9] * heDistance * vDistance
												* vDistance;

								double denominator = 1.0
										+ result[10] * heDistance
										+ result[11] * vDistance
										+ result[12] * heDistance * heDistance
										+ result[13] * vDistance * vDistance
										+ result[14] * heDistance * vDistance
										+ result[15] * heDistance * heDistance
												* heDistance
										+ result[16] * vDistance * vDistance
												* vDistance
										+ result[17] * heDistance * heDistance
												* vDistance
										+ result[18] * heDistance * vDistance
												* vDistance
										+ result[19] * heDistance * heDistance
												* heDistance * heDistance
										+ result[20] * vDistance * vDistance
												* vDistance * vDistance
										+ result[21] * heDistance * heDistance
												* heDistance * vDistance
										+ result[22] * heDistance * heDistance
												* vDistance * vDistance
										+ result[23] * heDistance * vDistance
												* vDistance * vDistance;

								double test = fabs(
										(concentration
												- (numerator / denominator))
												/ concentration);
								if (test > error) error = test;

//								std::cout << numHe << " " << numV << " "
//										<< (concentration
//												- (numerator / denominator))
//										<< std::endl;
							}
						}
					}

					if (error > 500.0) {
						padeVec[padeVec.size() - 3] = std::vector<double>();
					}
				} else {
					// Save an empty vector for pade 3 times to keep up with the indices
					padeVec.push_back(std::vector<double>());
					padeVec.push_back(std::vector<double>());
					padeVec.push_back(std::vector<double>());
				}
			}
		}
	}

	return;
}

/* ------------------------------------------------------------------- */

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "RHSFunction")
/*
 RHSFunction - Evaluates the right-hand-side of the nonlinear function defining the ODE

 Input Parameters:
 .  ts - the TS context
 .  ftime - the physical time at which the function is evaluated
 .  C - input vector
 .  ptr - optional user-defined context

 Output Parameter:
 .  F - function values
 */
/* ------------------------------------------------------------------- */
PetscErrorCode RHSFunction(TS ts, PetscReal ftime, Vec C, Vec F, void *) {
	// Start the RHSFunction Timer
	RHSFunctionTimer->start();

	PetscErrorCode ierr;

	// Get the local data vector from PETSc
	PetscFunctionBeginUser;
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);
	Vec localC;
	ierr = DMGetLocalVector(da, &localC);
	CHKERRQ(ierr);

	// Scatter ghost points to local vector, using the 2-step process
	// DMGlobalToLocalBegin(),DMGlobalToLocalEnd().
	// By placing code between these two statements, computations can be
	// done while messages are in transition.
	ierr = DMGlobalToLocalBegin(da, C, INSERT_VALUES, localC);
	CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da, C, INSERT_VALUES, localC);
	CHKERRQ(ierr);

	// Set the initial values of F
	ierr = VecSet(F, 0.0);
	CHKERRQ(ierr);

	// Compute the new concentrations
	auto& solverHandler = Solver::getSolverHandler();
	solverHandler.updateConcentration(ts, localC, F, ftime);

	// Stop the RHSFunction Timer
	RHSFunctionTimer->stop();

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "RHSJacobian")
/*
 Compute the Jacobian entries based on IFunction() and insert them into the matrix
 */
PetscErrorCode RHSJacobian(TS ts, PetscReal ftime, Vec C, Mat A, Mat J,
		void *) {
	// Start the RHSJacobian timer
	RHSJacobianTimer->start();

	PetscErrorCode ierr;

	// Get the matrix from PETSc
	PetscFunctionBeginUser;
	ierr = MatZeroEntries(J);
	CHKERRQ(ierr);
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);
	Vec localC;
	ierr = DMGetLocalVector(da, &localC);
	CHKERRQ(ierr);

	// Get the complete data array
	ierr = DMGlobalToLocalBegin(da, C, INSERT_VALUES, localC);
	CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da, C, INSERT_VALUES, localC);
	CHKERRQ(ierr);

	// Get the solver handler
	auto& solverHandler = Solver::getSolverHandler();

	/* ----- Compute the off-diagonal part of the Jacobian ----- */
	solverHandler.computeOffDiagonalJacobian(ts, localC, J, ftime);

	ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
	CHKERRQ(ierr);
	ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
	CHKERRQ(ierr);

	/* ----- Compute the partial derivatives for the reaction term ----- */
	solverHandler.computeDiagonalJacobian(ts, localC, J, ftime);

	ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
	CHKERRQ(ierr);
	ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
	CHKERRQ(ierr);

	if (A != J) {
		ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
		CHKERRQ(ierr);
		ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
		CHKERRQ(ierr);
	}

	// Stop the RHSJacobian timer
	RHSJacobianTimer->stop();

	PetscFunctionReturn(0);
}

PetscSolver::PetscSolver(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		Solver(registry) {
	RHSFunctionTimer = handlerRegistry->getTimer("RHSFunctionTimer");
	RHSJacobianTimer = handlerRegistry->getTimer("RHSJacobianTimer");
}

PetscSolver::~PetscSolver() {
}

void PetscSolver::setupMesh() {
}

void PetscSolver::initialize() {
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Initialize program
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	PetscInitialize(&numCLIArgs, &CLIArgs, (char*) 0, help);

	return;
}

void PetscSolver::solve(Options &options) {
	// Initial declarations
	PetscErrorCode ierr;
	double scale = 1.5;
	double **solutionArray = nullptr;
	int maxV = options.getMaxV(), maxHe = options.getMaxImpurity();
	std::map<std::string, int> idMap;
	std::vector<double> concVec;
	std::vector<std::vector<double> > padeVec;
	double maxAdaptTime = 0.0;
	auto heBounds = generateBounds(options.getGroupingMin(),
			options.getGroupingMin(), maxHe);
	auto vBounds = generateBounds(options.getGroupingMin(),
			options.getGroupingMin(), maxV);
	// Loop counter
	int loop = 0;

	// Initialize the random number generator
	std::srand(time(NULL));

	// Read the times if the information is in the HDF5 file
	auto fileName = options.getNetworkFilename();
	double time = 0.0, deltaTime = 1.0e-12;
	int tempTimeStep = -2;
	if (!fileName.empty()) {
		if (HDF5Utils::hasConcentrationGroup(fileName, tempTimeStep)) {
			HDF5Utils::readTimes(fileName, tempTimeStep, time, deltaTime);
		}
	}

	// Set the output precision for std::out
	std::cout.precision(16);

	// Initialiaze the converged reason
	TSConvergedReason reason = TS_CONVERGED_USER;
	while (reason == TS_CONVERGED_USER) {
		// Set the size of the network in the options
		options.setMaxV(maxV);
		options.setMaxImpurity(maxHe);

//		if (loop > 0) {
//			auto oldArgc = options.getPetscArgc();
//			auto oldArgv = options.getPetscArgv();
//
//			options.setPetscArgc(oldArgc + 3);
//
//			// The PETSc argv is an array of pointers to C strings.
//			auto argv = new char*[oldArgc + 4];
//			// Create the fake application name
//			std::string appName = "fakeXolotlApplicationNameForPETSc";
//			argv[0] = new char[appName.length() + 1];
//			strcpy(argv[0], appName.c_str());
//
//			// Now loop on the actual PETSc options
//			int idx = 0;
//			for (idx = 0; idx < oldArgc; idx++) {
//				argv[idx] = oldArgv[idx];
//			}
//			idx = oldArgc;
//			std::string addOp = "-snes_type";
//			argv[idx] = new char[addOp.length() + 1];
//			strcpy(argv[idx], addOp.c_str());
//			idx++;
//			addOp = "test";
//			argv[idx] = new char[addOp.length() + 1];
//			strcpy(argv[idx], addOp.c_str());
//			idx++;
//			addOp = "-snes_test_display";
//			argv[idx] = new char[addOp.length() + 1];
//			strcpy(argv[idx], addOp.c_str());
//			idx++;
//			argv[idx] = 0; // null-terminate the array
//
//			options.setPetscArgv(argv);
//			finalize();
//			setCommandLineOptions(options.getPetscArgc(),
//					options.getPetscArgv());
//			initialize();
//		}

		/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		 Create timestepping solver context
		 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
		TS ts;
		ierr = TSCreate(PETSC_COMM_WORLD, &ts);
		checkPetscError(ierr, "PetscSolver::solve: TSCreate failed.");
		ierr = TSSetType(ts, TSARKIMEX);
		checkPetscError(ierr, "PetscSolver::solve: TSSetType failed.");
		ierr = TSARKIMEXSetFullyImplicit(ts, PETSC_TRUE);
		checkPetscError(ierr,
				"PetscSolver::solve: TSARKIMEXSetFullyImplicit failed.");
		ierr = TSSetProblemType(ts, TS_NONLINEAR);
		checkPetscError(ierr, "PetscSolver::solve: TSSetProblemType failed.");
		ierr = TSSetRHSFunction(ts, NULL, RHSFunction, NULL);
		checkPetscError(ierr, "PetscSolver::solve: TSSetRHSFunction failed.");
		ierr = TSSetRHSJacobian(ts, NULL, NULL, RHSJacobian, NULL);
		checkPetscError(ierr, "PetscSolver::solve: TSSetRHSJacobian failed.");
		ierr = TSSetFromOptions(ts);
		checkPetscError(ierr, "PetscSolver::solve: TSSetFromOptions failed.");
		ierr = TSSetTime(ts, time);
		checkPetscError(ierr, "PetscSolver::solve: TSSetTime failed.");
		ierr = TSSetTimeStep(ts, deltaTime);
		checkPetscError(ierr, "PetscSolver::solve: TSSetTimeStep failed.");

		// Set the maximum adapt time
		TSAdapt adapt;
		ierr = TSGetAdapt(ts, &adapt);
		checkPetscError(ierr, "PetscSolver::solve: TSGetAdapt failed.");
		if (xolotlCore::equal(maxAdaptTime, 0.0)) {
			double minAdaptTime = 0.0;
			ierr = TSAdaptGetStepLimits(adapt, &minAdaptTime, &maxAdaptTime);
			checkPetscError(ierr,
					"PetscSolver::solve: TSAdaptGetStepLimits failed.");

		}
		ierr = TSAdaptSetStepLimits(adapt, 0.0, maxAdaptTime);
		checkPetscError(ierr,
				"PetscSolver::solve: TSAdaptSetStepLimits failed.");

		// Create the material factory
		auto materialFactory =
				xolotlFactory::IMaterialFactory::createMaterialFactory(
						options.getMaterial(), options.getDimensionNumber());

		// Initialize it with the options
		materialFactory->initializeMaterial(options);

		bool tempInitOK = xolotlFactory::initializeTempHandler(options);
		if (!tempInitOK) {
			throw std::runtime_error(
					"Unable to initialize temperature from inputs.");
		}
		// Access the temperature handler registry to get the temperature
		auto tempHandler = xolotlFactory::getTemperatureHandler();

		// Create the network handler factory
		auto networkFactory =
				xolotlFactory::IReactionHandlerFactory::createNetworkFactory(
						options.getMaterial());
		// Build a reaction network
		auto networkLoadTimer = handlerRegistry->getTimer("loadNetwork");
		networkLoadTimer->start();
		networkFactory->initializeReactionNetwork(options, handlerRegistry,
				heBounds, vBounds, padeVec, idMap);
		networkLoadTimer->stop();
		auto& network = networkFactory->getNetworkHandler();

		// Initialize and get the solver handler
		bool dimOK = xolotlFactory::initializeDimension(options, network);
		if (!dimOK) {
			throw std::runtime_error(
					"Unable to initialize dimension from inputs.");
		}
		solverHandler = &(xolotlFactory::getSolverHandler());
		getSolverHandler().initializeHandlers(materialFactory, tempHandler,
				options);

		// Create the solver context
		DM da;
		getSolverHandler().createSolverContext(da);

		/*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		 Extract global vector from DMDA to hold solution
		 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
		Vec C;
		ierr = DMCreateGlobalVector(da, &C);
		checkPetscError(ierr,
				"PetscSolver::solve: DMCreateGlobalVector failed.");
		ierr = TSSetDM(ts, da);
		checkPetscError(ierr, "PetscSolver::solve: TSSetDM failed.");
		ierr = TSSetSolution(ts, C);
		checkPetscError(ierr, "PetscSolver::solve: TSSetSolution failed.");

		// Switch on the number of dimensions to set the monitors
		int dim = getSolverHandler().getDimension();
		switch (dim) {
		case 0:
			// One dimension
			ierr = setupPetsc0DMonitor(ts);
			checkPetscError(ierr,
					"PetscSolver::solve: setupPetsc0DMonitor failed.");
			break;
		case 1:
			// One dimension
			ierr = setupPetsc1DMonitor(ts);
			checkPetscError(ierr,
					"PetscSolver::solve: setupPetsc1DMonitor failed.");
			break;
		case 2:
			// Two dimensions
			ierr = setupPetsc2DMonitor(ts);
			checkPetscError(ierr,
					"PetscSolver::solve: setupPetsc2DMonitor failed.");
			break;
		case 3:
			// Three dimensions
			ierr = setupPetsc3DMonitor(ts);
			checkPetscError(ierr,
					"PetscSolver::solve: setupPetsc3DMonitor failed.");
			break;
		default:
			throw std::string(
					"PetscSolver Exception: Wrong number of dimensions "
							"to set the monitors.");
		}

		/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		 Set initial conditions
		 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

		setupInitialConditions(da, C, concVec, idMap);

		/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		 Solve the ODE system
		 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
		if (ts != NULL && C != NULL) {
			ierr = TSSolve(ts, C);
			checkPetscError(ierr, "PetscSolver::solve: TSSolve failed.");
		} else {
			throw std::string(
					"PetscSolver Exception: Unable to solve! Data not configured properly.");
		}

		ierr = TSGetConvergedReason(ts, &reason);
		checkPetscError(ierr,
				"PetscSolver::solve: TSGetConvergedReason failed.");
		if (reason == TS_CONVERGED_USER) {
			// Get the solutionArray
			ierr = DMDAVecGetArrayDOF(da, C, &solutionArray);
			checkPetscError(ierr,
					"PetscSolver::solve: DMDAVecGetArrayDOF failed.");
			// Update the network with it
			network.updateConcentrationsFromArray(solutionArray[0]);
			// Get the cluster coordinates with the highest concentration
			auto coordinates = network.getHighestClusterCoordinates();
			// Print the values and the size of the current network
			std::cout << "Loop #" << loop << " done." << std::endl
					<< "Highest concentration at " << coordinates.first
					<< " He and " << coordinates.second << " V." << std::endl
					<< "Network size: " << maxHe << " x " << maxV << "."
					<< std::endl << std::endl;

			// Create the new bounds
			// Increase the size of the network for next loop
			maxV = maxV * scale;
			maxHe = maxHe * scale;
			maxAdaptTime = maxAdaptTime * scale;
			heBounds = network.generateBounds(0, maxHe);
			vBounds = network.generateBounds(1, maxV);

			// Create the new concentrations
			transformConcentrationVector(network, heBounds, vBounds, idMap,
					concVec, padeVec);

			// Restore the solutionArray
			ierr = DMDAVecRestoreArrayDOF(da, C, &solutionArray);
			checkPetscError(ierr,
					"PetscSolver::solve: DMDAVecRestoreArrayDOF failed.");

			// Get the current times
			ierr = TSGetTime(ts, &time);
			checkPetscError(ierr, "PetscSolver::solve: TSGetTime failed.");
			ierr = TSGetTimeStep(ts, &deltaTime);
			checkPetscError(ierr, "PetscSolver::solve: TSGetTimeStep failed.");
//			time = time - deltaTime;
		}

		/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		 Free work space.
		 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
		ierr = VecDestroy(&C);
		checkPetscError(ierr, "PetscSolver::solve: VecDestroy failed.");
		ierr = DMDestroy(&da);
		checkPetscError(ierr, "PetscSolver::solve: DMDestroy failed.");
		ierr = TSDestroy(&ts);
		checkPetscError(ierr, "PetscSolver::solve: TSDestroy failed.");

		// Increase the loop count
		loop++;
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Free work space.
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	return;
}

void PetscSolver::finalize() {
	PetscErrorCode ierr;

	ierr = PetscFinalize();
	checkPetscError(ierr, "PetscSolver::finalize: PetscFinalize failed.");

	return;
}

} /* end namespace xolotlSolver */
