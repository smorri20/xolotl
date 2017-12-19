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
	double scale = max / 25.0, tau = (double)(max - peak) / 7.0;
	std::vector<IReactant::SizeType> bounds;
	bounds.push_back(1);

	// Generate the bounds
	while (value < max) {

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
		std::map<std::string, int> &idMap, std::vector<double> & concVec) {
	// Initial declarations
	int i = 0;
	// Reset the concentration vector and id map
	concVec.clear();
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
							mom1 += concentration / ((double) numHe - mean1);
							mom2 += concentration / ((double) numV - mean2);
						}
					}
				}

				// Rescale the running quantities
				conc = conc
						/ (double) ((bounds1[j + 1] - bounds1[j])
								* (bounds2[k + 1] - bounds2[k]));
				mom1 = (double) (bounds1[j + 1] - bounds1[j] - 1) * mom1
						/ (double) ((bounds1[j + 1] - bounds1[j]) * 2
								* (bounds2[k + 1] - bounds2[k]));
				mom2 = (double) (bounds2[k + 1] - bounds2[k] - 1) * mom2
						/ (double) ((bounds1[j + 1] - bounds1[j]) * 2
								* (bounds2[k + 1] - bounds2[k]));

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
	double maxAdaptTime = 0.0;
	auto heBounds = generateBounds(options.getGroupingMin(),
			options.getGroupingMin(), maxHe);
	auto vBounds = generateBounds(options.getGroupingMin(),
			options.getGroupingMin(), maxV);
	// Loop counter
	int loop = 0;

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
				heBounds, vBounds);
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
					concVec);

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
