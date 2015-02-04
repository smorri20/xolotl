// Includes
#include "PetscSolver.h"
#include <xolotlPerf.h>
#include <VizHandlerRegistryFactory.h>
#include <PlotType.h>
#include <CvsXDataProvider.h>
#include <CvsXYDataProvider.h>
#include <LabelProvider.h>
#include <petscts.h>
#include <petscsys.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <memory>
#include <HDF5Utils.h>

namespace xolotlSolver {

/* ----- Error Handling Code ----- */

/**
 * This operation checks a Petsc error code and converts it to a bool.
 * @param errorCode The Petsc error code.
 * @return True if everything is OK, false otherwise.
 */
static inline bool checkPetscError(PetscErrorCode errorCode) {
	CHKERRQ(errorCode);
}

// Declaration of the functions defined in Monitor.cpp
extern PetscErrorCode monitorTime(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);
extern PetscErrorCode computeHeliumFluence(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);
extern PetscErrorCode monitorPerf(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);

// Declaration of the variables defined in Monitor.cpp
extern std::shared_ptr<xolotlViz::IPlot> perfPlot;
extern double previousTime;

//! How often HDF5 file is written
PetscInt hdf5Stride2D = 0;
//! HDF5 output file name
std::string hdf5OutputName2D = "xolotlStop.h5";
// Declare the vector that will store the Id of the helium clusters
std::vector<int> heIndices2D;
// Declare the vector that will store the weight of the helium clusters
// (their He composition)
std::vector<int> heWeights2D;
//! The pointer to the 2D plot used in MonitorSurface.
std::shared_ptr<xolotlViz::IPlot> surfacePlot2D;
//! The variable to store the interstitial flux at the previous time step.
std::vector<double> previousIFlux2D;
//! The variable to store the total number of interstitials going through the surface.
std::vector<double> nInterstitial2D;

/**
 * This is a monitoring method that will save an hdf5 file at each time step.
 * HDF5 is handling the parallel part, so no call to MPI here.
 */
PetscErrorCode startStop2D(TS ts, PetscInt timestep, PetscReal time, Vec solution,
		void *ictx) {
	PetscErrorCode ierr;
	double ***solutionArray, *gridPointSolution;
	Vec localSolution;
	int xs, xm, Mx, ys, ym, My;

	PetscFunctionBeginUser;

	// Don't do anything if it is not on the stride
	if (timestep % hdf5Stride2D != 0)
		PetscFunctionReturn(0);

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID (important when it is running in parallel)
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);

	// Get the local vector, which is capital when running in parallel,
	// and put it into solutionArray
	ierr = DMGetLocalVector(da, &localSolution);
	checkPetscError(ierr);
	ierr = DMGlobalToLocalBegin(da, solution, INSERT_VALUES, localSolution);
	checkPetscError(ierr);
	ierr = DMGlobalToLocalEnd(da, solution, INSERT_VALUES, localSolution);
	checkPetscError(ierr);
	ierr = DMDAVecGetArrayDOF(da, localSolution, &solutionArray);
	checkPetscError(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	checkPetscError(ierr);
	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	checkPetscError(ierr);

	// Get the solver handler
	auto solverHandler = PetscSolver::getSolverHandler();

	// Get the network
	auto network = solverHandler->getNetwork();

	// Network size
	const int networkSize = network->size();

	// Setup step size variable
	double h = solverHandler->getStepSize();

	// Open the already created HDF5 file
	xolotlCore::HDF5Utils::openFile(hdf5OutputName2D);

	// Get the current time step
	double currentTimeStep;
	ierr = TSGetTimeStep(ts, &currentTimeStep);
	checkPetscError(ierr);

	// Add a concentration sub group
	xolotlCore::HDF5Utils::addConcentrationSubGroup(timestep, networkSize, time,
			currentTimeStep);

	// Loop on the full grid
	for (int j = 0; j < My; j++) {
		for (int i = 0; i < Mx; i++) {
			// Size of the concentration that will be stored
			int concSize = -1;
			// Vector for the concentrations
			std::vector<std::vector<double> > concVector;

			// If it is the locally owned part of the grid
			if (i >= xs && i < xs + xm && j >= ys && j < ys + ym) {
				// Get the pointer to the beginning of the solution data for this grid point
				gridPointSolution = solutionArray[j][i];

				// Loop on the concentrations
				concVector.clear();
				for (int l = 0; l < networkSize; l++) {
					if (gridPointSolution[l] > 1.0e-16) {
						// Create the concentration vector for this cluster
						std::vector<double> conc;
						conc.push_back((double) l);
						conc.push_back(gridPointSolution[l]);

						// Add it to the main vector
						concVector.push_back(conc);
					}
				}

				// Send the size of the vector to the other processes
				concSize = concVector.size();
				// Loop on all the processes
				for (int l = 0; l < worldSize; l++) {
					// Skip its own
					if (l == procId)
						continue;

					// Send the size
					MPI_Send(&concSize, 1, MPI_INT, l, 0, MPI_COMM_WORLD);
				}
			}

			// Else: only receive the conc size
			else {
				MPI_Recv(&concSize, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE);
			}

			// Skip the grid point if the size is 0
			if (concSize == 0)
				continue;

			// All processes must create the dataset
			xolotlCore::HDF5Utils::addConcentrationDataset(concSize, i, j);

			// If it is the locally owned part of the grid
			if (i >= xs && i < xs + xm && j >= ys && j < ys + ym) {
				// Fill the dataset
				xolotlCore::HDF5Utils::fillConcentrations(concVector, i, j);
			}
		}
	}

	// Finalize the HDF5 file
	xolotlCore::HDF5Utils::closeFile();

	PetscFunctionReturn(0);
}

/**
 * This is a monitoring method that will compute the total helium fluence
 */
PetscErrorCode computeHeliumRetention2D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	PetscErrorCode ierr;
	int xs, xm, ys, ym;

	PetscFunctionBeginUser;

	// Get the solver handler
	auto solverHandler = PetscSolver::getSolverHandler();

	// Get the flux handler that will be used to compute fluxes.
	auto fluxHandler = solverHandler->getFluxHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	checkPetscError(ierr);

	// Setup step size variable
	double h = solverHandler->getStepSize();

	// Get the array of concentration
	PetscReal ***solutionArray;
	ierr = DMDAVecGetArrayDOF(da, solution, &solutionArray);
	checkPetscError(ierr);

	// Store the concentration over the grid
	double heConcentration = 0;

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal *gridPointSolution;

	// Loop on the grid
	for (int j = ys; j < ys + ym; j++) {
		for (int i = xs; i < xs + xm; i++) {
			// Get the pointer to the beginning of the solution data for this grid point
			gridPointSolution = solutionArray[j][i];

			// Loop on all the indices
			for (int l = 0; l < heIndices2D.size(); l++) {
				// Add the current concentration times the number of helium in the cluster
				// (from the weight vector)
				heConcentration += gridPointSolution[heIndices2D[l]] * heWeights2D[l] * h;
			}
		}
	}

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);
	// Get the current process ID
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	// Master process
	if (procId == 0) {
		// Loop on all the other processes
		for (int i = 1; i < worldSize; i++) {
			double otherConcentration = 0.0;

			// Receive the value from the other processes
			MPI_Recv(&otherConcentration, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);

			// Add them to the master one
			heConcentration += otherConcentration;
		}

		// Get the total size of the grid rescale the concentrations
		int Mx, My;
		ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE);
		checkPetscError(ierr);

		// Compute the total surface irradiated by the helium flux
		double surface = (My - 2) * h;

		// Rescale the concentration
		heConcentration = heConcentration / surface;

		// Get the fluence
		double heliumFluence = fluxHandler->getHeFluence();

		// Print the result
		std::cout << "\nTime: " << time << std::endl;
		std::cout << "Helium retention = "
				<< 100.0 * (heConcentration / heliumFluence) << " %"
				<< std::endl;
		std::cout << "Helium mean concentration = " << heConcentration << std::endl;
		std::cout << "Helium fluence = " << heliumFluence << "\n" << std::endl;

//		// Uncomment to write the retention and the fluence in a file
//		std::ofstream outputFile;
//		outputFile.open("retentionOut.txt", ios::app);
//		outputFile << heliumFluence << " "
//				<< 100.0 * (heConcentration / heliumFluence) << std::endl;
//		outputFile.close();
	}

	else {
		// Send the value of the timer to the master process
		MPI_Send(&heConcentration, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}

	PetscFunctionReturn(0);
}

/**
 * This is a monitoring method that will save 2D plots of the concentration of
 * a specific cluster at each grid point.
 */
PetscErrorCode monitorSurface2D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	PetscErrorCode ierr;
	double ***solutionArray, *gridPointSolution, x, y;
	Vec localSolution;
	int xs, xm, Mx, ys, ym, My;

	PetscFunctionBeginUser;

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);

	// Get the local vector, which is capital when running in parallel,
	// and put it into solutionArray
	ierr = DMGetLocalVector(da, &localSolution);
	checkPetscError(ierr);
	ierr = DMGlobalToLocalBegin(da, solution, INSERT_VALUES, localSolution);
	checkPetscError(ierr);
	ierr = DMGlobalToLocalEnd(da, solution, INSERT_VALUES, localSolution);
	checkPetscError(ierr);
	ierr = DMDAVecGetArrayDOF(da, localSolution, &solutionArray);
	checkPetscError(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	checkPetscError(ierr);
	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	checkPetscError(ierr);

	// Get the solver handler
	auto solverHandler = PetscSolver::getSolverHandler();

	// Get the network
	auto network = solverHandler->getNetwork();

	// Setup step size variable
	double h = solverHandler->getStepSize();

	// Choice of the cluster to be plotted
	int iCluster = 6;

	// Create a Point vector to store the data to give to the data provider
	// for the visualization
	auto myPoints = std::make_shared<std::vector<xolotlViz::Point> >();
	// Create a point here so that it is not created and deleted in the loop
	xolotlViz::Point thePoint;

	// Loop on the full grid
	for (int j = 0; j < My; j++) {
		for (int i = 0; i < Mx; i++) {
			// If it is the locally owned part of the grid
			if (i >= xs && i < xs + xm && j >= ys && j < ys + ym) {
				// Get the pointer to the beginning of the solution data for this grid point
				gridPointSolution = solutionArray[j][i];
				// Compute x and y
				x = i * h;
				y = j * h;

				// If it is procId 0 just store the value in the myPoints vector
				if (procId == 0) {
					// Modify the Point with the gridPointSolution[iCluster] as the value
					// and add it to myPoints
					thePoint.value = gridPointSolution[iCluster];
					thePoint.t = time;
					thePoint.x = x;
					thePoint.y = y;
					myPoints->push_back(thePoint);
				}
				// Else, the values must be sent to procId 0
				else {
					// Send the value of the local position to the master process
					MPI_Send(&x, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
					// Send the value of the local position to the master process
					MPI_Send(&y, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);

					// Send the value of the concentration to the master process
					MPI_Send(&gridPointSolution[iCluster], 1, MPI_DOUBLE, 0, 2,
							MPI_COMM_WORLD);
				}
			}
			// Else if it is NOT the locally owned part of the grid but still procId == 0,
			// it should receive the values for the point and add them to myPoint
			else if (procId == 0) {
				// Get the position
				MPI_Recv(&x, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE);
				MPI_Recv(&y, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE);

				// and the concentration
				double conc = 0.0;
				MPI_Recv(&conc, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE);

				// Modify the Point with the received values and add it to myPoints
				thePoint.value = conc;
				thePoint.t = time;
				thePoint.x = x;
				thePoint.y = y;
				myPoints->push_back(thePoint);
			}

			// Wait for everybody at each grid point
			MPI_Barrier(PETSC_COMM_WORLD);
		}
	}

	// Plot everything from procId == 0
	if (procId == 0) {
		// Get the data provider and give it the points
		surfacePlot2D->getDataProvider()->setPoints(myPoints);

		// Get the iCluster cluster to have access to its name
		auto reactants = network->getAll();
		auto cluster = (PSICluster *) reactants->at(iCluster);

		// Change the title of the plot and the name of the data
		std::stringstream title;
		title << cluster->getName();
		surfacePlot2D->getDataProvider()->setDataName(title.str());
		title << " concentration";
		surfacePlot2D->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::stringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		surfacePlot2D->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);
		checkPetscError(ierr);
		// Give the timestep to the label provider
		std::stringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
				<< "s";
		surfacePlot2D->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::stringstream fileName;
		fileName << cluster->getName() << "_surface_TS" << timestep << ".pnm";
		surfacePlot2D->write(fileName.str());
	}

	PetscFunctionReturn(0);
}

/**
 * This is a monitoring method that will compute the flux of interstitials
 * at the surface
 */
PetscErrorCode monitorInterstitial2D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	PetscErrorCode ierr;
	double ***solutionArray, *gridPointSolution, x;
	Vec localSolution;
	int xs, xm, xi, ys, ym, yj, Mx, My;

	PetscFunctionBeginUser;

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);

	// Get the local vector, which is capital when running in parallel,
	// and put it into solutionArray
	ierr = DMGetLocalVector(da, &localSolution);
	checkPetscError(ierr);
	ierr = DMGlobalToLocalBegin(da, solution, INSERT_VALUES, localSolution);
	checkPetscError(ierr);
	ierr = DMGlobalToLocalEnd(da, solution, INSERT_VALUES, localSolution);
	checkPetscError(ierr);
	ierr = DMDAVecGetArrayDOF(da, localSolution, &solutionArray);
	checkPetscError(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	checkPetscError(ierr);
	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	checkPetscError(ierr);

	// Get the solver handler
	auto solverHandler = PetscSolver::getSolverHandler();

	// Get the network
	auto network = solverHandler->getNetwork();
	// Get all the interstitial clusters
	auto interstitials = network->getAll("I");

	// Setup step size variable
	double h = solverHandler->getStepSize();
	double s = 1.0 / (h * h);

	// Initialize the boolean to know if the flux need to be reinitialized
	bool initializeFlux = false;

	// Loop on the possible yj
	for (yj = 1; yj < My - 1; yj++) {
		// Get the position of the surface at yj
		int surfacePos = solverHandler->getSurfacePosition(yj);
		xi = surfacePos + 1;

		// if xi, yj is on this process
		if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {
			// Get the concentrations at xi = surfacePos[0] + 1
			gridPointSolution = solutionArray[yj][xi];

			// Get the delta time from the previous timestep to this timestep
			double dt = time - previousTime;

			// Compute the total density of intersitials that escaped from the
			// surface since last timestep using the stored flux
			nInterstitial2D[yj] += previousIFlux2D[yj] * dt;


			// Initialize the value for the flux
			double newFlux = 0.0;

			// Loop on all the interstitial clusters
			for (int i = 0; i < interstitials.size(); i++) {
				// Get the cluster
				auto cluster = (PSICluster *) interstitials.at(i);
				// Get its id and concentration
				int id = cluster->getId() - 1;
				double conc = gridPointSolution[id];
				// Get its size and diffusion coefficient
				int size = cluster->getSize();
				double coef = cluster->getDiffusionCoefficient();

				// Compute the flux
				newFlux += (double) size * s * coef * conc;
			}

			previousIFlux2D[yj] = newFlux;

			// Send the information about nInterstitial2D and previousFlux2D
			// to the other processes
			// Loop on all the processes
			for (int i = 0; i < worldSize; i++) {
				// Skip this process
				if (i == procId) continue;

				// Send nInterstitial
				MPI_Send(&nInterstitial2D[yj], 1, MPI_DOUBLE, i, 4, MPI_COMM_WORLD);
				// Send previousFlux
				MPI_Send(&previousIFlux2D[yj], 1, MPI_DOUBLE, i, 4, MPI_COMM_WORLD);
			}
		}

		// xi, yj is not on this process, but the process needs to know what is the
		// flux and interstitial value from the other process
		else {
			// Receive nInterstitial
			MPI_Recv(&nInterstitial2D[yj], 1, MPI_DOUBLE, MPI_ANY_SOURCE, 4, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);
			// Receive previousFlux
			MPI_Recv(&previousIFlux2D[yj], 1, MPI_DOUBLE, MPI_ANY_SOURCE, 4, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);
		}

		// Wait for everybody at each grid point
		MPI_Barrier(PETSC_COMM_WORLD);

		// Now that all the processes have the same value of nInterstitials, compare
		// it to the threshold to now if we should move the surface

		// The density of tungsten is 62.8 atoms/nm3, thus the threshold is
		double threshold = 62.8 * h * h;
		if (nInterstitial2D[yj] > threshold) {
			// Compute the number of grid points to move the surface of
			int nGridPoints = (int) (nInterstitial2D[yj] / threshold);

			// Remove the number of interstitials we just transformed in new material
			// from nInterstitial2D
			nInterstitial2D[yj] = nInterstitial2D[yj] - threshold * (double) nGridPoints;

			// Compute the new surface position
			surfacePos -= nGridPoints;

			// Throw an exception if the position is negative
			if (surfacePos < 0) {
				throw std::string(
						"\nxolotlSolver::Monitor2D: The surface is trying to go outside of the grid!!");
			}

			// Printing information about the extension of the material
			if (procId == 0) {
				std::cout << "Adding " << nGridPoints << " points to the grid on yj = " << yj
						<< " at time: " << time << " s." << std::endl;
			}

			// Set it in the solver
			solverHandler->setSurfacePosition(surfacePos, yj);

			// Tell it the surface has moved
			solverHandler->changeSurfacePosition();

			// The flux will need to be initialized
			initializeFlux = true;
		}
	}

	// If we need to initialize the flux
	if (initializeFlux) {
		// Compute the mean value of the surface position
		int meanPosition = 0;
		for (int l = 1; l < My - 1; l++) {
			meanPosition += solverHandler->getSurfacePosition(l);
		}
		meanPosition =  meanPosition / (My - 2);

		// Get the flux handler to reinitialize it
		auto fluxHandler = solverHandler->getFluxHandler();

		// Initiliaze the flux with the new surface position
		fluxHandler->initializeFluxHandler(meanPosition, Mx, h, h);
	}

    PetscFunctionReturn(0);
}

/**
 * This operation sets up a monitor that will call monitorSolve
 * @param ts The time stepper
 * @return A standard PETSc error code
 */
PetscErrorCode setupPetsc2DMonitor(TS ts) {
	PetscErrorCode ierr;

	//! The xolotlViz handler registry
	auto vizHandlerRegistry = xolotlFactory::getVizHandlerRegistry();

	// Flags to launch the monitors or not
	PetscBool flagPerf, flagRetention, flagStatus, flag2DPlot;

	// Check the option -plot_perf
	ierr = PetscOptionsHasName(NULL, "-plot_perf", &flagPerf);
	checkPetscError(ierr);

	// Check the option -plot_2d
	ierr = PetscOptionsHasName(NULL, "-plot_2d", &flag2DPlot);
	checkPetscError(ierr);

	// Check the option -helium_retention
	ierr = PetscOptionsHasName(NULL, "-helium_retention", &flagRetention);
	checkPetscError(ierr);

	// Check the option -start_stop
	ierr = PetscOptionsHasName(NULL, "-start_stop", &flagStatus);
	checkPetscError(ierr);

	// Get the solver handler
	auto solverHandler = PetscSolver::getSolverHandler();

	// Get the network and its size
	auto network = solverHandler->getNetwork();
	const int networkSize = network->size();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);

	// Get the total size of the grid
	int Mx, My;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	checkPetscError(ierr);

	// Set the monitor to save performance plots (has to be in parallel)
	if (flagPerf) {
		// Create a ScatterPlot
		perfPlot = vizHandlerRegistry->getPlot("perfPlot",
				xolotlViz::PlotType::SCATTER);

		// Create and set the label provider
		auto labelProvider = std::make_shared<xolotlViz::LabelProvider>(
				"labelProvider");
		labelProvider->axis1Label = "Process ID";
		labelProvider->axis2Label = "Solver Time";

		// Give it to the plot
		perfPlot->setLabelProvider(labelProvider);

		// Create the data provider
		auto dataProvider = std::make_shared<xolotlViz::CvsXDataProvider>(
				"dataProvider");

		// Give it to the plot
		perfPlot->setDataProvider(dataProvider);

		// monitorPerf will be called at each timestep
		ierr = TSMonitorSet(ts, monitorPerf, NULL, NULL);
		checkPetscError(ierr);

	}

	// Set the monitor to compute the helium fluence for the retention calculation
	if (flagRetention) {
		// Get all the helium clusters
		auto heClusters = network->getAll(heType);

		// Get all the helium-vacancy clusters
		auto heVClusters = network->getAll(heVType);

		// Loop on the helium clusters
		for (int i = 0; i < heClusters.size(); i++) {
			auto cluster = (PSICluster *) heClusters[i];
			int id = cluster->getId() - 1;
			// Add the Id to the vector
			heIndices2D.push_back(id);
			// Add the number of heliums of this cluster to the weight
			heWeights2D.push_back(cluster->getSize());
		}

		// Loop on the helium-vacancy clusters
		for (int i = 0; i < heVClusters.size(); i++) {
			auto cluster = (PSICluster *) heVClusters[i];
			int id = cluster->getId() - 1;
			// Add the Id to the vector
			heIndices2D.push_back(id);
			// Add the number of heliums of this cluster to the weight
			auto comp = cluster->getComposition();
			heWeights2D.push_back(comp[heType]);
		}

		if (heIndices2D.size() == 0) {
			throw std::string(
					"PetscSolver Exception: Cannot compute the retention because there is "
					"no helium or helium-vacancy cluster in the network.");
		}

		// computeHeliumFluence will be called at each timestep
		ierr = TSMonitorSet(ts, computeHeliumFluence, NULL, NULL);
		checkPetscError(ierr);

		// computeHeliumRetention2D will be called at each timestep
		ierr = TSMonitorSet(ts, computeHeliumRetention2D, NULL, NULL);
		checkPetscError(ierr);

//		// Uncomment to clear the file where the retention will be written
//		std::ofstream outputFile;
//		outputFile.open("retentionOut.txt");
//		outputFile.close();
	}

	// Set the monitor to save the status of the simulation in hdf5 file
	if (flagStatus) {
		// Find the stride to know how often the HDF5 file has to be written
		PetscBool flag;
		ierr = PetscOptionsGetInt(NULL, "-start_stop", &hdf5Stride2D, &flag);
		checkPetscError(ierr);
		if (!flag)
			hdf5Stride2D = 1;

		// Initialize the HDF5 file for all the processes
		xolotlCore::HDF5Utils::initializeFile(hdf5OutputName2D, networkSize);

		// Get the solver handler
		auto solverHandler = PetscSolver::getSolverHandler();

		// Setup step size variable
		double h = solverHandler->getStepSize();

		// Get the refinement of the grid
		PetscInt refinement = 0;
		ierr = PetscOptionsGetInt(NULL, "-da_refine", &refinement, &flag);
		checkPetscError(ierr);
		if (!flag)
			refinement = 0;

		// Save the header in the HDF5 file
		xolotlCore::HDF5Utils::fillHeader(2, Mx, h, My, h);

		// Save the network in the HDF5 file
		xolotlCore::HDF5Utils::fillNetwork(network);

		// Finalize the HDF5 file
		xolotlCore::HDF5Utils::finalizeFile();

		// startStop2D will be called at each timestep
		ierr = TSMonitorSet(ts, startStop2D, NULL, NULL);
		checkPetscError(ierr);
	}

	// Set the monitor to save surface plots of clusters concentration
	if (flag2DPlot) {
		// Create a SurfacePlot
		surfacePlot2D = vizHandlerRegistry->getPlot("surfacePlot2D",
				xolotlViz::PlotType::SURFACE);

		// Create and set the label provider
		auto labelProvider = std::make_shared<xolotlViz::LabelProvider>(
				"labelProvider");
		labelProvider->axis1Label = "Depth (nm)";
		labelProvider->axis2Label = "Y (nm)";
		labelProvider->axis3Label = "Concentration";

		// Give it to the plot
		surfacePlot2D->setLabelProvider(labelProvider);

		// Create the data provider
		auto dataProvider = std::make_shared<xolotlViz::CvsXYDataProvider>(
				"dataProvider");

		// Give it to the plot
		surfacePlot2D->setDataProvider(dataProvider);

		// monitorSurface2D will be called at each timestep
		ierr = TSMonitorSet(ts, monitorSurface2D, NULL, NULL);
		checkPetscError(ierr);
	}

	// Initialize nInterstitial2D and previousIFlux2D before monitoring the
	// interstitial flux
	for (int l = 0; l < My; l++) {
		nInterstitial2D.push_back(0.0);
		previousIFlux2D.push_back(0.0);
	}

	// Set the monitor on the outgoing flux of interstitials at the surface
	// monitorInterstitial2D will be called at each timestep
	ierr = TSMonitorSet(ts, monitorInterstitial2D, NULL, NULL);
	checkPetscError(ierr);

	// Set the monitor to simply change the previous time to the new time
	// monitorTime will be called at each timestep
	ierr = TSMonitorSet(ts, monitorTime, NULL, NULL);
	checkPetscError(ierr);

	PetscFunctionReturn(0);
}

}

/* end namespace xolotlSolver */
