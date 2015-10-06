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

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "startStop2D")
/**
 * This is a monitoring method that will save an hdf5 file at each time step.
 * HDF5 is handling the parallel part, so no call to MPI here.
 */
PetscErrorCode startStop2D(TS ts, PetscInt timestep, PetscReal time, Vec solution,
		void *ictx) {
	PetscErrorCode ierr;
	const double ***solutionArray, *gridPointSolution;
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
	ierr = TSGetDM(ts, &da);CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);
	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);CHKERRQ(ierr);

	// Get the solver handler
	auto solverHandler = PetscSolver::getSolverHandler();

	// Get the network
	auto network = solverHandler->getNetwork();

	// Network size
	const int networkSize = network->size();

	// Get the vector of positions of the surface
	std::vector<int> surfaceIndices;
	for (int i = 0; i < My; i++) {
		surfaceIndices.push_back(solverHandler->getSurfacePosition(i));
	}

	// Open the already created HDF5 file
	xolotlCore::HDF5Utils::openFile(hdf5OutputName2D);

	// Get the current time step
	double currentTimeStep;
	ierr = TSGetTimeStep(ts, &currentTimeStep);CHKERRQ(ierr);

	// Add a concentration sub group
	xolotlCore::HDF5Utils::addConcentrationSubGroup(timestep, networkSize, time,
			currentTimeStep);

	// Write the surface positions in the concentration sub group
	xolotlCore::HDF5Utils::writeSurface2D(timestep, surfaceIndices);

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

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeHeliumRetention2D")
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
	ierr = TSGetDM(ts, &da);CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);

	// Get the physical grid in the x direction
	auto grid = solverHandler->getXGrid();

	// Setup step size variables
	double hy = solverHandler->getStepSizeY();

	// Get the array of concentration
	const double ***solutionArray, *gridPointSolution;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);CHKERRQ(ierr);

	// Store the concentration over the grid
	double heConcentration = 0.0;

	// Loop on the grid
	for (int j = ys; j < ys + ym; j++) {
		for (int i = xs; i < xs + xm; i++) {
			// Get the pointer to the beginning of the solution data for this grid point
			gridPointSolution = solutionArray[j][i];

			// Loop on all the indices
			for (int l = 0; l < heIndices2D.size(); l++) {
				// Add the current concentration times the number of helium in the cluster
				// (from the weight vector)
				heConcentration += gridPointSolution[heIndices2D[l]] * heWeights2D[l]
				                                     * (grid[i] - grid[i-1]) * hy;
			}
		}
	}

	// Get the current process ID
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	// Sum all the concentrations through MPI reduce
	double totalHeConcentration = 0.0;
	MPI_Reduce(&heConcentration, &totalHeConcentration, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	// Master process
	if (procId == 0) {
		// Get the total size of the grid rescale the concentrations
		int Mx, My;
		ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
		PETSC_IGNORE);CHKERRQ(ierr);

		// Compute the total surface irradiated by the helium flux
		double surface = (double) My * hy;

		// Rescale the concentration
		totalHeConcentration = totalHeConcentration / surface;

		// Get the fluence
		double heliumFluence = fluxHandler->getFluence();

		// Print the result
		std::cout << "\nTime: " << time << std::endl;
		std::cout << "Helium retention = "
				<< 100.0 * (totalHeConcentration / heliumFluence) << " %"
				<< std::endl;
		std::cout << "Helium mean concentration = " << totalHeConcentration << std::endl;
		std::cout << "Helium fluence = " << heliumFluence << "\n" << std::endl;

//		// Uncomment to write the retention and the fluence in a file
//		std::ofstream outputFile;
//		outputFile.open("retentionOut.txt", ios::app);
//		outputFile << heliumFluence << " "
//				<< 100.0 * (totalHeConcentration / heliumFluence) << std::endl;
//		outputFile.close();
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeHeliumConc2D")
/**
 * This is a monitoring method that will compute the helium concentrations
 */
PetscErrorCode computeHeliumConc2D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	PetscErrorCode ierr;
	int xs, xm, ys, ym;

	PetscFunctionBeginUser;

	// Get the number of processes
	int worldSize;
	MPI_Comm_size(PETSC_COMM_WORLD, &worldSize);

	// Gets the process ID (important when it is running in parallel)
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	// Get the solver handler
	auto solverHandler = PetscSolver::getSolverHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);

	// Get the physical grid in the x direction
	auto grid = solverHandler->getXGrid();

	// Setup step size variables
	double hy = solverHandler->getStepSizeY();

	// Get the total size of the grid rescale the concentrations
	int Mx, My;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);CHKERRQ(ierr);

	// Get the array of concentration
	const double ***solutionArray, *gridPointSolution;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);CHKERRQ(ierr);

	// Store the concentration over the grid
	std::vector<double> heConcLocal;
	std::vector<double> heConcentrations;
	// Initialize
	for (int i = 0; i < 1001; i++) {
		heConcLocal.push_back(0.0);
		heConcentrations.push_back(0.0);
	}

	// Open the file
	std::ofstream outputFile;
	if (procId == 0) {
		std::stringstream name;
		name << "heliumConc_" << timestep << ".dat";
		outputFile.open(name.str());
	}

	// Loop on the full grid
	for (int xi = 0; xi < Mx; xi++) {
		// Wait for everybody at each grid point
		MPI_Barrier(PETSC_COMM_WORLD);

		// Set x
		double x = grid[xi];

		// Loop on Y to integrate
		for (int yj = 0; yj < My; yj++) {
			// If we are on the right process
			if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {
				// Get the pointer to the beginning of the solution data for this grid point
				gridPointSolution = solutionArray[yj][xi];

				// Loop on all the indices
				for (int l = 0; l < heIndices2D.size(); l++) {
					// Add the current concentration
					heConcLocal[heWeights2D[l]] += gridPointSolution[heIndices2D[l]]
																		  * (grid[xi] - grid[xi-1]) * hy;
				}
			}
		}

		// Gather all the data
		for (int i = 0; i < 1001; i++) {
			MPI_Reduce(&heConcLocal[i], &heConcentrations[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		}

		// Print it from the main proc
		if (procId == 0) {
			for (int i = 0; i < 1001; i++) {
				if (heConcentrations[i] > 1.0e-16) {
					outputFile << x << " " << i << " " << heConcentrations[i] << std::endl;
				}
			}
		}

		// Reinitialize the concentrations
		for (int i = 0; i < 1001; i++) {
			heConcLocal[i] = 0.0;
			heConcentrations[i] = 0.0;
		}
	}

	// Close the file
	outputFile.close();

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "monitorSurface2D")
/**
 * This is a monitoring method that will save 2D plots of the concentration of
 * a specific cluster at each grid point.
 */
PetscErrorCode monitorSurface2D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	PetscErrorCode ierr;
	const double ***solutionArray, *gridPointSolution;
	int xs, xm, Mx, ys, ym, My;
	double x, y;

	PetscFunctionBeginUser;

	// Gets the process ID
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);
	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);CHKERRQ(ierr);

	// Get the solver handler
	auto solverHandler = PetscSolver::getSolverHandler();

	// Get the network
	auto network = solverHandler->getNetwork();

	// Get the physical grid in the x direction
	auto grid = solverHandler->getXGrid();

	// Setup step size variables
	double hy = solverHandler->getStepSizeY();

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
				x = grid[i];
				y = (double) j * hy;

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
					MPI_Send(&x, 1, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
					// Send the value of the local position to the master process
					MPI_Send(&y, 1, MPI_DOUBLE, 0, 11, MPI_COMM_WORLD);

					// Send the value of the concentration to the master process
					MPI_Send(&gridPointSolution[iCluster], 1, MPI_DOUBLE, 0, 12,
							MPI_COMM_WORLD);
				}
			}
			// Else if it is NOT the locally owned part of the grid but still procId == 0,
			// it should receive the values for the point and add them to myPoint
			else if (procId == 0) {
				// Get the position
				MPI_Recv(&x, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 10, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE);
				MPI_Recv(&y, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 11, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE);

				// and the concentration
				double conc = 0.0;
				MPI_Recv(&conc, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 12, MPI_COMM_WORLD,
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
		ierr = TSGetTimeStep(ts, &currentTimeStep);CHKERRQ(ierr);
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

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

/**
* This is a monitoring method that will compute the flux of interstitials
* at the surface and move the position of the surface if necessary.
*/
PetscErrorCode monitorInterstitial2D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	PetscErrorCode ierr;
	double ***solutionArray, *gridPointSolution;
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
	ierr = TSGetDM(ts, &da);CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);CHKERRQ(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);

	// Get the size of the total grid
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE);CHKERRQ(ierr);

	// Get the solver handler
	auto solverHandler = PetscSolver::getSolverHandler();

	// Get the network
	auto network = solverHandler->getNetwork();

	// Get all the interstitial clusters
	auto interstitials = network->getAll("I");
	// Get the single vacancy ID
	auto singleVacancyCluster = network->get(xolotlCore::vType, 1);
	int vacancyIndex = -1;
	if (singleVacancyCluster)
		vacancyIndex = singleVacancyCluster->getId() - 1;

	// Get the physical grid
	auto grid = solverHandler->getXGrid();

	// Setup step size variables
	double hy = solverHandler->getStepSizeY();

	// Get the initial vacancy concentration
	double initialVConc = solverHandler->getInitialVConc();

	// Get the delta time from the previous timestep to this timestep
	double dt = time - previousTime;

	// Loop on the possible yj
	for (yj = 0; yj < My; yj++) {
		// Get the position of the surface at yj
		int surfacePos = solverHandler->getSurfacePosition(yj);
		xi = surfacePos + 1;

		// Value to now on which processor is the location of the surface,
		// for MPI usage
		int surfaceProc = 0;

		// if xi, yj are on this process
		if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym) {
			// Get the concentrations at xi = surfacePos + 1
			gridPointSolution = solutionArray[yj][xi];

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

				// Factor for finite difference
				double hxLeft = grid[xi] - grid[xi-1];
				double hxRight = grid[xi+1] - grid[xi];
				double factor = 2.0 / (hxLeft * (hxLeft + hxRight));
				// Compute the flux going to the left
				newFlux += (double) size * factor * coef * conc;
			}

			// Update the previous flux at this position
			previousIFlux2D[yj] = newFlux;

			// Set the surface processor
			surfaceProc = procId;
		}

		// Get which processor will send the information
		int surfaceId = 0;
		MPI_Allreduce(&surfaceProc, &surfaceId, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		// Send the information about nInterstitial2D and previousFlux2D
		// to the other processes
		MPI_Bcast(&nInterstitial2D[yj], 1, MPI_DOUBLE, surfaceId, MPI_COMM_WORLD);
		MPI_Bcast(&previousIFlux2D[yj], 1, MPI_DOUBLE, surfaceId, MPI_COMM_WORLD);

		// Now that all the processes have the same value of nInterstitials, compare
		// it to the threshold to now if we should move the surface

		// The density of tungsten is 62.8 atoms/nm3, thus the threshold is
		double threshold = (62.8 - initialVConc) * (grid[xi] - grid[xi-1]);
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

			// Initialize the vacancy concentration on the new grid points
			// Loop on the new grid points
			while (nGridPoints > 0) {
				// Position of the newly created grid point
				xi = surfacePos + nGridPoints;

				// If xi and yj are on this process
				if (xi >= xs && xi < xs + xm && yj >= ys && yj < ys + ym && vacancyIndex > 0) {
					// Get the concentrations
					gridPointSolution = solutionArray[yj][xi];
					// Initialize the vacancy concentration
					gridPointSolution[vacancyIndex] = initialVConc;
				}

				// Decrease the number of grid points
				--nGridPoints;
			}
		}
	}

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

/**
 * This operation sets up a monitor that will call monitorSolve
 * @param ts The time stepper
 * @return A standard PETSc error code
 */
PetscErrorCode setupPetsc2DMonitor(TS ts) {
	PetscErrorCode ierr;

	// Get the process ID
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	// Get the xolotlViz handler registry
	auto vizHandlerRegistry = xolotlFactory::getVizHandlerRegistry();

	// Flags to launch the monitors or not
	PetscBool flagPerf, flagRetention, flagStatus, flag2DPlot, flagConc;

	// Check the option -plot_perf
	ierr = PetscOptionsHasName(NULL, "-plot_perf", &flagPerf);
	checkPetscError(ierr, "setupPetsc2DMonitor: PetscOptionsHasName (-plot_perf) failed.");

	// Check the option -plot_2d
	ierr = PetscOptionsHasName(NULL, "-plot_2d", &flag2DPlot);
	checkPetscError(ierr, "setupPetsc2DMonitor: PetscOptionsHasName (-plot_2d) failed.");

	// Check the option -helium_retention
	ierr = PetscOptionsHasName(NULL, "-helium_retention", &flagRetention);
	checkPetscError(ierr, "setupPetsc2DMonitor: PetscOptionsHasName (-helium_retention) failed.");

	// Check the option -start_stop
	ierr = PetscOptionsHasName(NULL, "-start_stop", &flagStatus);
	checkPetscError(ierr, "setupPetsc2DMonitor: PetscOptionsHasName (-start_stop) failed.");

	// Check the option -helium_conc
	ierr = PetscOptionsHasName(NULL, "-helium_conc", &flagConc);
	checkPetscError(ierr, "setupPetsc2DMonitor: PetscOptionsHasName (-helium_conc) failed.");

	// Get the solver handler
	auto solverHandler = PetscSolver::getSolverHandler();

	// Get the network and its size
	auto network = solverHandler->getNetwork();
	const int networkSize = network->size();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);CHKERRQ(ierr);
	checkPetscError(ierr, "setupPetsc2DMonitor: TSGetDM failed.");

	// Get the total size of the grid
	int Mx, My;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE);CHKERRQ(ierr);
	checkPetscError(ierr, "setupPetsc2DMonitor: DMDAGetInfo failed.");

	// Set the monitor to save performance plots (has to be in parallel)
	if (flagPerf) {
		// Only the master process will create the plot
		if (procId == 0) {
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
		}

		// monitorPerf will be called at each timestep
		ierr = TSMonitorSet(ts, monitorPerf, NULL, NULL);
		checkPetscError(ierr, "setupPetsc2DMonitor: TSMonitorSet (monitorPerf) failed.");
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
		checkPetscError(ierr, "setupPetsc2DMonitor: TSMonitorSet (computeHeliumFluence) failed.");

		// computeHeliumRetention2D will be called at each timestep
		ierr = TSMonitorSet(ts, computeHeliumRetention2D, NULL, NULL);
		checkPetscError(ierr, "setupPetsc2DMonitor: TSMonitorSet (computeHeliumRetention2D) failed.");

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
		checkPetscError(ierr, "setupPetsc2DMonitor: PetscOptionsGetInt (-start_stop) failed.");
		if (!flag)
			hdf5Stride2D = 1;

		// Initialize the HDF5 file for all the processes
		xolotlCore::HDF5Utils::initializeFile(hdf5OutputName2D, networkSize);

		// Get the solver handler
		auto solverHandler = PetscSolver::getSolverHandler();

		// Get the physical grid in the x direction
		auto grid = solverHandler->getXGrid();

		// Setup step size variables
		double hy = solverHandler->getStepSizeY();

		// Save the header in the HDF5 file
		xolotlCore::HDF5Utils::fillHeader(Mx, grid[1] - grid[0], My, hy);

		// Save the network in the HDF5 file
		xolotlCore::HDF5Utils::fillNetwork(network);

		// Finalize the HDF5 file
		xolotlCore::HDF5Utils::finalizeFile();

		// startStop2D will be called at each timestep
		ierr = TSMonitorSet(ts, startStop2D, NULL, NULL);
		checkPetscError(ierr, "setupPetsc2DMonitor: TSMonitorSet (startStop2D) failed.");
	}

	// Set the monitor to save surface plots of clusters concentration
	if (flag2DPlot) {
		// Only the master process will create the plot
		if (procId == 0) {
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
		}

		// monitorSurface2D will be called at each timestep
		ierr = TSMonitorSet(ts, monitorSurface2D, NULL, NULL);
		checkPetscError(ierr, "setupPetsc2DMonitor: TSMonitorSet (monitorSurface2D) failed.");
	}

	// Set the monitor to compute the helium concentrations
	if (flagConc) {
		// computeHeliumConc2D will be called at each timestep
		ierr = TSMonitorSet(ts, computeHeliumConc2D, NULL, NULL);
		checkPetscError(ierr, "setupPetsc2DMonitor: TSMonitorSet (computeHeliumConc2D) failed.");
	}

	// If the user wants the surface to be able to move
	if (solverHandler->moveSurface()) {
		// Initialize nInterstitial2D and previousIFlux2D before monitoring the
		// interstitial flux
		for (int j = 0; j < My; j++) {
			nInterstitial2D.push_back(0.0);
			previousIFlux2D.push_back(0.0);
		}

		// Set the monitor on the outgoing flux of interstitials at the surface
		// monitorInterstitial2D will be called at each timestep
		ierr = TSMonitorSet(ts, monitorInterstitial2D, NULL, NULL);
		checkPetscError(ierr, "setupPetsc2DMonitor: TSMonitorSet (monitorInterstitial2D) failed.");
	}

	// Set the monitor to simply change the previous time to the new time
	// monitorTime will be called at each timestep
	ierr = TSMonitorSet(ts, monitorTime, NULL, NULL);
	checkPetscError(ierr, "setupPetsc2DMonitor: TSMonitorSet (monitorTime) failed.");

	PetscFunctionReturn(0);
}

}

/* end namespace xolotlSolver */
