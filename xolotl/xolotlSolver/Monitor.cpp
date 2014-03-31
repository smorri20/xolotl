// Includes
#include "../xolotlViz/plot/Plot.h"
#include "../xolotlViz/dataprovider/Point.h"
#include "../xolotlViz/plot/ScatterPlot.h"
#include "../xolotlViz/dataprovider/CvsXDataProvider.h"
#include "../xolotlViz/labelprovider/LabelProvider.h"
#include "PetscSolver.h"
#include <petscts.h>
#include <petscsys.h>
#include <sstream>
#include <vector>
#include <memory>

using namespace xolotlSolver;

/* ----- Error Handling Code ----- */

/**
 * This operation checks a Petsc error code and converts it to a bool.
 * @param errorCode The Petsc error code.
 * @return True if everything is OK, false otherwise.
 */
static inline bool checkPetscError(PetscErrorCode errorCode) {
	CHKERRQ(errorCode);
}

//! The pointer to the plot that will be used to visualize the data.
std::shared_ptr<xolotlViz::Plot> plot;

/**
 * This is a monitoring method that the user has to change to plot the data he/she
 * wants to monitor.
 */
static PetscErrorCode monitorSolve(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	// Network size
	const int size = PetscSolver::getNetwork()->size();
	// The array of cluster names
	std::vector<std::string> names(size);
	// Get the processor id
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);
	// Header and output string streams
	std::stringstream header, outputData;
	// Create a stream for naming the file
	std::stringstream outputFileNameStream;
	outputFileNameStream << "xolotl_out_" << procId << "_" << timestep;
	PetscErrorCode ierr;
	PetscViewer viewer;
	PetscReal *solutionArray, *gridPointSolution, x, hx;
	PetscInt xs, xm, Mx;
	int xi, i;

	PetscFunctionBeginUser;

	// Create the PETScViewer and get the data
	VecGetArray(solution, &solutionArray);
	PetscViewerASCIIOpen(PETSC_COMM_WORLD, outputFileNameStream.str().c_str(),
			&viewer);

	// Create the header for the file
	auto reactants = PetscSolver::getNetwork()->getAll();
	std::shared_ptr<PSICluster> cluster;
	header << "# t x ";
	for (int i = 0; i < size; i++) {
		// Get the cluster from the list, its id and composition
		cluster = std::dynamic_pointer_cast<PSICluster>(reactants->at(i));
		int id = cluster->getId() - 1;
		auto composition = cluster->getComposition();
		// Make the header entry
		std::stringstream name;
		name << (cluster->getName()).c_str() << "_(" << composition["He"] << ","
				<< composition["V"] << "," << composition["I"] << ") ";
		// Push the header entry on the array
		name >> names[id];
	}
	for (int i = 0; i < size; i++) {
		header << names[i] << " ";
	}
	header << "\n";
	PetscViewerASCIIPrintf(viewer, header.str().c_str());

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE);
	checkPetscError(ierr);
	// Setup some step size variables
	hx = 8.0 / (PetscReal) (Mx - 1);
	checkPetscError(ierr);

	// Create a Point vector to store the data to give to the data provider
	// for the visualization
	std::shared_ptr< std::vector<xolotlViz::Point> > myPoints(
			new (std::vector<xolotlViz::Point>));

	// Print the solution data
	for (xi = xs; xi < xs + xm; xi++) {
		// Dump x
		x = xi * hx;
		outputData << timestep << " " << x << " ";
		// Get the pointer to the beginning of the solution data for this grid point
		gridPointSolution = solutionArray + size * xi;
		// Update the concentrations in the network to have physics results
		// (non negative)
		PetscSolver::getNetwork()->updateConcentrationsFromArray(gridPointSolution);
		// Get the concentrations from the network
		double concentrations[size];
		double * concentration = &concentrations[0];
		PetscSolver::getNetwork()->fillConcentrationsArray(concentration);
		// Dump the data to the stream
		for (i = 0; i < size; i++) {
			outputData << concentration[i] << " ";
		}

		// Create a Point and add it to myPoints
		xolotlViz::Point aPoint;
		aPoint.value = concentration[2];
		aPoint.t = time; aPoint.x = x;
		myPoints->push_back(aPoint);

		// End the line
		outputData << "\n";
	}

	// Get the data provider
	auto dataProvider = plot->getDataProvider();

	// Give it the points
	dataProvider->setPoints(myPoints);
//	plot->setDataProvider(dataProvider);

	// Change the title of the plot
	std::stringstream title;
	title << "Concentration at t = " << time;
	plot->plotLabelProvider->titleLabel = title.str();

	// Render
	plot->render();

	// Dump the data to file
	PetscViewerASCIIPrintf(viewer, outputData.str().c_str());
	// Restore the array and kill the viewer
	VecRestoreArray(solution, &solutionArray);
	PetscViewerDestroy(&viewer);

	PetscFunctionReturn(0);
}

/**
 * This operation sets up a monitor that will call monitorSolve
 * @param ts The time stepper
 * @return A standard PETSc error code
 */
PetscErrorCode setupPetscMonitor(TS ts) {
	PetscErrorCode ierr;
	PetscBool flg;

	PetscFunctionBeginUser;
	ierr = PetscOptionsHasName(NULL, "-mymonitor", &flg);
	checkPetscError(ierr);
	if (!flg)
		PetscFunctionReturn(0);

	// Create a ScatterPlot
	plot = std::shared_ptr<xolotlViz::ScatterPlot> (
			new xolotlViz::ScatterPlot());

	// Create and set the label provider
	std::shared_ptr<xolotlViz::LabelProvider> labelProvider(
			new xolotlViz::LabelProvider());
	labelProvider->axis1Label = "x Position on the Grid";
	labelProvider->axis2Label = "Concentration";

	// Give it to the plot
	plot->setLabelProvider(labelProvider);

	// Create the data provider
	std::shared_ptr<xolotlViz::CvsXDataProvider> dataProvider(
			new xolotlViz::CvsXDataProvider());

	// Give it to the plot
	plot->setDataProvider(dataProvider);

	// monitorSolve will be called at each timestep
	ierr = TSMonitorSet(ts, monitorSolve, NULL, NULL);
	checkPetscError(ierr);
	PetscFunctionReturn(0);
}
