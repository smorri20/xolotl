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
	PetscErrorCode ierr;
	PetscReal *solutionArray, *gridPointSolution, x, hx;
	PetscInt xs, xm, Mx;
	int xi, i;

	PetscFunctionBeginUser;

	// Fill the array of clusters name because the Id is not the same as
	// reactants->at(i)
	auto reactants = PetscSolver::getNetwork()->getAll();
	std::shared_ptr<PSICluster> cluster;
	for (int i = 0; i < size; i++) {

		// Get the cluster from the list, its id and composition
		cluster = std::dynamic_pointer_cast<PSICluster>(reactants->at(i));
		int id = cluster->getId() - 1;
		auto composition = cluster->getComposition();

		// Create the name
		std::stringstream name;
		name << (cluster->getName()).c_str() << "(" << composition["He"] << ","
				<< composition["V"] << "," << composition["I"] << ") ";

		// Push the header entry on the array
		name >> names[id];
	}

	// Get the data
	VecGetArray(solution, &solutionArray);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);

	// Get the corners of the grid
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr);
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE);
	checkPetscError(ierr);
	// Setup some step size variables
	hx = 8.0 / (PetscReal) (Mx - 1);

	// Create a Point vector to store the data to give to the data provider
	// for the visualization
	std::shared_ptr< std::vector<xolotlViz::Point> > myPoints(
			new (std::vector<xolotlViz::Point>));

	// Choice of the cluster to be plotted
	int iCluster = 2;

	// Print the solution data
	for (xi = xs; xi < xs + xm; xi++) {
		// Dump x
		x = xi * hx;
		// Get the pointer to the beginning of the solution data for this grid point
		gridPointSolution = solutionArray + size * xi;
		// Update the concentrations in the network to have physics results
		// (non negative)
		PetscSolver::getNetwork()->updateConcentrationsFromArray(gridPointSolution);
		// Get the concentrations from the network
		double concentrations[size];
		double * concentration = &concentrations[0];
		PetscSolver::getNetwork()->fillConcentrationsArray(concentration);

		// Create a Point with the concentration[iCluster] as the value
		// and add it to myPoints
		xolotlViz::Point aPoint;
		aPoint.value = concentration[iCluster];
		aPoint.t = time; aPoint.x = x;
		myPoints->push_back(aPoint);
	}

	// Get the data provider
	auto dataProvider = plot->getDataProvider();

	// Give it the points
	dataProvider->setPoints(myPoints);

	// Change the title of the plot
	std::stringstream title;
	title << names[iCluster] << "TS" << timestep << ".pnm";
	plot->plotLabelProvider->titleLabel = title.str();

	// Render
	plot->render();
	// Restore the array
	VecRestoreArray(solution, &solutionArray);

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

}

/* end namespace xolotlSolver */
