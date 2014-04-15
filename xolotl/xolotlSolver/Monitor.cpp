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
#include "../xolotlPerf/HandlerRegistryFactory.h"

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

	PetscFunctionBeginUser;

	// Get the number of processes
	int size;
	MPI_Comm_size(PETSC_COMM_WORLD, &size);

	// Print a warning if only one process
	if (size == 1) {
		cout << "You are trying to plot things that don't have any sense!! "
				<< "\nRemove -mymonitor or run in parallel." << endl;
		PetscFunctionReturn(0);
	}

	// Get the current process ID
	int procId;
	MPI_Comm_rank(PETSC_COMM_WORLD, &procId);

	// Get the solve timer
	auto solverTimer = xolotlPerf::getHandlerRegistry()->getTimer("solve");

	// Stop it to access its value
	solverTimer->stop();

	// Master process
	if (procId == 0) {

		// Create a Point vector to store the data to give to the data provider
		// for the visualization
		auto myPoints = std::make_shared<std::vector<xolotlViz::Point> >();

		// Give it the value for procId = 0
		xolotlViz::Point aPoint;
		aPoint.value = solverTimer->getValue();
		aPoint.t = time;
		aPoint.x = procId;
		myPoints->push_back(aPoint);

		// Loop on all the other processes
		for (int i = 1; i < size; i++) {
			double counter = 0.;

			// Receive the value from the other processes
			MPI_Recv(&counter, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			// Give it the value for procId = i
			aPoint.value = counter;
			aPoint.t = time;
			aPoint.x = i;
			myPoints->push_back(aPoint);
		}

		// Get the data provider and give it the points
		plot->getDataProvider()->setPoints(myPoints);

		// Change the title of the plot
		std::stringstream title;
		title << "timeTS" << timestep << ".pnm";
		plot->plotLabelProvider->titleLabel = title.str();

		// Render
		plot->render();
	}

	else {
		double counter = solverTimer->getValue();

		// Send the value of the timer to the master process
		MPI_Send(&counter, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	// Restart the timer
	solverTimer->start();

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
	plot = std::make_shared<xolotlViz::ScatterPlot>();

	// Create and set the label provider
	auto labelProvider = std::make_shared<xolotlViz::LabelProvider>();
	labelProvider->axis1Label = "x Position on the Grid";
	labelProvider->axis2Label = "Concentration";

	// Give it to the plot
	plot->setLabelProvider(labelProvider);

	// Create the data provider
	auto dataProvider = std::make_shared<xolotlViz::CvsXDataProvider>();

	// Give it to the plot
	plot->setDataProvider(dataProvider);

	// monitorSolve will be called at each timestep
	ierr = TSMonitorSet(ts, monitorSolve, NULL, NULL);
	checkPetscError(ierr);
	PetscFunctionReturn(0);
}

}

/* end namespace xolotlSolver */
