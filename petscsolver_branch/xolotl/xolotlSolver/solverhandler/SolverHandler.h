#ifndef SOLVERHANDLER_H
#define SOLVERHANDLER_H

// Includes
#include "ISolverHandler.h"

namespace xolotlSolver {

/**
 * This class and its subclasses realize the ISolverHandler interface to solve the
 * advection-diffusion-reaction problem with currently supported solvers.
 */
class SolverHandler: public ISolverHandler {
protected:

	//! The name of the network file
	std::string networkName;

	//! The original network created from the network loader.
	xolotlCore::PSIClusterReactionNetwork *network;

	//! Vector storing the grid in the x direction
	std::vector<double> grid;

	//! The grid step size in the y direction.
	double hY;

	//! The grid step size in the z direction.
	double hZ;

	//! The initial vacancy concentration.
	double initialVConc;

	//! The original flux handler created.
	xolotlCore::IFluxHandler *fluxHandler;

	//! The original temperature handler created.
	xolotlCore::ITemperatureHandler *temperatureHandler;

	//! The original diffusion handler created.
	xolotlCore::IDiffusionHandler *diffusionHandler;

	//! The original advection handler created.
	xolotlCore::IAdvectionHandler *advectionHandler;

	//! The original modified trap-mutation handler created.
	xolotlCore::ITrapMutationHandler *mutationHandler;

	//! The number of dimensions for the problem.
	int dimension;

	//! If the user wants to use a regular grid.
	bool useRegularGrid;

	//! Method generating the grid in the x direction
	void generateGrid(int nx, double hx) {
		// Clear the grid
		grid.clear();

		// Check if the user wants a regular grid
		if (useRegularGrid) {
			// The grid will me made of nx points separated by hx nm
			for (int l = 0; l < nx; l++){
				grid.push_back((double) l * hx);
			}
		}
		// If it is not regular do a fine mesh close to the surface and
		// increase the step size when away from the surface
		else {
			// Initialize the value of the previous point
			double previousPoint = 0.0;
			// The first grid point will be at x = 0.0
			grid.push_back(0.0);

			// The loop starts at 1 because the first grid point was
			// already added to the grid vector
			for (int l = 1; l < nx; l++) {
				// 0.1nm step near the surface (x < 2.5nm)
				if (l < 26) {
					grid.push_back(previousPoint + 0.1);
					previousPoint += 0.1;
				}
				// Then 0.25nm (2.5nm < x < 5.0nm)
				else if (l < 36) {
					grid.push_back(previousPoint + 0.25);
					previousPoint += 0.25;
				}
				// Then 0.5nm (5.0nm < x < 7.5nm)
				else if (l < 41) {
					grid.push_back(previousPoint + 0.5);
					previousPoint += 0.5;
				}
				// 1.0nm step size for all the other ones
				// (7.5nm < x)
				else {
					grid.push_back(previousPoint + 1.0);
					previousPoint += 1.0;
				}
			}
		}

		return;
	}

public:

	/**
	 * Initialize all the physics handlers that are needed to solve the ADR equations.
     * \see ISolverHandler.h
	 */
	void initializeHandlers(std::shared_ptr<xolotlFactory::IMaterialFactory> material,
			std::shared_ptr<xolotlCore::ITemperatureHandler> tempHandler,
			xolotlCore::Options &options) {

		// Set the flux handler
		fluxHandler = (xolotlCore::IFluxHandler *) material->getFluxHandler().get();

		// Set the temperature handler
		temperatureHandler = (xolotlCore::ITemperatureHandler *) tempHandler.get();

		// Set the diffusion handler
		diffusionHandler = (xolotlCore::IDiffusionHandler *) material->getDiffusionHandler().get();

		// Set the advection handler
		advectionHandler = (xolotlCore::IAdvectionHandler *) material->getAdvectionHandler().get();

		// Set the bubble bursting handler
		mutationHandler = (xolotlCore::ITrapMutationHandler *) material->getTrapMutationHandler().get();

		// Set the initial vacancy concentration
		initialVConc = options.getInitialVConcentration();

		// Set the number of dimension
		dimension = options.getDimensionNumber();

		// Look at if the user wants to use a regular grid in the x direction
		useRegularGrid = options.useRegularXGrid();

		return;
	}

	/**
	 * Initialize the network and network file name.
	 * \see ISolverHandler.h
	 */
	void initializeNetwork(const std::string& fileName,
			xolotlCore::PSIClusterReactionNetwork *net) {

		// Set the network loader
		networkName = fileName;

		// Set the network
		network = net;

		return;
	}

	/**
	 * Get the grid in the x direction.
	 * \see ISolverHandler.h
	 */
	std::vector<double> getXGrid() const {return grid;}

	/**
	 * Get the step size in the y direction.
	 * \see ISolverHandler.h
	 */
	double getStepSizeY() const {return hY;}

	/**
	 * Get the step size in the z direction.
	 * \see ISolverHandler.h
	 */
	double getStepSizeZ() const {return hZ;}

	/**
	 * Get the number of dimensions of the problem.
	 * \see ISolverHandler.h
	 */
	int getDimension() const {return dimension;}

	/**
	 * Get the flux handler.
	 * \see ISolverHandler.h
	 */
	xolotlCore::IFluxHandler *getFluxHandler() const {return fluxHandler;}

	/**
	 * Get the network.
	 * \see ISolverHandler.h
	 */
	xolotlCore::PSIClusterReactionNetwork *getNetwork() const {return network;}

}; //end class SolverHandler

} /* end namespace xolotlSolver */
#endif
