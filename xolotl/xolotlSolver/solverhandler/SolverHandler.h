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

	//! The vector of advection handlers.
	std::vector<xolotlCore::IAdvectionHandler *> advectionHandlers;

	//! The original modified trap-mutation handler created.
	xolotlCore::ITrapMutationHandler *mutationHandler;

	//! The original bubble bursting handler created.
	xolotlCore::IBubbleBurstingHandler *burstingHandler;

	//! The number of dimensions for the problem.
	int dimension;

	//! The portion of void at the beginning of the problem.
	double portion;

	//! If the user wants to use a regular grid.
	bool useRegularGrid;

	//! If the user wants to move the surface.
	bool movingSurface;

	//! Method generating the grid in the x direction
	void generateGrid(int nx, double hx, int surfacePos) {
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
			if (nx != 38) throw std::string("Wrong size of grid here!! ");
			else {
				grid = {0.0, 2.0, 4.0, 7.0, 12.0, 17.0, 22.0, 27.0, 32.0, 37.0,
						42.0, 47.0, 52.0, 57.0, 62.0, 67.0, 72.0, 77.0, 82.0, 87.0,
						92.0, 97.0, 102.0, 107.0, 112.0, 117.0, 128.0, 158.0, 188.0,
						218.0, 248.0, 278.0, 350.0, 420.0, 490.0, 496.0, 498.0, 500.0};
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

		// Set the advection handlers
		auto handlers = material->getAdvectionHandler();
		for (int i = 0; i < handlers.size(); i++) {
			advectionHandlers.push_back((xolotlCore::IAdvectionHandler *) handlers[i].get());
		}


		// Set the modified trap-mutation handler
		mutationHandler = (xolotlCore::ITrapMutationHandler *) material->getTrapMutationHandler().get();

		// Set the bubble bursting handler
		burstingHandler = (xolotlCore::IBubbleBurstingHandler *) material->getBubbleBurstingHandler().get();

		// Set the initial vacancy concentration
		initialVConc = options.getInitialVConcentration();

		// Set the number of dimension
		dimension = options.getDimensionNumber();

		// Set the void portion
		portion = options.getVoidPortion();

		// Look at if the user wants to use a regular grid in the x direction
		useRegularGrid = options.useRegularXGrid();

		// Should we be able to move the surface?
		auto map = options.getProcesses();
		movingSurface = map["movingSurface"];

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
	 * Get the initial vacancy concentration.
	 * \see ISolverHandler.h
	 */
	double getInitialVConc() const {return initialVConc;}

	/**
	 * To know if the surface should be able to move.
	 * \see ISolverHandler.h
	 */
	bool moveSurface() const {return movingSurface;}

	/**
	 * Get the flux handler.
	 * \see ISolverHandler.h
	 */
	xolotlCore::IFluxHandler *getFluxHandler() const {return fluxHandler;}

	/**
	 * Get the advection handler.
	 * \see ISolverHandler.h
	 */
	xolotlCore::IAdvectionHandler *getAdvectionHandler() const {return advectionHandlers[0];}

	/**
	 * Get the advection handlers.
	 * \see ISolverHandler.h
	 */
	std::vector<xolotlCore::IAdvectionHandler *> getAdvectionHandlers() const {return advectionHandlers;}

	/**
	 * Get the modified trap-mutation handler.
	 * \see ISolverHandler.h
	 */
	xolotlCore::ITrapMutationHandler *getMutationHandler() const {return mutationHandler;}

	/**
	 * Get the bubble bursting handler.
	 * \see ISolverHandler.h
	 */
	xolotlCore::IBubbleBurstingHandler *getBurstingHandler() const {return burstingHandler;}

	/**
	 * Get the network.
	 * \see ISolverHandler.h
	 */
	xolotlCore::PSIClusterReactionNetwork *getNetwork() const {return network;}

	/**
	 * Get the network name.
	 * \see ISolverHandler.h
	 */
	std::string getNetworkName() const {return networkName;}

}; //end class SolverHandler

} /* end namespace xolotlSolver */
#endif
