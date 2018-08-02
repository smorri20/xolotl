#include <iostream>
#include <cassert>
#include <SolverHandlerFactory.h>
#include <PetscSolver0DHandler.h>
#include <PetscSolver1DHandler.h>
#include <PetscSolver2DHandler.h>
#include <PetscSolver3DHandler.h>

namespace xolotlFactory {

// Create the desired type of handler registry.
SolverHandlerFactory::SolverHandlerFactory(const xolotlCore::Options &options,
		                                xolotlCore::IReactionNetwork& network) {

	// Get the wanted dimension
	int dim = options.getDimensionNumber();

	// Switch on the dimension
	// TODO Once we have widespread C++14 support, use std::make_unique
	// instead of this two-step construction.
	xolotlSolver::ISolverHandler* rawSolverHandler = nullptr;
	switch (dim) {
	case 0:
		rawSolverHandler = new xolotlSolver::PetscSolver0DHandler(network);
		break;
	case 1:
		rawSolverHandler = new xolotlSolver::PetscSolver1DHandler(network);
		break;
	case 2:
		rawSolverHandler = new xolotlSolver::PetscSolver2DHandler(network);
		break;
	case 3:
		rawSolverHandler = new xolotlSolver::PetscSolver3DHandler(network);
		break;
	default:
		// The requested dimension is not supported.
        std::ostringstream estr;
        estr << "Invalid dimension " << dim << " requested when creating solver handler";
        throw std::runtime_error(estr.str());
	}
	assert(rawSolverHandler != nullptr);
	handler = std::unique_ptr<xolotlSolver::ISolverHandler>(rawSolverHandler);
}

} // end namespace xolotlFactory

