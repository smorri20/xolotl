#ifndef SOLVERHANDLERFACTORY_H
#define SOLVERHANDLERFACTORY_H

#include <memory>
#include <ISolverHandler.h>
#include <Options.h>

namespace xolotlFactory {


class SolverHandlerFactory {
private:
    // Our singleton solver handler.
    std::unique_ptr<xolotlSolver::ISolverHandler> handler;

public:
    /**
     * Construct a SolverHandlerFactory.
     * Default and copy constructors explicitly disallowed.
     */
    SolverHandlerFactory(void) = delete;
    SolverHandlerFactory(const SolverHandlerFactory& other) = delete;


    /**
     * Construct a SolverHandlerFactory.
     * Builds the desired type of solver.
     * Throws an exception if the solver handler can't be created 
     * based on the given options.
     *
     * @param options Options for the program.
     * @param network The reaction network we will use when solving.
     */
    SolverHandlerFactory(const xolotlCore::Options &options,
                            xolotlCore::IReactionNetwork& network);

    /**
     * Access the created solver handler.
     *
     * @return The solver handler.
     */
    xolotlSolver::ISolverHandler& getSolverHandler() { return *handler; }
    const xolotlSolver::ISolverHandler& getSolverHandler() const { return *handler; }
};

} // namespace xolotlFactory

#endif /* SOLVERHANDLERFACTORY_H */
