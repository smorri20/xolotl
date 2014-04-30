#include "FeFitFluxHandler.h"

using namespace xolotlSolver;

FeFitFluxHandler::FeFitFluxHandler(){

}

FeFitFluxHandler::~FeFitFluxHandler() {

}

double FeFitFluxHandler::getIncidentFlux(std::vector<int> compositionVec,
				std::vector<double> position, double currentTime){

	return 0.0;
}

void FeFitFluxHandler::setOutgoingFlux(std::vector<int> compositionVec,
			std::vector<int> position, double time, double outgoingFlux){

}
