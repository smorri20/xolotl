#include "WFitFluxHandler.h"

using namespace xolotlSolver;

WFitFluxHandler::WFitFluxHandler(){

}

WFitFluxHandler::~WFitFluxHandler() {

}

double WFitFluxHandler::getIncidentFlux(std::vector<int> compositionVec,
				std::vector<double> position, double currentTime){

	double incidentFlux = 0.0e-16;

	incidentFlux = 0.0006 * position[0] * position[0] * position[0] - 0.0087 * position[0] * position[0] + 0.0300 * position[0];

	return incidentFlux;
}

void WFitFluxHandler::setOutgoingFlux(std::vector<int> compositionVec,
			std::vector<int> position, double time, double outgoingFlux){

}
