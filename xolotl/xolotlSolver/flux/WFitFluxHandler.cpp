#include "WFitFluxHandler.h"
#include <iostream>

using namespace xolotlSolver;

WFitFluxHandler::WFitFluxHandler() : stepSize(0.0e-16) {

}

void WFitFluxHandler::initializeFluxHandler(int numGridpoints, double step){

	// Set the step size
	stepSize = step;

	for(int i = 0; i < numGridpoints; i++)
	{
		auto x = i * stepSize;
		auto incidentFlux = 0.0006 * x * x * x - 0.0087 * x * x + 0.0300 * x;
		incidentFluxVec.push_back(incidentFlux);
	}
//	std::cout << "\n\nincidentFluxVec: " << std::endl;
//	for(int i = 0; i < numGridpoints; i++)
//		std::cout << incidentFluxVec[i] << std::endl;

}

double WFitFluxHandler::getIncidentFlux(std::vector<int> compositionVec,
				std::vector<double> position, double currentTime){

	auto i = position[0] / stepSize;

	return incidentFluxVec[i];
}

void WFitFluxHandler::setOutgoingFlux(std::vector<int> compositionVec,
			std::vector<int> position, double time, double outgoingFlux){

}
