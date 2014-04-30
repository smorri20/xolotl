#include "FeFitFluxHandler.h"
#include <iostream>

using namespace xolotlSolver;

FeFitFluxHandler::FeFitFluxHandler() : stepSize(0.0e-16) {

}

void FeFitFluxHandler::initializeFluxHandler(int numGridpoints, double step){

	// Set the step size
	stepSize = step;

	//	// Single He production rate = 0.0029593 * incidentFlux  ( in 1 / nm^3 /sec)
	//	// where incidentFlux is shown below by Single He flux
	//	//
	//	// Single He flux
	//	// calculated at each depth (x <= 128nm, if x > 128nm, no production)
	//	a0 =  -0.0007309d0
	//	a1 =   -0.002933d0
	//	b1 =    0.003981d0
	//	a2 =   -0.004196d0
	//	b2 =    0.008418d0
	//	a3 =   -0.001564d0
	//	b3 =    0.009943d0
	//	a4 =    0.002591d0
	//	b4 =    0.006301d0
	//	a5 =    0.003863d0
	//	b5 =     0.00088d0
	//	a6 =    0.002226d0
	//	b6 =   -0.001758d0
	//	a7 =   0.0005369d0
	//	b7 =   -0.001357d0
	//	a8 =  1.092d-5
	//	b8 =  -0.0003591d0
	//	w =     0.01388d0
	//
	//	incidentFlux = a0 + a1*cos(x*w) + b1*sin(x*w) + a2*cos(2.d0*x*w) + b2*sin(2.d0*x*w) +
	//				a3*cos(3.d0*x*w) + b3*sin(3.d0*x*w) + a4*cos(4.d0*x*w) + b4*sin(4.d0*x*w) +
	//				a5*cos(5.d0*x*w) + b5*sin(5.d0*x*w) + a6*cos(6.d0*x*w) + b6*sin(6.d0*x*w) +
	//				a7*cos(7.d0*x*w) + b7*sin(7.d0*x*w) + a8*cos(8.d0*x*w) + b8*sin(8.d0*x*w);

	for(int i = 0; i < numGridpoints; i++)
	{
		auto x = i * stepSize;
		auto incidentFlux = 0.0;
		incidentFluxVec.push_back(incidentFlux);
	}
//	std::cout << "\n\nincidentFluxVec: " << std::endl;
//	for(int i = 0; i < numGridpoints; i++)
//		std::cout << incidentFluxVec[i] << std::endl;

}

double FeFitFluxHandler::getIncidentFlux(std::vector<int> compositionVec,
				std::vector<double> position, double currentTime){

	auto i = position[0] / stepSize;

	return incidentFluxVec[i];
}

void FeFitFluxHandler::setOutgoingFlux(std::vector<int> compositionVec,
			std::vector<int> position, double time, double outgoingFlux){

}
