// Includes
#include "ScatterPlot.h"

using namespace xolotlViz;

ScatterPlot::ScatterPlot() {
}

ScatterPlot::~ScatterPlot() {
}

void ScatterPlot::render() {
}

std::vector< std::vector<double> > ScatterPlot::getPoints() const {
	//TODO Auto-generated method stub
}

std::string ScatterPlot::getAxis1Label() const {
	return plotLabelProvider.axis1Label;
}

std::string ScatterPlot::getAxis2Label() const {
	return plotLabelProvider.axis2Label;
}
