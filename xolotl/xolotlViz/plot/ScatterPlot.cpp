// Includes
#include "ScatterPlot.h"
#include <iostream>

using namespace xolotlViz;

ScatterPlot::ScatterPlot() {
}

ScatterPlot::~ScatterPlot() {
}

void ScatterPlot::render() {
	if (!plotLabelProvider){
		std::cout << "The LabelProvider is not set!!" << std::endl;
		return;
	}
	if (!plotDataProvider){
		std::cout << "The DataProvider is not set!!" << std::endl;
		return;
	}

	std::cout << "The x axis label is: " << plotLabelProvider->axis1Label
			<< std::endl;

	std::cout << "The y axis label is: " << plotLabelProvider->axis2Label
			<< std::endl;

	std::cout << "The title is: " << plotLabelProvider->titleLabel
			<< std::endl;

	std::cout << "The unit is: " << plotLabelProvider->unitLabel
			<< std::endl;

	auto xVector = plotDataProvider->getAxis1Vector();
	auto yVector = plotDataProvider->getAxis2Vector();

	for (int i = 0; i < xVector.size(); i++){
		std::cout << xVector.at(i) << " " << yVector.at(i) << std::endl;
	}

	std::cout << std::endl;
	std::cout << std::endl;

	return;
}
