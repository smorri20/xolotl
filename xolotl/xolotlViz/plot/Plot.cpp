// Includes
#include "Plot.h"

using namespace xolotlViz;

Plot::Plot() {
}

Plot::~Plot() {
}

void Plot::write(std::string fileName) {
	render();
	return;
}

void Plot::setPlottingStyle(PlottingStyle style) {
	plotStyle = style;
	return;
}

PlottingStyle Plot::getPlottingStyle() {
	return plotStyle;
}

void Plot::setDataProvider(std::shared_ptr<DataProvider> dataProvider) {
	plotDataProvider = dataProvider;
	return;
}

std::shared_ptr<DataProvider> Plot::getDataProvider() const {
	return plotDataProvider;
}

void Plot::setLabelProvider(std::shared_ptr<LabelProvider> labelProvider) {
	plotLabelProvider = labelProvider;
	return;
}

std::shared_ptr<LabelProvider> Plot::getLabelProvider() const {
	return plotLabelProvider;
}

void Plot::showLegend(bool legendShow) {
	enableLegend = legendShow;
	return;
}

std::string Plot::getLegend() const {
	return " ";
}

void Plot::setLogScale(bool logScale) {
	enableLogScale = logScale;
	return;
}
