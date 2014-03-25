// Includes
#include "Plot.h"

using namespace xolotlViz;

Plot::Plot() {
}

Plot::~Plot() {
}

std::string Plot::getUnit() const {
	return plotUnit;
}

std::string Plot::getTitle() const {
	return plotTitle;
}

void Plot::showLegend(bool legendShow) {
	enableLegend = legendShow;
	return;
}

std::string Plot::getLegend() const {
	return " ";
}

void Plot::write(std::string fileName) {
}

void Plot::setTitle(std::string title) {
	plotTitle = title;
	return;
}

void Plot::setPlottingStyle(PlottingStyle style) {
	plotStyle = style;
	return;
}

void Plot::setUnit(std::string unit) {
	plotUnit = unit;
	return;
}

void Plot::setDataProvider(DataProvider dataProvider) {
}

DataProvider Plot::getDataProvider() const {
	//TODO Auto-generated method stub
}

void Plot::setLabelProvider(LabelProvider labelProvider) {
}

LabelProvider Plot::getLabelProvider() const {
	//TODO Auto-generated method stub
}
