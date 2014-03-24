// Includes
#include "Plot.h"

using namespace xolotlViz;

Plot::Plot() {
}

Plot::~Plot() {
}

std::string Plot::getUnit() {
	return " ";
}

std::string Plot::getTitle() {
	return " ";
}

void Plot::showLegend(bool legendShow) {
}

std::string Plot::getLegend() {
	return " ";
}

void Plot::write(std::string fileName) {
}

void Plot::setTitle(std::string title) {
}

void Plot::setPlottingStyle(PlottingStyle style) {
}

void Plot::setUnit(std::string unit) {
}

void Plot::setDataProvider(DataProvider dataProvider) {
}

DataProvider Plot::getDataProvider() {
	//TODO Auto-generated method stub
}

void Plot::setLabelProvider(LabelProvider labelProvider) {
}

LabelProvider Plot::getLabelProvider() {
	//TODO Auto-generated method stub
}
