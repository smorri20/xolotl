// Includes
#include "DataProvider.h"

using namespace xolotlViz;

DataProvider::DataProvider() {
}

DataProvider::~DataProvider() {
}

std::vector<Point> DataProvider::getDataPoints() const {
	return data;
}

void DataProvider::setPoints(std::vector<Point> points) {

	// Loop on all the points in the points vector
	for (auto it = points.begin();
			it != points.end(); it++) {

		// Add the current Point to the data vector
		data.push_back(*it);
	}
	return;
}

double DataProvider::getDataMean() const {
	// The size of the data vector
	int size = data.size();

	// Use to add the value of each Point
	double valueSum = 0.;

	// Loop on all the points in the data vector
	for (auto it = data.begin();
			it != data.end(); it++) {

		// Add the current value to the sum
		valueSum += (*it).value;
	}

	// Result equals the sum divided by the size
	double result = (double) valueSum / size;

	return result;
}
