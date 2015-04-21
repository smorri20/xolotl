// Includes
#include "VideoPlot.h"

using namespace xolotlViz;

VideoPlot::VideoPlot(const std::string& name) : Plot(name) {
}

VideoPlot::~VideoPlot() {
}

void VideoPlot::render(const std::string& fileName) {
}

void VideoPlot::setFrameRate(double fRate) {
	frameRate = fRate;
	return;
}

double VideoPlot::getFrameRate() const {
	return frameRate;
}
