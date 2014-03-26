// Includes
#include "VideoPlot.h"

using namespace xolotlViz;

VideoPlot::VideoPlot() {
}

VideoPlot::~VideoPlot() {
}

void VideoPlot::render() {
}

void VideoPlot::setFrameRate(double fRate) {
	frameRate = fRate;
	return;
}

double VideoPlot::getFrameRate() const {
	return frameRate;
}
