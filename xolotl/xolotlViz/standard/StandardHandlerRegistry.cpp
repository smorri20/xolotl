#include "StandardHandlerRegistry.h"
#include <Plot.h>
#include <ScatterPlot.h>
#include <SeriesPlot.h>
#include <SurfacePlot.h>
#include <VideoPlot.h>
#include <DataProvider.h>
#include <LabelProvider.h>

namespace xolotlViz {

StandardHandlerRegistry::StandardHandlerRegistry(OutputType otype)
    : outputType(otype)
{
}

StandardHandlerRegistry::~StandardHandlerRegistry() {
}

std::shared_ptr<IPlot> StandardHandlerRegistry::getPlot(std::string name, PlotType type) {
    bool isRaster = (outputType==png);
//    std::cerr << "isRaster="<<isRaster<<std::endl;
	switch(type) {
        case PlotType::SCATTER: return std::make_shared <ScatterPlot> (name, isRaster);
        case PlotType::SERIES: return std::make_shared <SeriesPlot> (name, isRaster);
        case PlotType::SURFACE: return std::make_shared <SurfacePlot> (name, isRaster);
        case PlotType::VIDEO: return std::make_shared <VideoPlot> (name, isRaster);
        default: return std::make_shared <Plot> (name, isRaster);
	}
}

}    //end namespace xolotlViz

