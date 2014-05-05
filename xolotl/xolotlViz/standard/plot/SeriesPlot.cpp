// Includes
#include "SeriesPlot.h"
#include "eavl.h"
#include "eavlDataSet.h"
#include "eavlColor.h"
#include "eavlRenderSurfaceOSMesa.h"
#include "eavlScene.h"
#include "eavl1DWindow.h"
#include <iostream>

using namespace xolotlViz;

#define W_WIDTH 1024
#define W_HEIGHT 1024

SeriesPlot::SeriesPlot(std::string name) : Plot(name) {
	plotDataProviders = std::make_shared< std::vector< std::shared_ptr<DataProvider> > > ();
}

SeriesPlot::~SeriesPlot() {
}

void SeriesPlot::render() {

	// Check if the label provider is set
	if (!plotLabelProvider){
		std::cout << "The LabelProvider is not set!!" << std::endl;
		return;
	}

	// Check if the data provider is set
	if (plotDataProviders->empty()){
		std::cout << "No DataProvider!!" << std::endl;
		return;
	}

    // Create an offscreen render surface
    eavlRenderSurface *surface = new eavlRenderSurfaceOSMesa;

    // Pick a background color
    eavlColor bg(1,1,1,1);

    eavlColor lineColor[6] =
    {eavlColor::magenta, eavlColor::blue, eavlColor::red,
    eavlColor::cyan, eavlColor::green, eavlColor::yellow};

    // Create a window
    eavlScene *scene = new eavl1DGLScene();
    eavl1DWindow *window = new eavl1DWindow(bg, surface, scene);
    window->Initialize();
    window->Resize(W_WIDTH,W_HEIGHT);

    // Set the log scale
    if (enableLogScale) window->view.view2d.logy = true;

    // Loop on all the data providers to plot the different series
    for (int i = 0; i < getDataProviderNumber(); i++){

    	// Get the value that will be plotted on X and Y
    	auto xVector = plotDataProviders->at(i)->getAxis1Vector();
    	auto yVector = plotDataProviders->at(i)->getAxis2Vector();

    	// Create the eavlDataSet
    	eavlDataSet *data = new eavlDataSet();
    	data->SetNumPoints(xVector.size());

    	// Give it the xVector
    	std::vector< std::vector<double> > coords;
    	coords.push_back(xVector);
    	std::vector<std::string> coordNames;
    	coordNames.push_back("xcoord");
    	AddRectilinearMesh(data, coords, coordNames, true, "RectilinearGridCells");

    	// Give the yVector to the axisValues
    	eavlArray *axisValues = new eavlFloatArray("coords", 1);
    	axisValues->SetNumberOfTuples(data->GetNumPoints());
    	for (int i = 0; i < yVector.size(); i++){
    		axisValues->SetComponentFromDouble(i, 0, yVector.at(i));
    	}

    	// Add the axisValues to a field of the data set
    	eavlField *field = new eavlField(1, axisValues, eavlField::ASSOC_POINTS);
    	data->AddField(field);

    	// Set up a plot for the data set
    	eavlRenderer *plot;
    	plot = new eavlCurveRenderer(data, NULL,
    			lineColor[i%6],
    			"",
    			"coords");

    	// Add the plot to the scene
    	scene->plots.push_back(plot);
    }

    // Set the view
    scene->ResetView(window);

    // Paint
    window->Paint();

    // Save the final buffer as an image
    char fn[25];
    sprintf(fn, (plotLabelProvider->titleLabel).c_str());
    window->SaveWindowAsPNM(fn);

	return;
}

void SeriesPlot::addDataProvider(std::shared_ptr<DataProvider> dataProvider){
	plotDataProviders->push_back(dataProvider);
	return;
}

std::shared_ptr<DataProvider> SeriesPlot::getDataProvider(int i) const {
	return plotDataProviders->at(i);
}

int SeriesPlot::getDataProviderNumber() const {
	return plotDataProviders->size();
}
