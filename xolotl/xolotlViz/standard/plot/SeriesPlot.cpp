// Includes
#include "SeriesPlot.h"
#include "eavl.h"
#include "eavlDataSet.h"
#include "eavlColor.h"
#include "XolotlConfigViz.h"
#ifdef HAVE_OSMESA
#include <GL/gl_mangle.h>
#include "eavlRenderSurfaceOSMesa.h"
#include "eavlSceneRendererGL.h"
#include "eavlWorldAnnotatorGL.h"
#endif
#include "eavlRenderSurfacePS.h"
#include "eavlSceneRendererPS.h"
#include "eavlWorldAnnotatorPS.h"
#include "eavlScene.h"
#include "eavl1DWindow.h"
#include <iostream>

using namespace xolotlViz;

#define W_WIDTH 1024
#define W_HEIGHT 1024

SeriesPlot::SeriesPlot(std::string name) : Plot(name) {
	plotDataProviders = std::make_shared< std::vector< std::shared_ptr<IDataProvider> > > ();
}

SeriesPlot::~SeriesPlot() {
}

void SeriesPlot::render(std::string fileName) {

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

    // Pick a background color
    eavlColor bg(1,1,1,1);

    eavlColor lineColor[18] =
    {eavlColor::magenta, eavlColor::blue, eavlColor::red,
    eavlColor::cyan, eavlColor::green, eavlColor::yellow,
    eavlColor(0.5, 0.0, 1.0), eavlColor(1.0, 0.0, 0.5), eavlColor(0.5, 0.0, 0.5),
    eavlColor(0.0, 0.5, 0.5), eavlColor(1.0, 0.5, 0.0), eavlColor(0.5, 0.0, 0.0),
    eavlColor(0.0, 0.5, 0.0), eavlColor(0.0, 0.0, 0.5), eavlColor(0.5, 0.5, 0.5),
    eavlColor(0.2, 0.8, 0.2), eavlColor(0.8, 0.2, 0.2), eavlColor(0.2, 0.2, 0.1)};

    // Create a window
    eavlScene *scene = new eavl1DScene();
#ifdef HAVE_OSMESA
    ///\todo: get OpenGL mode set some proper way
    bool OpenGL_Mode = true;
#else
    bool OpenGL_Mode = false;
#endif
    eavlRenderSurface *surface;
    eavlSceneRenderer *renderer = NULL;
    eavlWorldAnnotator *annotator = NULL;
#ifdef HAVE_OSMESA
    if (OpenGL_Mode)
    {
        surface = new eavlRenderSurfaceOSMesa;
        renderer = new eavlSceneRendererGL;
        annotator = new eavlWorldAnnotatorGL;
    }
    else
#endif
    {
        surface = new eavlRenderSurfacePS;
        renderer = new eavlSceneRendererPS;
        annotator = new eavlWorldAnnotatorPS;
    }
    eavl1DWindow *window = new eavl1DWindow(bg, surface, scene, renderer, annotator);
    window->Initialize();
    window->Resize(W_WIDTH,W_HEIGHT);

    // Print the title
    auto titleAnnotation = new eavlScreenTextAnnotation(window, plotLabelProvider->titleLabel,
    		eavlColor::black, 0.065, 0.2, 0.8);
    titleAnnotation->SetAlignment(eavlTextAnnotation::HCenter, eavlTextAnnotation::VCenter);
    window->AddAnnotation(titleAnnotation);

    // Print the axis labels
    auto axis1Annotation = new eavlScreenTextAnnotation(window, plotLabelProvider->axis1Label,
    		eavlColor::black, 0.05, 0.2, -0.8);
    axis1Annotation->SetAlignment(eavlTextAnnotation::HCenter, eavlTextAnnotation::VCenter);
    window->AddAnnotation(axis1Annotation);
    auto axis2Annotation = new eavlScreenTextAnnotation(window, plotLabelProvider->axis2Label,
    		eavlColor::black, 0.05, -0.65, 0.0, 90.0);
    axis2Annotation->SetAlignment(eavlTextAnnotation::HCenter, eavlTextAnnotation::VCenter);
    window->AddAnnotation(axis2Annotation);

    // Add the time information
    auto timeAnnotation = new eavlScreenTextAnnotation(window, plotLabelProvider->timeLabel,
    		eavlColor::black, 0.055, -0.8, -0.85);
    timeAnnotation->SetAlignment(eavlTextAnnotation::HCenter, eavlTextAnnotation::VCenter);
    window->AddAnnotation(timeAnnotation);
    auto timeStepAnnotation = new eavlScreenTextAnnotation(window, plotLabelProvider->timeStepLabel,
    		eavlColor::black, 0.055, -0.8, -0.91);
    timeStepAnnotation->SetAlignment(eavlTextAnnotation::HCenter, eavlTextAnnotation::VCenter);
    window->AddAnnotation(timeStepAnnotation);

    // Set the log scale
    if (enableLogScale) window->view.view2d.logy = true;

    // Translate the viewport not to have the legend on top of the y axis
    window->view.vl = -0.45;
    window->view.vr = 0.85;

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
    	eavlArray *axisValues = new eavlFloatArray(plotDataProviders->at(i)->getDataName(), 1);
    	axisValues->SetNumberOfTuples(data->GetNumPoints());
    	for (int i = 0; i < yVector.size(); i++){
    		axisValues->SetComponentFromDouble(i, 0, yVector.at(i));
    	}

    	// Add the axisValues to a field of the data set
    	eavlField *field = new eavlField(1, axisValues, eavlField::ASSOC_POINTS);
    	data->AddField(field);

    	// Set up a plot for the data set
    	eavlPlot *plot;
        plot = new eavl1DPlot(data, "RectilinearGridCells");
        plot->SetSingleColor(lineColor[i%18]);
        plot->SetField(plotDataProviders->at(i)->getDataName());

    	// Add the plot to the scene
    	scene->plots.push_back(plot);
    }

    // Set the view
    scene->ResetView(window);

    // Paint
    window->Paint();

    // Save the final buffer as an image
    char fn[25];
    sprintf(fn, "%s", (fileName).c_str());
    ///\todo: file extension currently hard-coded by caller
    surface->SaveAs(fn, OpenGL_Mode ? eavlRenderSurface::PNG : eavlRenderSurface::EPS);

	return;
}

void SeriesPlot::addDataProvider(std::shared_ptr<IDataProvider> dataProvider){
	plotDataProviders->push_back(dataProvider);
	return;
}

std::shared_ptr<IDataProvider> SeriesPlot::getDataProvider(int i) const {
	return plotDataProviders->at(i);
}

int SeriesPlot::getDataProviderNumber() const {
	return plotDataProviders->size();
}
