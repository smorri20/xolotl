// Includes
#include "SurfacePlot.h"
#include "eavl.h"
#include "eavlDataSet.h"
#include "eavlColor.h"
///\todo: Get HAVE_OSMESA set up properly in cmake.
#define HAVE_OSMESA
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
#include "eavl2DWindow.h"
#include <iostream>

using namespace xolotlViz;

#define W_WIDTH 1024
#define W_HEIGHT 1024

SurfacePlot::SurfacePlot(std::string name) : Plot(name) {
}

SurfacePlot::~SurfacePlot() {
}

void SurfacePlot::render(std::string fileName) {

	// Check if the label provider is set
	if (!plotLabelProvider){
		std::cout << "The LabelProvider is not set!!" << std::endl;
		return;
	}

	// Check if the data provider is set
	if (!plotDataProvider){
		std::cout << "The DataProvider is not set!!" << std::endl;
		return;
	}

	// Get the value that will be plotted on X, Y, and Z
	auto xVector = plotDataProvider->getAxis1Vector();
	auto yVector = plotDataProvider->getAxis2Vector();
	auto zVector = plotDataProvider->getAxis3Vector();

	// Create the eavlDataSet
    eavlDataSet *data = new eavlDataSet();
    data->SetNumPoints(zVector.size());

    // Give it the xVector and yVector
	std::vector< std::vector<double> > coords;
	coords.push_back(xVector);
	coords.push_back(yVector);
	std::vector<std::string> coordNames;
	coordNames.push_back("xcoord");
	coordNames.push_back("ycoord");
	AddRectilinearMesh(data, coords, coordNames, true, "RectilinearGridCells");

	// Give the zVector to the axisValues
	eavlArray *axisValues = new eavlFloatArray(plotDataProvider->getDataName(), 1);
	axisValues->SetNumberOfTuples(data->GetNumPoints());
	for (int i = 0; i < zVector.size(); i++){
		axisValues->SetComponentFromDouble(i, 0, zVector.at(i));
	}

	// Add the axisValues to a field of the data set
	eavlField *field = new eavlField(1, axisValues, eavlField::ASSOC_POINTS);
	data->AddField(field);

    // Pick a background color
//    eavlColor bg(0.15, 0.05, 0.1, 1.0);
    eavlColor bg(0.5, 0.5, 0.5, 0.5);

    // Create a window
    eavlScene *scene = new eavl2DScene();
    ///\todo: get OpenGL mode set some proper way
    bool OpenGL_Mode = false;
    eavlRenderSurface *surface;
    eavlSceneRenderer *renderer = NULL;
    eavlWorldAnnotator *annotator = NULL;
    if (OpenGL_Mode)
    {
        surface = new eavlRenderSurfaceOSMesa;
        renderer = new eavlSceneRendererGL;
        annotator = new eavlWorldAnnotatorGL;
    }
    else
    {
        surface = new eavlRenderSurfacePS;
        renderer = new eavlSceneRendererPS;
        annotator = new eavlWorldAnnotatorPS;
    }
    eavl2DWindow *window = new eavl2DWindow(bg, surface, scene, renderer, annotator);
    window->Initialize();
    window->Resize(W_WIDTH,W_HEIGHT);

    // Print the title
    auto titleAnnotation = new eavlScreenTextAnnotation(window, plotLabelProvider->titleLabel,
    		eavlColor::white, 0.065, 0.0, 0.96);
    titleAnnotation->SetAlignment(eavlTextAnnotation::HCenter, eavlTextAnnotation::VCenter);
    window->AddAnnotation(titleAnnotation);

    // Print the axis labels
    auto axis1Annotation = new eavlScreenTextAnnotation(window, plotLabelProvider->axis1Label,
    		eavlColor::white, 0.05, 0.0, -0.35);
    axis1Annotation->SetAlignment(eavlTextAnnotation::HCenter, eavlTextAnnotation::VCenter);
    window->AddAnnotation(axis1Annotation);
    auto axis2Annotation = new eavlScreenTextAnnotation(window, plotLabelProvider->axis2Label,
    		eavlColor::white, 0.05, -0.9, 0.0, 90.0);
    axis2Annotation->SetAlignment(eavlTextAnnotation::HCenter, eavlTextAnnotation::VCenter);
    window->AddAnnotation(axis2Annotation);
    auto axis3Annotation = new eavlScreenTextAnnotation(window, plotLabelProvider->axis3Label,
    		eavlColor::white, 0.05, 0.0, 0.78);
    axis3Annotation->SetAlignment(eavlTextAnnotation::HCenter, eavlTextAnnotation::VCenter);
    window->AddAnnotation(axis3Annotation);

    // Add the time information
    auto timeAnnotation = new eavlScreenTextAnnotation(window, plotLabelProvider->timeLabel,
    		eavlColor::white, 0.055, 0.8, -0.85);
    timeAnnotation->SetAlignment(eavlTextAnnotation::HCenter, eavlTextAnnotation::VCenter);
    window->AddAnnotation(timeAnnotation);
    auto timeStepAnnotation = new eavlScreenTextAnnotation(window, plotLabelProvider->timeStepLabel,
    		eavlColor::white, 0.055, 0.8, -0.91);
    timeStepAnnotation->SetAlignment(eavlTextAnnotation::HCenter, eavlTextAnnotation::VCenter);
    window->AddAnnotation(timeStepAnnotation);

    // Set up a plot for the data set
    eavlPlot *plot;
    plot = new eavlPlot(data, data->GetCellSet(0)->GetName());
    plot->SetField(plotDataProvider->getDataName());
    plot->SetColorTableName("temperature");
    scene->plots.push_back(plot);

    // Set the view
    scene->ResetView(window);

    // Paint
    window->Paint();

    // Save the final buffer as an image
    char fn[25];
    sprintf(fn, "%s", (fileName).c_str());
    ///\todo: file extension currently hard-coded by caller
    surface->SaveAs(fn, OpenGL_Mode ? eavlRenderSurface::PNM : eavlRenderSurface::EPS);

	return;
}
