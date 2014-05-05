// Includes
#include "SurfacePlot.h"
#include "eavl.h"
#include "eavlDataSet.h"
#include "eavlColor.h"
#include "eavlRenderSurfaceOSMesa.h"
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
	eavlArray *axisValues = new eavlFloatArray(plotDataProvider->dataName, 1);
	axisValues->SetNumberOfTuples(data->GetNumPoints());
	for (int i = 0; i < zVector.size(); i++){
		axisValues->SetComponentFromDouble(i, 0, zVector.at(i));
	}

	// Add the axisValues to a field of the data set
	eavlField *field = new eavlField(1, axisValues, eavlField::ASSOC_POINTS);
	data->AddField(field);

    // Create an offscreen render surface
    eavlRenderSurface *surface = new eavlRenderSurfaceOSMesa;

    // Pick a background color
//    eavlColor bg(0.15, 0.05, 0.1, 1.0);
    eavlColor bg(0.5, 0.5, 0.5, 0.5);

    // Create a 2D scene
    eavlScene *scene = new eavl2DGLScene();

    // Create the window
    eavl2DWindow *window = new eavl2DWindow(bg, surface, scene);
    window->Initialize();
    window->Resize(W_WIDTH,W_HEIGHT);

    // Print the title
    window->AddWindowAnnotation(plotLabelProvider->titleLabel, 0,1, .5,1, 0.08);

    // Print the axis
    window->AddViewportAnnotation(plotLabelProvider->axis1Label, 0,-1, 0,-.1, .5,1, 0.05);
    window->AddViewportAnnotation(plotLabelProvider->axis2Label, -1,0, -.1,0, .5,0, 0.05, 90);
    window->AddWindowAnnotation(plotLabelProvider->axis3Label, 0,1, .5,5, 0.05);

    // Set up a plot for the data set
    eavlRenderer *plot;
    plot = new eavlPseudocolorRenderer(data, NULL,
                                       "temperature",
                                       false,
                                       data->GetCellSet(0)->GetName(),
                                       plotDataProvider->dataName);
    scene->plots.push_back(plot);

    // Set the view
    scene->ResetView(window);

    // Paint
    window->Paint();

    // Save the final buffer as an image
    char fn[25];
    sprintf(fn, (fileName).c_str());
    window->SaveWindowAsPNM(fn);

	return;
}
