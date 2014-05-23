#ifndef HDF5UTILS_H
#define HDF5UTILS_H

#include <PSIClusterReactionNetwork.h>
#include <memory>

namespace xolotlCore {

namespace HDF5Utils {

	/**
	 * Create the HDF5 file with the needed structure.
	 * @param timeStep The number of the time step.
	 * @param networkSize The total number of cluster in the network.
	 * @param gridSize The total number of grid points.
	 */
	void initializeFile(int timeStep, int networkSize, int gridSize);

	/**
	 * Fill the header.
	 * @param physicalDim The physical length of the material on which one is solving the ADR equation.
	 * @param refinement The refinement of the grid.
	 * @param time The physical time at this time step.
	 * @param deltaTime The physical length of the time step.
	 */
	void fillHeader(int physicalDim, int refinement, double time, double deltaTime);

	/**
	 * Fill the network dataset.
	 * @param network The network of clusters.
	 */
	void fillNetwork(std::shared_ptr<PSIClusterReactionNetwork> network);

	/**
	 * Fill the concentration dataset at a specific grid point.
	 * @param concArray The vector of concentration at a grid point.
	 * @param position The physical position on the grid.
	 */
	void fillConcentrations(double * concArray, int index, double position);

	/**
	 * Add the data to the file and close it.
	 */
	void finalizeFile();

	/**
	 * Read the header of a HDF5 file.
	 * @param fileName The name of the file to read from.
	 * @param physicalDim The physical length of the material to be changed.
	 * @param time The physical time to be changed.
	 * @param deltaTime The time step length to be changed.
	 */
	void readHeader(std::string fileName, int & physicalDim, double & time, double & deltaTime);

	/**
	 * Read the network from a HDF5 file.
	 * @param fileName The name of the file to read from.
	 * @return The vector of vector which contain the network dataset.
	 */
	std::vector< std::vector <double> > readNetwork(std::string fileName);

	/**
	 * Read the i-th grid point concentrations from a HDF5 file.
	 * @param fileName The name of the file to read from.
	 * @param networkSize The size of the network.
	 * @param i The index of the grid point.
	 * @param concentrations The array of concentrations.
	 */
	void readGridPoint(std::string fileName, int networkSize, int i, double * concentrations);

};

} /* namespace xolotlCore */
#endif
