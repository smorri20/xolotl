#include "HDF5Utils.h"
#include <iostream>
#include <sstream>
#include <hdf5.h>

using namespace xolotlCore;

hid_t fileId, concGroupId, concSId, networkGroupId, networkSId, headerGroupId;
herr_t status;

void HDF5Utils::initializeFile(int timeStep, int networkSize) {
	// Set the name of the file
	std::stringstream fileName;
	fileName << "xolotlStop_" << timeStep << ".h5";

	// Create the file
	fileId = H5Fcreate(fileName.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
			H5P_DEFAULT);

	// Create the group where the header will be stored
	headerGroupId = H5Gcreate2(fileId, "headerGroup", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Create the group where the network will be stored
	networkGroupId = H5Gcreate2(fileId, "networkGroup", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Create the dataspace for the network with dimension dims
	hsize_t dims[2];
	dims[0] = networkSize;
	dims[1] = 8;
	networkSId = H5Screate_simple(2, dims, NULL);

	// Create the group where the concentrations will be stored
	concGroupId = H5Gcreate2(fileId, "concentrationsGroup", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Create the dataspace for the concentrations with dimension dim
	hsize_t dim[1];
	dim[0] = networkSize;
	concSId = H5Screate_simple(1, dim, NULL);

	return;
}

void HDF5Utils::fillHeader(int physicalDim, int refinement, double time, double deltaTime) {
	// Create, write, and close the physicalDim attribute
	hid_t dimSId = H5Screate(H5S_SCALAR);
	hid_t dimAId = H5Acreate2 (headerGroupId, "physicalDim", H5T_STD_I32LE, dimSId,
	                             H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(dimAId, H5T_STD_I32LE, &physicalDim);
	status = H5Aclose(dimAId);

	// Create, write, and close the refinement attribute
	hid_t refineSId = H5Screate(H5S_SCALAR);
	hid_t refineAId = H5Acreate2 (headerGroupId, "refinement", H5T_STD_I32LE, refineSId,
	                             H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(refineAId, H5T_STD_I32LE, &refinement);
	status = H5Aclose(refineAId);

	// Create, write, and close the absolute time attribute
	hid_t timeSId = H5Screate(H5S_SCALAR);
	hid_t timeAId = H5Acreate2 (headerGroupId, "absoluteTime", H5T_IEEE_F64LE, timeSId,
	                             H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(timeAId, H5T_IEEE_F64LE, &time);
	status = H5Aclose(timeAId);

	// Create, write, and close the timestep time attribute
	hid_t deltaSId = H5Screate(H5S_SCALAR);
	hid_t deltaAId = H5Acreate2 (headerGroupId, "deltaTime", H5T_IEEE_F64LE, deltaSId,
	                             H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(deltaAId, H5T_IEEE_F64LE, &deltaTime);
	status = H5Aclose(deltaAId);

	return;
}

void HDF5Utils::fillNetwork(std::shared_ptr<PSIClusterReactionNetwork> network) {
	// Create the array that will store the network
	int networkSize = network->size();
	double networkArray [networkSize][8];

	// Get all the reactants
	auto reactants = network->getAll();

	// Loop on them
	for (int i = 0; i < networkSize; i++) {
		// Get the i-th reactant
		std::shared_ptr<PSICluster> reactant =
				std::static_pointer_cast<PSICluster>(reactants->at(i));

		// Get the reactant Id to keep the same order as the input file
		int id = reactant->getId() - 1;

		// Get its composition to store it
		auto composition = reactant->getComposition();
		networkArray[id][0] = composition["He"];
		networkArray[id][1] = composition["V"];
		networkArray[id][2] = composition["I"];

		// Get its binding energies to store them
		auto bindingEnergies = reactant->getBindingEnergies();
		networkArray[id][3] = bindingEnergies.at(0); // Helium binding energy
		networkArray[id][4] = bindingEnergies.at(1); // Vacancy binding energy
		networkArray[id][5] = bindingEnergies.at(2); // Interstitial binding energy

		// Get its migration energy to store it
		double migrationEnergy = reactant->getMigrationEnergy();
		networkArray[id][6] = migrationEnergy;

		// Get its diffusion factor to store it
		double diffusionFactor = reactant->getDiffusionFactor();
		networkArray[id][7] = diffusionFactor;
	}

	// Create the dataset for the network
	hid_t datasetId = H5Dcreate2(networkGroupId, "network", H5T_IEEE_F64LE, networkSId,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write networkArray in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, &networkArray);

	// Create the attribute for the network size
	hid_t networkSizeSId = H5Screate(H5S_SCALAR);
	hid_t networkSizeAId = H5Acreate2 (datasetId, "networkSize", H5T_STD_I32LE, networkSizeSId,
	                             H5P_DEFAULT, H5P_DEFAULT);

	// Write it
	status = H5Awrite(networkSizeAId, H5T_STD_I32LE, &networkSize);

	// Close everything
	status = H5Aclose(networkSizeAId);
	status = H5Dclose(datasetId);

	return;
}

void HDF5Utils::fillConcentrations(double * concArray, int index, double position) {
	// Set the dataset name
	std::stringstream datasetName;
	datasetName << "position_" << index;

	// Create the dataset of concentrations for this position
	hid_t datasetId = H5Dcreate2(concGroupId, datasetName.str().c_str(), H5T_IEEE_F64LE, concSId,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write concArray in the dataset
	status = H5Dwrite(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
	H5P_DEFAULT, concArray);

	// Create the attribute for the physical position
	hid_t posSId = H5Screate(H5S_SCALAR);
	hid_t posAId = H5Acreate2 (datasetId, "physicalPos", H5T_IEEE_F64LE, posSId,
	                             H5P_DEFAULT, H5P_DEFAULT);

	// Write it
	status = H5Awrite(posAId, H5T_IEEE_F64LE, &position);

	// Close everything
	status = H5Aclose(posAId);
	status = H5Dclose(datasetId);

	return;
}

void HDF5Utils::finalizeFile() {
	// Close everything
	status = H5Gclose(headerGroupId);
	status = H5Gclose(networkGroupId);
	status = H5Gclose(concGroupId);
	status = H5Fclose(fileId);

	return;
}

void HDF5Utils::readHeader(std::string fileName, int & physicalDim, double & time, double & deltaTime) {
	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

	// Open the header group
	hid_t groupId = H5Gopen(fileId, "/headerGroup", H5P_DEFAULT);

	// Open and read the physicalDim attribute
	hid_t physicalDimAId = H5Aopen(groupId, "physicalDim", H5P_DEFAULT);
	status = H5Aread(physicalDimAId, H5T_STD_I32LE, &physicalDim);
	status = H5Aclose(physicalDimAId);

	// Open and read the absoluteTime attribute
	hid_t timeAId = H5Aopen(groupId, "absoluteTime", H5P_DEFAULT);
	status = H5Aread(timeAId, H5T_IEEE_F64LE, &time);
	status = H5Aclose(timeAId);

	// Open and read the deltaTime attribute
	hid_t deltaTimeAId = H5Aopen(groupId, "deltaTime", H5P_DEFAULT);
	status = H5Aread(deltaTimeAId, H5T_IEEE_F64LE, &deltaTime);
	status = H5Aclose(deltaTimeAId);

	// Close everything
	status = H5Gclose(groupId);
	status = H5Fclose(fileId);

	return;
}

std::vector< std::vector <double> > HDF5Utils::readNetwork(std::string fileName) {
	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

	// Open the dataset
	hid_t datasetId = H5Dopen(fileId, "/networkGroup/network", H5P_DEFAULT);

	// Open and read the networkSize attribute
	hid_t networkSizeAId = H5Aopen(datasetId, "networkSize", H5P_DEFAULT);
	int networkSize = 0;
	status = H5Aread(networkSizeAId, H5T_STD_I32LE, &networkSize);
	status = H5Aclose(networkSizeAId);

	// Create the array that will receive the network
	double networkArray[networkSize][8];

	// Read the data set
	status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &networkArray);

	// Fill the vector to return with the dataset
	std::vector< std::vector <double> > networkVector;
	// Loop on the size of the network
	for (int i = 0; i < networkSize; i++) {
		// Create the line to give to the vector
		std::vector <double> line;
		for (int j = 0; j < 8; j++) {
			line.push_back(networkArray[i][j]);
		}
		networkVector.push_back(line);
	}

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Fclose(fileId);

	return networkVector;
}

void HDF5Utils::readGridPoint(std::string fileName, int networkSize, int i, double * concentrations) {
	// Open the given HDF5 file with read only access
	fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

	// Set the dataset name
	std::stringstream datasetName;
	datasetName << "concentrationsGroup/position_" << i;

	// Open the dataset
	hid_t datasetId = H5Dopen(fileId, datasetName.str().c_str(), H5P_DEFAULT);

	// Create the array that will receive the network
	double conc[networkSize];

	// Read the data set
	status = H5Dread(datasetId, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &conc);

	// Loop on the size of the network
	for (int i = 0; i < networkSize; i++) {
		concentrations[i] = conc[i];
	}

	// Close everything
	status = H5Dclose(datasetId);
	status = H5Fclose(fileId);

	return;
}
