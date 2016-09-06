// Includes
#include <PetscSolverHandler.h>
#include <HDF5Utils.h>
#include <MathUtils.h>
#include <Constants.h>

namespace xolotlSolver {

void PetscSolverHandler::getDiagonalFill(PetscInt *diagFill,
		int diagFillSize) {
	// Get all the super clusters
	auto superClusters = network->getAll("Super");

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network->size() + 2 * superClusters.size();
	const int diagSize = dof * dof;

	// Declarations for the loop
	std::vector<int> connectivity;
	int connectivityLength, id, index;

	// Fill the diagonal block if the sizes match up
	if (diagFillSize == diagSize) {
		// Get the connectivity for each reactant
		for (int i = 0; i < dof - 2 * superClusters.size(); i++) {
			// Get the reactant and its connectivity
			auto reactant = allReactants->at(i);
			connectivity = reactant->getConnectivity();
			connectivityLength = connectivity.size();
			// Get the reactant id so that the connectivity can be lined up in
			// the proper column
			id = reactant->getId() - 1;
			// Create the vector that will be inserted into the dFill map
			std::vector<int> columnIds;
			// Add it to the diagonal fill block
			for (int j = 0; j < connectivityLength; j++) {
				// The id starts at j*connectivity length and is always offset
				// by the id, which denotes the exact column.
				index = id * dof + j;
				diagFill[index] = connectivity[j];
				// Add a column id if the connectivity is equal to 1.
				if (connectivity[j] == 1) {
					columnIds.push_back(j);
				}
			}
			// Update the map
			dFillMap[id] = columnIds;
		}
		// Get the connectivity for each moment
		for (int i = 0; i < superClusters.size(); i++) {
			// Get the reactant and its connectivity
			auto reactant = superClusters[i];
			connectivity = reactant->getConnectivity();
			connectivityLength = connectivity.size();
			// Get the helium momentum id so that the connectivity can be lined up in
			// the proper column
			id = reactant->getHeMomentumId() - 1;

			// Create the vector that will be inserted into the dFill map
			std::vector<int> columnIds;
			// Add it to the diagonal fill block
			for (int j = 0; j < connectivityLength; j++) {
				// The id starts at j*connectivity length and is always offset
				// by the id, which denotes the exact column.
				index = (id) * dof + j;
				diagFill[index] = connectivity[j];
				// Add a column id if the connectivity is equal to 1.
				if (connectivity[j] == 1) {
					columnIds.push_back(j);
				}
			}
			// Update the map
			dFillMap[id] = columnIds;

			// Get the vacancy momentum id so that the connectivity can be lined up in
			// the proper column
			id = reactant->getVMomentumId() - 1;

			// Add it to the diagonal fill block
			for (int j = 0; j < connectivityLength; j++) {
				// The id starts at j*connectivity length and is always offset
				// by the id, which denotes the exact column.
				index = (id) * dof + j;
				diagFill[index] = connectivity[j];
			}
			// Update the map
			dFillMap[id] = columnIds;
		}

		// Debug output
//		std::cout << "Number of degrees of freedom = " << dof
//				<< std::endl;
//		printf("\n");
//		for (i = 0; i < dof; i++) {
//			for (j = 0; j < dof; j++) {
//				printf("%d ", dfill[i * dof + j]);
//			}
//			printf("\n");
//		}
//		printf("\n");
	} else {
		std::string err =
				"PetscSolverHandler Exception: Invalid diagonal block size!\n";
		throw std::string(err);
	}

	return;
}

} /* end namespace xolotlSolver */
