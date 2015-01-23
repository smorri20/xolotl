// Includes
#include <PetscSolver2DHandler.h>
#include <HDF5Utils.h>
#include <MathUtils.h>
#include <Constants.h>

namespace xolotlSolver {

/**
 * This operation checks a Petsc error code and converts it to a bool.
 *
 * @param errorCode The Petsc error code.
 * @return True if everything is OK, false otherwise.
 */
inline bool checkPetscError(PetscErrorCode errorCode) {
	CHKERRQ(errorCode);
}

void PetscSolver2DHandler::createSolverContext(DM &da, int nx, double hx, int ny,
		double hy, int nz, double hz) {
	PetscErrorCode ierr;

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network->size();

	// Initialize the all reactants pointer
	allReactants = network->getAll();

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create distributed array (DMDA) to manage parallel grid and vectors
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED,
	DMDA_STENCIL_STAR, nx, ny, PETSC_DECIDE, PETSC_DECIDE, dof, 1, NULL, NULL, &da);
	checkPetscError(ierr);

	// Set the step size
	h = hx;

	// Set the size of the partial derivatives vectors
	clusterPartials.resize(dof, 0.0);
	reactingPartialsForCluster.resize(dof, 0.0);

	// Set the last temperature to 0
	lastTemperature = 0.0;

	/*  The only spatial coupling in the Jacobian is due to diffusion.
	 *  The ofill (thought of as a dof by dof 2d (row-oriented) array represents
	 *  the nonzero coupling between degrees of freedom at one point with degrees
	 *  of freedom on the adjacent point to the left or right. A 1 at i,j in the
	 *  ofill array indicates that the degree of freedom i at a point is coupled
	 *  to degree of freedom j at the adjacent point.
	 *  In this case ofill has only a few diagonal entries since the only spatial
	 *  coupling is regular diffusion.
	 */
	PetscInt *ofill, *dfill;
	ierr = PetscMalloc(dof * dof * sizeof(PetscInt), &ofill);
	checkPetscError(ierr);
	ierr = PetscMalloc(dof * dof * sizeof(PetscInt), &dfill);
	checkPetscError(ierr);
	ierr = PetscMemzero(ofill, dof * dof * sizeof(PetscInt));
	checkPetscError(ierr);
	ierr = PetscMemzero(dfill, dof * dof * sizeof(PetscInt));
	checkPetscError(ierr);

	// Fill ofill, the matrix of "off-diagonal" elements that represents diffusion
	diffusionHandler->initializeOFill(network, ofill);

	// Get the diagonal fill
	getDiagonalFill(dfill, dof * dof);

	// Load up the block fills
	ierr = DMDASetBlockFills(da, dfill, ofill);
	checkPetscError(ierr);

	// Free the temporary fill arrays
	ierr = PetscFree(ofill);
	checkPetscError(ierr);
	ierr = PetscFree(dfill);
	checkPetscError(ierr);

	return;
}

void PetscSolver2DHandler::initializeConcentration(DM &da, Vec &C) const {
	PetscErrorCode ierr;

	// Pointer for the concentration vector
	PetscScalar ***concentrations;
	ierr = DMDAVecGetArrayDOF(da, C, &concentrations);
	checkPetscError(ierr);

	// Get the local boundaries
	PetscInt xs, xm, ys, ym;
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	checkPetscError(ierr);

	// Get the last time step written in the HDF5 file
	int tempTimeStep = -2;
	bool hasConcentrations = xolotlCore::HDF5Utils::hasConcentrationGroup(
			networkName, tempTimeStep);

	// Get the total size of the grid for the boundary conditions
	PetscInt Mx, My;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	checkPetscError(ierr);

	// Initialize the flux handler
	fluxHandler->initializeFluxHandler(Mx, h, h);

	// Initialize the advection handler
	advectionHandler->initialize(network);

	// Pointer for the concentration vector at a specific grid point
	PetscScalar *concOffset;

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network->size();

	// Get the single vacancy ID.
	int vacancyIndex = (network->get(xolotlCore::vType, 1)->getId()) - 1;

	// Loop on all the grid points
	for (int j = ys; j < ys + ym; j++) {
		for (int i = xs; i < xs + xm; i++) {
			concOffset = concentrations[j][i];

			// Loop on all the clusters to initialize at 0.0
			for (int n = 0; n < dof; n++) {
				concOffset[n] = 0.0;
			}

			// Initialize the vacancy concentration
			if (i > 0 && i < Mx - 1 && j > 0 && j < My - 1) {
				concOffset[vacancyIndex] = initialVConc / h;
			}
		}
	}

	// If the concentration must be set from the HDF5 file
	if (hasConcentrations) {
		// Loop on the full grid
		for (int j = 0; j < My; j++) {
			for (int i = 0; i < Mx; i++) {
				// Read the concentrations from the HDF5 file
				auto concVector = xolotlCore::HDF5Utils::readGridPoint(networkName,
						tempTimeStep, i, j);

				// Change the concentration only if we are on the locally owned part of the grid
				if (i >= xs && i < xs + xm && j >= ys && j < ys + ym) {
					concOffset = concentrations[j][i];
					// Loop on the concVector size
					for (int l = 0; l < concVector.size(); l++) {
						concOffset[(int) concVector.at(l).at(0)] =
								concVector.at(l).at(1);
					}
				}
			}
		}
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOF(da, C, &concentrations);
	checkPetscError(ierr);

	return;
}

void PetscSolver2DHandler::updateConcentration(TS &ts, Vec &localC, Vec &F,
		PetscReal ftime, bool &temperatureChanged) {
	PetscErrorCode ierr;

	// Get the local data vector from petsc
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);

	// Get the total size of the grid for the boundary conditions
	PetscInt Mx, My;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	checkPetscError(ierr);

	// Pointers to the Petsc arrays that start at the beginning (xs) of the
	// local array!
	PetscScalar ***concs, ***updatedConcs;
	// Get pointers to vector data
	ierr = DMDAVecGetArrayDOF(da, localC, &concs);
	checkPetscError(ierr);
	ierr = DMDAVecGetArrayDOF(da, F, &updatedConcs);
	checkPetscError(ierr);

	//Get local grid boundaries
	PetscInt xs, xm, ys, ym;
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	checkPetscError(ierr);

	// The following pointers are set to the first position in the conc or
	// updatedConc arrays that correspond to the beginning of the data for the
	// current grid point. They are accessed just like regular arrays.
	PetscScalar *concOffset, *updatedConcOffset;

	// Set some step size variable
	double s = 1.0 / (h * h);

	// Get the incident flux vector
	auto incidentFluxVector = fluxHandler->getIncidentFluxVec(ftime);

	// Declarations for variables used in the loop
	int reactantIndex;
	double flux;
	auto heCluster = (xolotlCore::PSICluster *) network->get(xolotlCore::heType, 1);
	xolotlCore::PSICluster *cluster = NULL;
	double **concVector = new double*[5];
	std::vector<double> gridPosition = { 0.0, 0.0, 0.0 };

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network->size();

	// Loop over grid points computing ODE terms for each grid point
	for (int yj = ys; yj < ys + ym; yj++) {
		for (int xi = xs; xi < xs + xm; xi++) {
			//xi = 1; // Uncomment this line for debugging in a single cell.

			// Compute the old and new array offsets
			concOffset = concs[yj][xi];
			updatedConcOffset = updatedConcs[yj][xi];

			// Fill the concVector with the pointer to the middle, left, right, bottom, and top grid points
			concVector[0] = concOffset; // middle
			concVector[1] = concs[yj][xi - 1]; // left
			concVector[2] = concs[yj][xi + 1]; // right
			concVector[3] = concs[yj - 1][xi]; // bottom
			concVector[4] = concs[yj + 1][xi]; // top

			// Boundary conditions
			if (xi == 0 || xi == Mx - 1 || yj == 0 || yj == My - 1) {
				for (int i = 0; i < dof; i++) {
					updatedConcOffset[i] = 1.0 * concOffset[i];
				}

				continue;
			}

			double x = xi * h;
			double y = yj * h;

			// Vector representing the physical position
			gridPosition[0] = x;
			gridPosition[1] = y;
			auto temperature = temperatureHandler->getTemperature(gridPosition,
					ftime);

			// Update the network if the temperature changed
			if (!xolotlCore::equal(temperature, lastTemperature)) {
				network->setTemperature(temperature);
				lastTemperature = temperature;
				temperatureChanged = true;
			}

			// Copy data into the PSIClusterReactionNetwork so that it can
			// compute the fluxes properly. The network is only used to compute the
			// fluxes and hold the state data from the last time step. I'm reusing
			// it because it cuts down on memory significantly (about 400MB per
			// grid point) at the expense of being a little tricky to comprehend.
			network->updateConcentrationsFromArray(concOffset);

			// ----- Account for flux of incoming He by computing forcing that
			// produces He of cluster size 1 -----
			if (heCluster) {
				reactantIndex = heCluster->getId() - 1;
				// Update the concentration of the cluster
				updatedConcOffset[reactantIndex] += incidentFluxVector[xi];
			}

			// ---- Compute diffusion over the locally owned part of the grid -----
			diffusionHandler->computeDiffusion(network, s, concVector,
					updatedConcOffset);

			// ---- Compute advection over the locally owned part of the grid -----

			// Change the second value in the concVector because the advection handler
			// excepts the middle point for the first value and the right point for the
			// second one
			concVector[1] = concs[yj][xi + 1]; // right

			advectionHandler->computeAdvection(network, h, gridPosition, concVector,
					updatedConcOffset);

			// ----- Compute all of the new fluxes -----
			for (int i = 0; i < dof; i++) {
				cluster = (xolotlCore::PSICluster *) allReactants->at(i);
				// Compute the flux
				flux = cluster->getTotalFlux();
				// Update the concentration of the cluster
				reactantIndex = cluster->getId() - 1;
				updatedConcOffset[reactantIndex] += flux;
			}

			// Uncomment this line for debugging in a single cell.
			//break;
		}
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOF(da, localC, &concs);
	checkPetscError(ierr);
	ierr = DMDAVecRestoreArrayDOF(da, F, &updatedConcs);
	checkPetscError(ierr);
	ierr = DMRestoreLocalVector(da, &localC);
	checkPetscError(ierr);

	return;
}

void PetscSolver2DHandler::computeOffDiagonalJacobian(TS &ts, Vec &localC, Mat &J) const {
	PetscErrorCode ierr;

	// Get the distributed array
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);

	// Get the total size of the grid for the boundary conditions
	PetscInt Mx, My;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	checkPetscError(ierr);

	// Setup some step size variables
	double s = 1.0 / (h * h);

	// Get pointers to vector data
	PetscScalar ***concs;
	ierr = DMDAVecGetArrayDOF(da, localC, &concs);
	checkPetscError(ierr);

	// Get local grid boundaries
	PetscInt xs, xm, ys, ym;
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	checkPetscError(ierr);

	// The degree of freedom is the size of the network
	const int dof = network->size();

	// Pointer to the concentrations at a given grid point
	PetscScalar *concOffset;

	// Get the total number of diffusing clusters
	const int nDiff = diffusionHandler->getNumberOfDiffusing();

	// Get the total number of advecting clusters
	const int nAdvec = advectionHandler->getNumberOfAdvecting();

	// Arguments for MatSetValuesStencil called below
	MatStencil row, cols[5];
	PetscScalar vals[5 * nDiff];
	PetscInt indices[nDiff];
	std::vector<double> gridPosition = { 0.0, 0.0, 0.0 };

	/*
	 Loop over grid points computing Jacobian terms for diffusion and advection
	 at each grid point
	 */
	for (int yj = ys; yj < ys + ym; yj++) {
		for (int xi = xs; xi < xs + xm; xi++) {
			//xi = 1; // Uncomment this line for debugging in a single cell

			// Boundary conditions
			if (xi == 0 || xi == Mx - 1 || yj == 0 || yj == My - 1) continue;

			// Set the grid position
			gridPosition[0] = xi * h;
			gridPosition[1] = yj * h;

			// Copy data into the PSIClusterReactionNetwork so that it can
			// compute the new concentrations.
			concOffset = concs[yj][xi];
			network->updateConcentrationsFromArray(concOffset);

			// Get the partial derivatives for the diffusion
			diffusionHandler->computePartialsForDiffusion(network, s, vals, indices);

			// Loop on the number of diffusion cluster to set the values in the Jacobian
			for (int i = 0; i < nDiff; i++) {
				// Set grid coordinate and component number for the row
				row.i = xi;
				row.j = yj;
				row.c = indices[i];

				// Set grid coordinates and component numbers for the columns
				// corresponding to the middle, left, right, bottom and top grid points
				cols[0].i = xi; // middle
				cols[0].j = yj;
				cols[0].c = indices[i];
				cols[1].i = xi - 1; // left
				cols[1].j = yj;
				cols[1].c = indices[i];
				cols[2].i = xi + 1; // right
				cols[2].j = yj;
				cols[2].c = indices[i];
				cols[3].i = xi; // bottom
				cols[3].j = yj - 1;
				cols[3].c = indices[i];
				cols[4].i = xi; // top
				cols[4].j = yj + 1;
				cols[4].c = indices[i];

				ierr = MatSetValuesStencil(J, 1, &row, 5, cols, vals + (5 * i), ADD_VALUES);
				checkPetscError(ierr);
			}

			// Get the partial derivatives for the advection
			advectionHandler->computePartialsForAdvection(network, h, vals,
					indices, gridPosition);

			// Loop on the number of advecting cluster to set the values in the Jacobian
			for (int i = 0; i < nAdvec; i++) {
				// Set grid coordinate and component number for the row
				row.i = xi;
				row.j = yj;
				row.c = indices[i];

				// Set grid coordinates and component numbers for the columns
				// corresponding to the left and middle grid points
				cols[0].i = xi - 1; // left
				cols[0].j = yj;
				cols[0].c = indices[i];
				cols[1].i = xi; // middle
				cols[1].j = yj;
				cols[1].c = indices[i];

				// Update the matrix
				ierr = MatSetValuesStencil(J, 1, &row, 2, cols, vals + (2 * i), ADD_VALUES);
				checkPetscError(ierr);
			}

			//break;   // Uncomment this line for debugging in a single cell.
		}
	}

	return;
}

void PetscSolver2DHandler::computeDiagonalJacobian(TS &ts, Vec &localC, Mat &J) {
	PetscErrorCode ierr;

	// Get the distributed array
	DM da;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);

	// Get the total size of the grid for the boundary conditions
	PetscInt Mx, My;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
	PETSC_IGNORE);
	checkPetscError(ierr);

	// Get pointers to vector data
	PetscScalar ***concs;
	ierr = DMDAVecGetArrayDOF(da, localC, &concs);
	checkPetscError(ierr);

	// Get local grid boundaries
	PetscInt xs, xm, ys, ym;
	ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);
	checkPetscError(ierr);

	// The degree of freedom is the size of the network
	const int dof = network->size();

	// Pointer to the concentrations at a given grid point
	PetscScalar *concOffset;

	// Arguments for MatSetValuesStencil called below
	MatStencil rowId;
	MatStencil colIds[dof];
	int pdColIdsVectorSize = 0;

	// Declarations for variables used in the loop
	int reactantIndex;

	// Loop over the grid points
	for (int yj = ys; yj < ys + ym; yj++) {
		for (int xi = xs; xi < xs + xm; xi++) {
			//xi = 1; // Uncomment this line for debugging in a single cell

			// Boundary conditions
			if (xi == 0 || xi == Mx - 1 || yj == 0 || yj == My - 1) continue;

			// Copy data into the PSIClusterReactionNetwork so that it can
			// compute the new concentrations.
			concOffset = concs[yj][xi];
			network->updateConcentrationsFromArray(concOffset);

			// Update the column in the Jacobian that represents each reactant
			for (int i = 0; i < dof; i++) {
				auto reactant = allReactants->at(i);
				// Get the reactant index
				reactantIndex = reactant->getId() - 1;

				// Set grid coordinate and component number for the row
				rowId.i = xi;
				rowId.j = yj;
				rowId.c = reactantIndex;

				// Get the partial derivatives
				reactant->getPartialDerivatives(clusterPartials);
				// Get the list of column ids from the map
				auto pdColIdsVector = dFillMap.at(reactantIndex);
				// Number of partial derivatives
				pdColIdsVectorSize = pdColIdsVector.size();
				// Loop over the list of column ids
				for (int j = 0; j < pdColIdsVectorSize; j++) {
					// Set grid coordinate and component number for a column in the list
					colIds[j].i = xi;
					colIds[j].j = yj;
					colIds[j].c = pdColIdsVector[j];
					// Get the partial derivative from the array of all of the partials
					reactingPartialsForCluster[j] =
							clusterPartials[pdColIdsVector[j]];
					// Reset the cluster partial value to zero. This is much faster
					// than using memset.
					clusterPartials[pdColIdsVector[j]] = 0.0;
				}
				// Update the matrix
				ierr = MatSetValuesStencil(J, 1, &rowId, pdColIdsVectorSize, colIds,
						reactingPartialsForCluster.data(), ADD_VALUES);
				checkPetscError(ierr);
			}

			// Uncomment this line for debugging in a single cell.
			//break;
		}
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArrayDOF(da, localC, &concs);
	checkPetscError(ierr);
	ierr = DMRestoreLocalVector(da, &localC);
	checkPetscError(ierr);

	return;
}

} /* end namespace xolotlSolver */
