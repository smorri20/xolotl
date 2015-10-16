// Includes
#include "Bubble.h"
#include <Constants.h>

using namespace xolotlCore;

Bubble::Bubble(int point, double he, double size) {
	// Clear everything
	gridPointList.clear();
	interfaceList.clear();

	// Set the helium reference value
	heliumReference = 4.0 * size * tungstenDensity;

	// Add the point to the list of points
	gridPointList.emplace_front(point);

	// Initialize the helium quantity to share at the interface
	double interHelium = 0.0;

	// Check if the helium quantity is bigger than the reference
	if (he > heliumReference) {
		// We are in the case where the bubble is already saturated
		// The helium will be split between the heliumQuantity and the interfaces
		heliumQuantity = heliumReference;
		// Compute the helium quantity to have at the interface
		interHelium = (he - heliumReference) / 2.0;
	}
	else {
		// There won't be any helium on the interface yet
		heliumQuantity = he;
	}

	// Add the corresponding interfaces (two in 1D)
	// Left one
	Interface * leftInter = new Interface(point, point - 1, interHelium);
	interfaceList.emplace_front(leftInter);
	// Right one
	Interface * rightInter = new Interface(point, point + 1, interHelium);
	interfaceList.emplace_front(rightInter);

	return;
}

double Bubble::getHeliumQuantity() {
	// Initialize the helium quantity
	double quantity = heliumQuantity;

	// Loop on all the interfaces and add their helium quantities
	for (auto it = interfaceList.begin();
			it != interfaceList.end(); it++) {
		quantity += (*it)->heliumQuantity;
	}

	return quantity;
}

std::forward_list<int> * Bubble::getGridPointList() {
	return &gridPointList;
}

std::forward_list<Bubble::Interface *> * Bubble::getInterfaceList() {
	return &interfaceList;
}

void Bubble::addGridPoint(int point) {
	gridPointList.emplace_front(point);
	// Add the helium reference value to the helium quantity
	heliumQuantity += heliumReference;

	return;
}

void Bubble::addHelium(double he, double dt, double size) {
	// he is only accounting for the external flux here,
	// we have to add the flux from diffusion through the interface to it
	// Loop on the interfaces
	for (std::forward_list<Interface *>::iterator inter = interfaceList.begin();
			inter != interfaceList.end(); inter++) {
		he += (*inter)->heliumFlux * dt * size;
	}

	// Check the size of the bubble
	int bubbleSize = std::distance(gridPointList.begin(), gridPointList.end());
	if (bubbleSize == 1) {
		// Check if the bubble is saturated yet
		if (heliumQuantity < heliumReference) {
			// The bubble is not saturated, we can add the helium to it
			heliumQuantity += he;

			// Check if it is saturated now
			if (heliumQuantity > heliumReference) {
				// It is saturated and will have to spread the excess on the interfaces
				he = heliumQuantity - heliumReference;
				// The quantity cannot be more than the reference
				heliumQuantity = heliumReference;
			}
			else {
				// There isn't any leftover helium to add to the interfaces
				he = 0.0;
			}
		}
	}

	// Divide the helium quantity on all the interfaces
	int nInter = std::distance(interfaceList.begin(), interfaceList.end());
	he = he / (double) nInter;
	// Loop on the interfaces
	for (std::forward_list<Interface *>::iterator inter = interfaceList.begin();
			inter != interfaceList.end(); inter++) {
		// Add the helium quantity
		(*inter)->heliumQuantity += he;
	}

	return;
}

void Bubble::resetHeliumQuantity(double he) {
	// Get the number of point covered by the bubble
	int nPoints = std::distance(gridPointList.begin(), gridPointList.end());

	// Check if the bubble is saturated
	if (he > (double) nPoints * heliumReference) {
		// The bubble is saturated, there will be helium on the interfaces
		heliumQuantity = (double) nPoints * heliumReference;
		// Keep the leftover for interfaces
		he = he - heliumQuantity;
		// Divide the helium quantity on all the interfaces
		int nInter = std::distance(interfaceList.begin(), interfaceList.end());
		he = he / (double) nInter;
		// Loop on the interfaces
		for (std::forward_list<Interface *>::iterator inter = interfaceList.begin();
				inter != interfaceList.end(); inter++) {
			// Add the helium quantity
			(*inter)->heliumQuantity = he;
		}
	}
	// Else the bubble is not saturated, there is no helium on the interface
	else {
		heliumQuantity = he;
		// Loop on the interfaces
		for (std::forward_list<Interface *>::iterator inter = interfaceList.begin();
				inter != interfaceList.end(); inter++) {
			// Add the helium quantity
			(*inter)->heliumQuantity = 0.0;
		}
	}

	return;
}
