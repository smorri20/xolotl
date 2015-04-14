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

	// Compute the helium quantity to have at the interface
	double interHelium = (he - heliumReference) / 2.0;

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
	double quantity = (double) std::distance(gridPointList.begin(), gridPointList.end())
	* heliumReference;

	// Loop on all the interfaces and add their helium quantities
	for (std::forward_list<Interface *>::iterator it = interfaceList.begin();
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

	return;
}
