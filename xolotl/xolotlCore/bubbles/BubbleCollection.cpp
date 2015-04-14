// Includes
#include "BubbleCollection.h"

using namespace xolotlCore;

void BubbleCollection::addBubble(Bubble * bubble) {
	bubbleList.emplace_front(bubble);

	return;
}

std::forward_list<Bubble *> * BubbleCollection::getBubbleList() {
	return &bubbleList;
}

bool BubbleCollection::isPointABubble(int point) {
	// Loop on all the bubbles
	for (std::forward_list<Bubble *>::iterator it = bubbleList.begin();
			it != bubbleList.end(); it++) {
		// Get the list of grid points of the bubble
		auto gridPoints = (*it)->getGridPointList();
		// Loop on the grid points
		for (std::forward_list<int>::iterator pt = gridPoints->begin();
				pt != gridPoints->end(); pt++) {
			if ((*pt) == point) return true;
		}
	}

	// By default return false
	return false;
}

double BubbleCollection::getHeliumQuantity() {
	// Initialize the helium quantity
	double quantity = 0.0;

	// Loop on all the bubbles
	for (std::forward_list<Bubble *>::iterator it = bubbleList.begin();
				it != bubbleList.end(); it++) {
		// Add the helium quantity
		quantity += (*it)->getHeliumQuantity();
	}

	return quantity;
}

bool BubbleCollection::getNeighbor(Bubble * firstBubble, Bubble * secondBubble) {
	// Get the grid point list from both bubbles
	auto firstGridPointList = firstBubble->getGridPointList();
	auto secondGridPointList = secondBubble->getGridPointList();

	// Loop on the first list
	for (std::forward_list<int>::iterator firstPt = firstGridPointList->begin();
			firstPt != firstGridPointList->end(); firstPt++) {
		// Loop on the second list
		for (std::forward_list<int>::iterator secondPt = secondGridPointList->begin();
					secondPt != secondGridPointList->end(); secondPt++) {
			if ((*firstPt) == (*secondPt) - 1 || (*firstPt) == (*secondPt) + 1) {
				// They are next to each other
				return true;
			}
		}
	}

	// By default return false
	return false;
}
