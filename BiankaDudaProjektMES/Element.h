#pragma once
#include "stdafx.h"
#include "Node.h"
#include "UniversalElement.h"
#include <iostream>
using namespace std;

struct Element {
	int elementID;
	int nodesIDtab[4];

	double jMatrix[4][4];
	double detJMatrix[4];
	double jacobianMatrix[4][4];

	void calculateJacobianMatrixes(UniversalElement elementU);
	Element();
};

Element::Element() {
	elementID = 0;
	nodesIDtab[0] =  0;
	nodesIDtab[1] = 0;
	nodesIDtab[2] = 0;
	nodesIDtab[3] = 0;
}
