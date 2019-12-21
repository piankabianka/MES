#pragma once
#include "stdafx.h"
#include "Node.h"
#include <iostream>
using namespace std;

struct Element {
	int elementID;
	int nodesIDtab[4];
	Element();
};

Element::Element() {
	elementID = 0;
	nodesIDtab[0] =  0;
	nodesIDtab[1] = 0;
	nodesIDtab[2] = 0;
	nodesIDtab[3] = 0;
}
