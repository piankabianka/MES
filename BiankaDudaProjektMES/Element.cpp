#include "stdafx.h"
#include "Node.h"
#include "Element.h"
#include "UniversalElement.h"
#include <iostream>
using namespace std;

Element::Element() {
	elementID = 0;
	nodesIDtab[0] = 0;
	nodesIDtab[1] = 0;
	nodesIDtab[2] = 0;
	nodesIDtab[3] = 0;
}
