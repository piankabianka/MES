// BiankaDudaProjektMES.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include "Node.h"
#include "Element.h"
#include "GlobalData.h"
#include "Grid.h"
#include "UniversalElement.h"
using namespace std;


int main()
{
	//dzialanie wezla
	Node wezel;
	wezel.BC = true;
	cout << wezel.x << "\t" << wezel.y << "\t" << wezel.BC << endl;

	//dzialanie elementu
	Element element;
	element.nodesIDtab[3] = 2;
	cout << element.elementID << "\t" << element.nodesIDtab[3] << endl;

	//dzialanie wczytywania
	GlobalData data;
	data.showData();

	//dzialanie siatki
	Grid grid;
	
	grid.showGrid();
	grid.showDetailedGrid();

	//dzialanie elementu uniwersalnego

	UniversalElement elementU;
	elementU.showShapeFunctionsMatrix();
	elementU.showKsiDerivativeMatrix();
	elementU.showEtaDerivativeMatrix();

		 


	system("pause");
    return 0;
}

