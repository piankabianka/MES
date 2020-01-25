#pragma once
#include "stdafx.h"
#include "Node.h"
#include "Element.h"
#include "GlobalData.h"
#include <iostream>
using namespace std;
struct Grid {
	Node* nodesTab;
	Element* elementsTab;
	int nodesNumber;
	int elementNumber;
	UniversalElement universalElement;
	double jMatrix[4][4];
	double detJ[4];
	double jMatrixReverse[4][4];
	double dNdX[4][4];
	double dNdY[4][4];
	int conductivity = 30;
	int convection= 25;
	int specyficHeat = 700;
	int density = 7800;
	double matrixH[4][4];
	double matrixC[4][4];
	

	void showNode(int number);
	void showElement(int number);
	void showGrid();
	void showDetailedGrid();

	void calculateJacobianForEveryElement();
	void calculateMatrixH();
	void calculateMatrixC();
	void calculateBorderConditionMatrix();


	Grid();
};