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
	int nodesNumberH;	//number of nodes upright
	int nodesNumberW;	//number of nodes horizontally
	UniversalElement universalElement;
	double jMatrix[4][4];
	double detJ[4];
	double jMatrixReverse[4][4];
	double dNdX[4][4];
	double dNdY[4][4];
	int conductivity;
	int convection;
	int specyficHeat;
	int density;
	int ambientTemp;
	int timeStep;
	int simulationTime;
	double initialTemp;
	double matrixH[4][4];
	double matrixC[4][4];
	double matrixP[4];
	double localBorderMatrix[4][4]{ 0 };
	double **matrixHGlobal;
	double **matrixCGlobal;
	double *matrixPGlobal;
	double *t0Vector;
	double *t1Vector;

	void showNode(int number);
	void showElement(int number);
	void showGrid();
	void showDetailedGrid();

	void calculateJacobianForEveryElement(int numberOfElement);
	void calculateMatrixH();
	void calculateMatrixC();
	void calculateMatrixP(int numberOfElement);
	void calculateBorderConditionMatrix(int numberOfElement);
	void addMatrixHandBC();
	void initGlobalMatrixes();
	void calculateGlobalMatrixes();
	void calculations();
	void calculateSolution();


	Grid();
};