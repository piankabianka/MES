#include "stdafx.h"
#include "Node.h"
#include "Element.h"
#include "GlobalData.h"
#include "Grid.h"
#include <iostream>
using namespace std;



Grid::Grid() {
	GlobalData data;

	nodesNumber = data.totalNodesNumber;
	elementNumber = data.totalElementNumber;

	nodesTab = new Node[data.totalNodesNumber];
	elementsTab = new Element[data.totalElementNumber];

	int nodeNumber = 1;

	for (int i = 0; i < data.totalElementNumber; i++) {
		if (nodeNumber % (data.nodesNumberH + 1) == 0 && i != 0) {
			nodeNumber++;
		}

		this->elementsTab[i].elementID = i + 1;

		this->elementsTab[i].nodesIDtab[0] = nodeNumber;
		this->elementsTab[i].nodesIDtab[1] = nodeNumber + 1;
		this->elementsTab[i].nodesIDtab[2] = nodeNumber + data.nodesNumberH;
		this->elementsTab[i].nodesIDtab[3] = nodeNumber + data.nodesNumberH + 1;
		nodeNumber++;
	}

	double tmp1, tmp2;
	double app = 0.000001;
	double x = 0;				//beginning of coordinate system for nodes
	double y = 0;

	for (int i = 0; i < data.totalNodesNumber; i++) {
		if ((i + 1) % data.nodesNumberH == 0 && i != 0) {
			this->nodesTab[i].nodeID = i + 1;
			this->nodesTab[i].x = x;
			this->nodesTab[i].y = y;
			x += data.dx;
			y = 0;
		}
		else {
			this->nodesTab[i].x = x;
			this->nodesTab[i].y = y;
			y += data.dy;
		}
		tmp1 = nodesTab[i].x;
		tmp2 = nodesTab[i].y;

		if (tmp1 <= app) {
			nodesTab[i].BC = true;
		}
		if (tmp1 >= data.W - app) {
			nodesTab[i].BC = true;
		}
		if (tmp2 <= app) {
			nodesTab[i].BC = true;
		}
		if (tmp2 >= data.H - app) {
			nodesTab[i].BC = true;
		}
	}
}

void Grid::showNode(int number) {
	cout << "\tNode number: " << number << endl;
	cout << "\t\tx:" << nodesTab[number - 1].x << endl;
	cout << "\t\ty:" << nodesTab[number - 1].y << endl;
	cout << "\t\tBC:" << nodesTab[number - 1].BC << endl;
}

void Grid::showElement(int number) {
	cout << "........................" << endl;
	cout << "Element number: " << number << endl;
	showNode(elementsTab[number - 1].nodesIDtab[0]);
	showNode(elementsTab[number - 1].nodesIDtab[1]);
	showNode(elementsTab[number - 1].nodesIDtab[2]);
	showNode(elementsTab[number - 1].nodesIDtab[3]);

}

void Grid::showGrid() {
	cout << "----------------------------" << endl;
	cout << "GRID WITH ELEMENTS:\n" << endl;
	for (int i = 0; i < elementNumber; i++) {
		cout << "Element " << i + 1 << " nodes are:\t\t(" << elementsTab[i].nodesIDtab[0] << "," << elementsTab[i].nodesIDtab[1] << "," << elementsTab[i].nodesIDtab[2] << "," << elementsTab[i].nodesIDtab[3] << ")" << endl;
	}
}

void Grid::showDetailedGrid() {
	cout << "----------------------------" << endl;
	cout << "DETAILED GRID: \n" << endl;
	for (int i = 0; i < elementNumber; i++) {
		showElement(i + 1);
	}
}

void Grid::calculateJacobianForEveryElement() {

	double tabX[4] = {0, 0.025, 0.025, 0};
	double tabY[4] = { 0, 0, 0.025, 0.025};

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == 0) {
				jMatrix[i][j] = universalElement.dNdKsiDerivativeMatrix[i][0] * tabX[0] + universalElement.dNdKsiDerivativeMatrix[i][1] * tabX[1] + universalElement.dNdKsiDerivativeMatrix[i][2] * tabX[2] + universalElement.dNdKsiDerivativeMatrix[i][3] * tabX[3];
			}
			if (i == 1) {
				jMatrix[i][j] = universalElement.dNdKsiDerivativeMatrix[i][0] * tabY[0] + universalElement.dNdKsiDerivativeMatrix[i][1] * tabY[1] + universalElement.dNdKsiDerivativeMatrix[i][2] * tabY[2] + universalElement.dNdKsiDerivativeMatrix[i][3] * tabY[3];
			}
			if (i == 2) {
				jMatrix[i][j] = universalElement.dNdEtaDerivativeMatrix[i][0] * tabX[0] + universalElement.dNdEtaDerivativeMatrix[i][1] * tabX[1] + universalElement.dNdEtaDerivativeMatrix[i][2] * tabX[2] + universalElement.dNdEtaDerivativeMatrix[i][3] * tabX[3];
			}
			if (i == 3) {
				jMatrix[i][j] = universalElement.dNdEtaDerivativeMatrix[i][0] * tabY[0] + universalElement.dNdEtaDerivativeMatrix[i][1] * tabY[1] + universalElement.dNdEtaDerivativeMatrix[i][2] * tabY[2] + universalElement.dNdEtaDerivativeMatrix[i][3] * tabY[3];
			}

		}
	}

	jMatrix[0][0] = universalElement.dNdKsiDerivativeMatrix[0][0] * tabX[0] + universalElement.dNdKsiDerivativeMatrix[0][1] * tabX[1] + universalElement.dNdKsiDerivativeMatrix[0][2] * tabX[2] + universalElement.dNdKsiDerivativeMatrix[0][3] * tabX[3];

	cout << "\nJMATRIX : " << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << jMatrix[i][j] << "\t";
		}
		cout << endl;
	}

	
	double tmp = 0;

	for (int i = 0; i < 4; i++) {
		detJ[i] = jMatrix[0][i] * jMatrix[3][i] - jMatrix[1][i] * jMatrix[2][i];
	}
	cout << "\nDET J: " << endl;
	for (int i = 0; i < 4; i++) {
		cout << detJ[i] << "\t";
	}
	cout << endl;

	int k = 3;


	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {

			jMatrixReverse[i][j] = jMatrix[k][j] / detJ[i];
		}
		k--;
	}

	cout << "\nJMATRIX REVERSE: " << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << jMatrixReverse[i][j] << "\t";
		}
		cout << endl;
	}



	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			dNdX[i][j] = jMatrixReverse[0][i] * universalElement.dNdKsiDerivativeMatrix[i][j] + jMatrixReverse[1][i] * universalElement.dNdEtaDerivativeMatrix[i][j];
			dNdY[i][j] = jMatrixReverse[2][i] * universalElement.dNdKsiDerivativeMatrix[i][j] + jMatrixReverse[3][i] * universalElement.dNdEtaDerivativeMatrix[i][j];

		}
	}



	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << dNdX[i][j] << "\t\t";
		}
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << dNdY[i][j] << "\t\t";
		}
		cout << endl;
	}

}

void Grid::calculateMatrixH() {
	double dNdXP1[4][4];
	double dNdXP2[4][4];
	double dNdXP3[4][4];
	double dNdXP4[4][4];

	double dNdYP1[4][4];
	double dNdYP2[4][4];
	double dNdYP3[4][4];
	double dNdYP4[4][4];


	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			dNdXP1[i][j] = dNdX[0][i]*dNdX[0][j]*detJ[0] * conductivity;
			dNdXP2[i][j] = dNdX[1][i] * dNdX[1][j] * detJ[1] * conductivity;
			dNdXP3[i][j] = dNdX[2][i] * dNdX[2][j]*detJ[2] * conductivity;
			dNdXP4[i][j] = dNdX[3][i] * dNdX[3][j]*detJ[3] * conductivity;

			dNdYP1[i][j]= dNdY[0][i] * dNdY[0][j] * detJ[0] * conductivity;
			dNdYP2[i][j] = dNdY[1][i] * dNdY[1][j] * detJ[1] * conductivity;
			dNdYP3[i][j] = dNdY[2][i] * dNdY[2][j] * detJ[2] * conductivity;
			dNdYP4[i][j] = dNdY[3][i] * dNdY[3][j]* detJ[3] * conductivity;

			matrixH[i][j] = dNdXP1[i][j] + dNdYP1[i][j] + dNdXP2[i][j] + dNdYP2[i][j] + dNdXP3[i][j] + dNdYP3[i][j] + dNdXP4[i][j] + dNdYP4[i][j];
		}
	}
	
	
	cout << "\nH MATRIX" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << matrixH[i][j] << "\t";
		}
		cout << endl;
	}
}

void Grid::calculateMatrixC() {
	double p1Matrix[4][4];
	double p2Matrix[4][4];
	double p3Matrix[4][4];
	double p4Matrix[4][4];
	UniversalElement universalElement;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			p1Matrix[i][j]=universalElement.shapeFunctionsMatrix[0][i]*universalElement.shapeFunctionsMatrix[0][j] *detJ[0]*specyficHeat*density;
			p2Matrix[i][j]= universalElement.shapeFunctionsMatrix[1][i] * universalElement.shapeFunctionsMatrix[1][j] * detJ[0] * specyficHeat*density;
			p3Matrix[i][j]=universalElement.shapeFunctionsMatrix[2][i] * universalElement.shapeFunctionsMatrix[2][j] * detJ[0] * specyficHeat*density;
			p4Matrix[i][j]= universalElement.shapeFunctionsMatrix[3][i] * universalElement.shapeFunctionsMatrix[3][j] * detJ[0] * specyficHeat*density;
			matrixC[i][j] = p1Matrix[i][j] + p2Matrix[i][j] + p3Matrix[i][j] + p4Matrix[i][j];
		}
	}

	cout << "\nC MATRIX" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << matrixC[i][j] << "\t";
		}
		cout << endl;
	}
}

void Grid::calculateBorderConditionMatrix() {
	
	UniversalElement uE;
	double localBorderMatrix[4][4]{ 0 };
	
	int number = 1;
	double tabX[4];
	double tabY[4];
	int bc[4];
	double J = 0;
	GlobalData gD;
	double dx=gD.dx;
	double dy=gD.dy;
	for (int i = 0; i < 4; i++) {
		tabX[i] = nodesTab[elementsTab[number - 1].nodesIDtab[i] - 1].x;
		tabY[i] = nodesTab[elementsTab[number - 1].nodesIDtab[i] - 1].y;
		bc[i] = nodesTab[elementsTab[number - 1].nodesIDtab[i] - 1].BC;
		
	}
	
	if (bc[0] == 1 && bc[1] == 1) {
		J = dx / 2;
		localBorderMatrix[0][0] += uE.calculateshapeFunctionMatrixValue(1, -1/sqrt(3), -1)*uE.calculateshapeFunctionMatrixValue(1, -1 / sqrt(3), -1)*J*convection;
		localBorderMatrix[0][1] += uE.calculateshapeFunctionMatrixValue(1, -1 / sqrt(3), -1)*uE.calculateshapeFunctionMatrixValue(2, -1 / sqrt(3), -1)*J*convection;
		localBorderMatrix[1][0] += uE.calculateshapeFunctionMatrixValue(2, -1 / sqrt(3), -1)*uE.calculateshapeFunctionMatrixValue(1, -1 / sqrt(3), -1)*J*convection;
		localBorderMatrix[1][1] += uE.calculateshapeFunctionMatrixValue(2,- 1 / sqrt(3), -1)*uE.calculateshapeFunctionMatrixValue(2, -1 / sqrt(3), -1)*J*convection;

		localBorderMatrix[0][0] += uE.calculateshapeFunctionMatrixValue(1, 1 / sqrt(3), -1)*uE.calculateshapeFunctionMatrixValue(1, 1 / sqrt(3), -1)*J*convection;
		localBorderMatrix[0][1] += uE.calculateshapeFunctionMatrixValue(1, 1 / sqrt(3), -1)*uE.calculateshapeFunctionMatrixValue(2,1 / sqrt(3), -1)*J*convection;
		localBorderMatrix[1][0] += uE.calculateshapeFunctionMatrixValue(2, 1 / sqrt(3), -1)*uE.calculateshapeFunctionMatrixValue(1,1 / sqrt(3), -1)*J*convection;
		localBorderMatrix[1][1] += uE.calculateshapeFunctionMatrixValue(2, 1 / sqrt(3), -1)*uE.calculateshapeFunctionMatrixValue(2,1 / sqrt(3), -1)*J*convection;
	}

	if (bc[1] == 1 && bc[2] == 1) {
		J = dy / 2;
		localBorderMatrix[1][1] += uE.calculateshapeFunctionMatrixValue(2, 1, -1 / sqrt(3))*uE.calculateshapeFunctionMatrixValue(2, 1, -1 / sqrt(3))*J*convection;
		localBorderMatrix[1][2] += uE.calculateshapeFunctionMatrixValue(2, 1, -1 / sqrt(3))*uE.calculateshapeFunctionMatrixValue(3, 1, -1 / sqrt(3))*J*convection;
		localBorderMatrix[2][1] += uE.calculateshapeFunctionMatrixValue(3, 1, -1 / sqrt(3))*uE.calculateshapeFunctionMatrixValue(2, 1, -1 / sqrt(3))*J*convection;
		localBorderMatrix[2][2] += uE.calculateshapeFunctionMatrixValue(3, 1, -1 / sqrt(3))*uE.calculateshapeFunctionMatrixValue(3, 1, -1 / sqrt(3))*J*convection;

		localBorderMatrix[1][1] += uE.calculateshapeFunctionMatrixValue(2, 1, 1 / sqrt(3))*uE.calculateshapeFunctionMatrixValue(2, 1, 1 / sqrt(3))*J*convection;
		localBorderMatrix[1][2] += uE.calculateshapeFunctionMatrixValue(2, 1, 1 / sqrt(3))*uE.calculateshapeFunctionMatrixValue(3, 1, 1 / sqrt(3))*J*convection;
		localBorderMatrix[2][1] += uE.calculateshapeFunctionMatrixValue(3, 1, 1 / sqrt(3))*uE.calculateshapeFunctionMatrixValue(2, 1, 1 / sqrt(3))*J*convection;
		localBorderMatrix[2][2] += uE.calculateshapeFunctionMatrixValue(3, 1, 1 / sqrt(3))*uE.calculateshapeFunctionMatrixValue(3, 1, 1 / sqrt(3))*J*convection;
	}

	if (bc[2] == 1 && bc[3] == 1) {
		J = dx / 2;
		localBorderMatrix[2][2] += uE.calculateshapeFunctionMatrixValue(3, 1 / sqrt(3), 1)*uE.calculateshapeFunctionMatrixValue(3, 1 / sqrt(3), 1)*J*convection;
		localBorderMatrix[2][3] += uE.calculateshapeFunctionMatrixValue(3, 1 / sqrt(3), 1)*uE.calculateshapeFunctionMatrixValue(4, 1 / sqrt(3), 1)*J*convection;
		localBorderMatrix[3][2] += uE.calculateshapeFunctionMatrixValue(4, 1 / sqrt(3), 1)*uE.calculateshapeFunctionMatrixValue(3, 1 / sqrt(3), 1)*J*convection;
		localBorderMatrix[3][3] += uE.calculateshapeFunctionMatrixValue(4,1 / sqrt(3), 1)*uE.calculateshapeFunctionMatrixValue(4, 1 / sqrt(3), 1)*J*convection;

		localBorderMatrix[2][2] += uE.calculateshapeFunctionMatrixValue(3, -1 / sqrt(3), 1)*uE.calculateshapeFunctionMatrixValue(3,- 1 / sqrt(3), 1)*J*convection;
		localBorderMatrix[2][3] += uE.calculateshapeFunctionMatrixValue(3, -1 / sqrt(3), 1)*uE.calculateshapeFunctionMatrixValue(4, -1 / sqrt(3), 1)*J*convection;
		localBorderMatrix[3][2] += uE.calculateshapeFunctionMatrixValue(4, -1 / sqrt(3), 1)*uE.calculateshapeFunctionMatrixValue(3, -1 / sqrt(3), 1)*J*convection;
		localBorderMatrix[3][3] += uE.calculateshapeFunctionMatrixValue(4, -1 / sqrt(3), 1)*uE.calculateshapeFunctionMatrixValue(4, -1 / sqrt(3), 1)*J*convection;
	}
	if (bc[3] == 1 && bc[0] == 1) {
		J = dy / 2;
		localBorderMatrix[0][0] += uE.calculateshapeFunctionMatrixValue(1, -1, 1 / sqrt(3))*uE.calculateshapeFunctionMatrixValue(1, -1, 1 / sqrt(3))*J*convection;
		localBorderMatrix[0][3] += uE.calculateshapeFunctionMatrixValue(1, -1, 1 / sqrt(3))* uE.calculateshapeFunctionMatrixValue(4, -1, 1 / sqrt(3))*J*convection;
		localBorderMatrix[3][0] += uE.calculateshapeFunctionMatrixValue(4, -1, 1 / sqrt(3))*uE.calculateshapeFunctionMatrixValue(1, -1,1 / sqrt(3))*J*convection;
		localBorderMatrix[3][3] += uE.calculateshapeFunctionMatrixValue(4, -1, 1 / sqrt(3))*uE.calculateshapeFunctionMatrixValue(4, -1,1 / sqrt(3))*J*convection;

		localBorderMatrix[0][0] += uE.calculateshapeFunctionMatrixValue(1, -1, -1 / sqrt(3))*uE.calculateshapeFunctionMatrixValue(1, -1, -1 / sqrt(3))*J*convection;
		localBorderMatrix[0][3] += uE.calculateshapeFunctionMatrixValue(1, -1, -1 / sqrt(3))* uE.calculateshapeFunctionMatrixValue(4, -1, -1 / sqrt(3))*J*convection;
		localBorderMatrix[3][0] += uE.calculateshapeFunctionMatrixValue(4, -1, -1 / sqrt(3))*uE.calculateshapeFunctionMatrixValue(1, -1, -1 / sqrt(3))*J*convection;
		localBorderMatrix[3][3] += uE.calculateshapeFunctionMatrixValue(4, -1, -1 / sqrt(3))*uE.calculateshapeFunctionMatrixValue(4, -1, -1 / sqrt(3))*J*convection;
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << localBorderMatrix[i][j] << "\t";
		}
		cout << endl;
	}
	
}