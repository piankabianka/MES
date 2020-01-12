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
	
	
	void showNode(int number);
	void showElement(int number);
	void showGrid();
	void showDetailedGrid();

	void calculateJacobianForEveryElement();
	

	Grid();
};

Grid::Grid() {
	GlobalData data;

	nodesNumber = data.totalNodesNumber;
	elementNumber = data.totalElementNumber;

	nodesTab = new Node[data.totalNodesNumber];
	elementsTab = new Element[data.totalElementNumber];

	int nodeNumber = 1;

	for (int i = 0; i < data.totalElementNumber; i++) {
		if (nodeNumber % (data.nodesNumberH + 1)==0 && i != 0) {
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
		if (tmp1 >=data.W-app) {
			nodesTab[i].BC = true;
		}
		if (tmp2 <= app) {
			nodesTab[i].BC = true;
		}
		if (tmp2 >= data.H-app) {
			nodesTab[i].BC = true;
		}
	}
}

void Grid::showNode(int number) {
		cout << "\tNode number: " << number << endl;
		cout << "\t\tx:" << nodesTab[number - 1].x << endl;
		cout << "\t\ty:" << nodesTab[number - 1].y << endl;
		cout<<"\t\tBC:"<< nodesTab[number - 1].BC << endl;
}

void Grid::showElement(int number) {
		cout << "........................" << endl;
		cout << "Element number: " << number << endl;
		showNode(elementsTab[number-1].nodesIDtab[0]);
		showNode(elementsTab[number-1].nodesIDtab[1]);
		showNode(elementsTab[number-1].nodesIDtab[2]);
		showNode(elementsTab[number-1].nodesIDtab[3]);
	
}

void Grid::showGrid() {
	cout << "----------------------------" << endl;
	cout << "GRID WITH ELEMENTS:\n" << endl;
	for (int i = 0; i < elementNumber; i++) {
		cout << "Element " << i + 1 << " nodes are:\t\t(" << elementsTab[i].nodesIDtab[0] << "," << elementsTab[i].nodesIDtab[1] << "," << elementsTab[i].nodesIDtab[2] << "," << elementsTab[i].nodesIDtab[3] <<")"<< endl;
	}
}

void Grid::showDetailedGrid() {
	cout << "----------------------------" << endl;
	cout << "DETAILED GRID: \n" << endl;
	for (int i = 0; i < elementNumber; i++) {
		showElement(i + 1);
	}
}

void Grid::calculateJacobianForEveryElement(){

	double tabX[4];
	double tabY[4];

	for (int i = 0; i < 4; i++) {
		if (i == 0) {
			tabX[i] = 0;
			tabY[i] = 0;
		}
		if (i == 1) {
			tabX[i] = 0.025;
			tabY[i] = 0;
		}
		if (i == 2) {
			tabX[i] = 0.025;
			tabY[i] = 0.025;
		}
		if (i == 3) {
			tabX[i] = 0;
			tabY[i] = 0.025;
		}
	}

	double jMatrix[4][4];

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == 0) {
				jMatrix[i][j] = universalElement.dNdKsiDerivativeMatrix[i][0]*tabX[0] + universalElement.dNdKsiDerivativeMatrix[i][1] *tabX[1] + universalElement.dNdKsiDerivativeMatrix[i][2] *tabX[2] + universalElement.dNdKsiDerivativeMatrix[i][3] *tabX[3];
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

	double detJ[4];
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

	double jMatrixReverse[4][4];
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

	double dNdX[4][4];
	double dNdY[4][4];

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