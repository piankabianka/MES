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
	
	Grid();
	void showNode(int number);
	void showElement(int number);
	void showGrid();
	void showDetailedGrid();
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