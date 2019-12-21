#pragma once
#include "stdafx.h"
#include <iostream>
#include <fstream>
using namespace std;

struct GlobalData {
	double H;			//height
	double W;			//width
	int nodesNumberH;	//number of nodes upright
	int nodesNumberW;	//number of nodes horizontally
	int totalNodesNumber;
	int totalElementNumber;
	double dx;
	double dy;
	GlobalData();
	void showData();
};

GlobalData::GlobalData() {
	fstream file;
	file.open("data.txt");

	if (file.good()) {
		while (!file.eof()) {
			file >> H;
			file >> W;
			file >> nodesNumberH;
			file >> nodesNumberW;
		}

		totalNodesNumber = nodesNumberH * nodesNumberW;
		totalElementNumber = (nodesNumberH - 1)*(nodesNumberW - 1);

		dx = (W / (nodesNumberW - 1));
		dy = (H / (nodesNumberH - 1));

	}

	else
		cout << "Plik sie wyjebal" << endl;
	file.close();
}

void GlobalData::showData() {
	cout << "----------------------------" << endl;
	cout << "GLOBAL DATA: \n" << endl;
	cout << "Height:\t\t\t\t" << H << endl;
	cout << "Width:\t\t\t\t" << W << endl;
	cout << "Number of nodes upright:\t" << nodesNumberH << endl;
	cout << "Number of nodes horizontally:\t" << nodesNumberW<< endl;
	cout << "Total number of nodes:\t\t" << totalNodesNumber << endl;
	cout << "Total number of elements:\t" << totalElementNumber << endl;
	cout << "dx:\t\t\t\t" << dx << endl;
	cout << "dy:\t\t\t\t" << dy << endl;
}