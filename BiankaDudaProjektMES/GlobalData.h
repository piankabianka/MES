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
