#pragma once
#include "stdafx.h"
#include <iostream>
using namespace std;

struct Node {
	double x;
	double y;
	int nodeID;
	bool BC;
	int t;
	Node();
};

Node::Node() {
	x = 0;
	y = 0;
	nodeID = 0;
	BC = false;
	t = 0;
}