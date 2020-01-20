#pragma once
#include "stdafx.h"
#include <iostream>
using namespace std;

struct UniversalElement {
	double ksi;
	double eta;

	double weight[2];
	double integralPoints[2];

	double shapeFunctionsMatrix[4][4];
	double dNdKsiDerivativeMatrix[4][4];
	double dNdEtaDerivativeMatrix[4][4];

	void showShapeFunctionsMatrix();
	void showdNdKsiDerivativeMatrix();
	void showdNdEtaDerivativeMatrix();
	
	UniversalElement();
};
