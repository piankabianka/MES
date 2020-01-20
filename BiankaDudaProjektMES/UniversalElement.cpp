#include "stdafx.h"
#include "UniversalElement.h"
#include <iostream>
using namespace std;


UniversalElement::UniversalElement() {
	weight[0] = -1;
	weight[1] = 1;

	integralPoints[0] = -1 / sqrt(3);
	integralPoints[1] = 1 / sqrt(3);

	for (int i = 0; i < 4; i++) {
		if (i == 0) {
			ksi = integralPoints[0];
			eta = integralPoints[0];
		}
		else if (i == 1) {
			ksi = integralPoints[1];
			eta = integralPoints[0];
		}
		else if (i == 2) {
			ksi = integralPoints[1];
			eta = integralPoints[1];
		}
		else if (i == 3) {
			ksi = integralPoints[0];
			eta = integralPoints[1];
		}

		shapeFunctionsMatrix[i][0] = (0.25 * (1 - ksi) * (1 - eta));
		shapeFunctionsMatrix[i][1] = (0.25 * (1 + ksi) * (1 - eta));
		shapeFunctionsMatrix[i][2] = (0.25 * (1 + ksi) * (1 + eta));
		shapeFunctionsMatrix[i][3] = (0.25 * (1 - ksi) * (1 + eta));

		dNdKsiDerivativeMatrix[i][0] = (-0.25 * (1 - eta));
		dNdKsiDerivativeMatrix[i][1] = (0.25 * (1 - eta));
		dNdKsiDerivativeMatrix[i][2] = (0.25 * (1 + eta));
		dNdKsiDerivativeMatrix[i][3] = (-0.25 * (1 + eta));

		dNdEtaDerivativeMatrix[i][0] = (-0.25 * (1 - ksi));
		dNdEtaDerivativeMatrix[i][1] = (-0.25 * (1 + ksi));
		dNdEtaDerivativeMatrix[i][2] = (0.25 * (1 + ksi));
		dNdEtaDerivativeMatrix[i][3] = (0.25 * (1 - ksi));

	}

}

void UniversalElement::showShapeFunctionsMatrix() {
	cout << "----------------------------" << endl;
	cout << "SHAPE FUNCTION MATRIX:\n" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << shapeFunctionsMatrix[i][j] << "\t";
		}
		cout << endl;
	}
}

void UniversalElement::showdNdKsiDerivativeMatrix() {
	cout << "----------------------------" << endl;
	cout << "KSI DERIVATIVE MATRIX:\n" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << dNdKsiDerivativeMatrix[i][j] << "\t";
		}
		cout << endl;
	}
}

void UniversalElement::showdNdEtaDerivativeMatrix() {
	cout << "----------------------------" << endl;
	cout << "ETA DERIVATIVE MATRIX:\n" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << dNdEtaDerivativeMatrix[i][j] << "\t";
		}
		cout << endl;
	}
}