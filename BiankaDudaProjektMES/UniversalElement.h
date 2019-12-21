#pragma once
#include "stdafx.h"
#include <iostream>
using namespace std;

struct UniversalElement {
	double ksi;
	double eta;
	double interval[2];
	double integralPoints[2];			// 2-point integral scheme
	double shapeFunctionsMatrix[4][4];
	double ksiDerivativeMatrix[4][4];		// dN/dKsi
	double etaDerivativeMatrix[4][4];		//dN/dEta

	double returnKsi(int number);
	double returnEta(int number);

	double calculateShapeFunctionsMatrixValue(double ksi, double eta, int number);
	double calculateKsiDerivativeMatrixValue(double eta, int number);
	double calculateEtaDerivativeMatrixValue(double ksi, int number);

	void showShapeFunctionsMatrix();
	void showKsiDerivativeMatrix();
	void showEtaDerivativeMatrix();

	UniversalElement();
};



UniversalElement::UniversalElement() {
	interval[0] = -1;
	interval[1] = 1;

	integralPoints[0] = -1 / sqrt(3);
	integralPoints[1] = 1 / sqrt(3);

	ksi = 0;
	eta = 0;

	for (int i = 0; i < 4; i++) {
		
		ksi = returnKsi(i);
		eta = returnEta(i);

		for (int j = 0; j < 4; j++) {
			shapeFunctionsMatrix[i][j] = calculateShapeFunctionsMatrixValue(ksi, eta, j + 1);
			ksiDerivativeMatrix[i][j] = calculateKsiDerivativeMatrixValue(eta, j + 1);
			etaDerivativeMatrix[i][j] = calculateEtaDerivativeMatrixValue(ksi, j + 1);
		}
	}

	
}

double UniversalElement::returnKsi(int number) {
	if (number == 1 || number==4) {
		return integralPoints[0];
	}
	else {
		return integralPoints[1];
	}
	
}

double UniversalElement::returnEta(int number) {
	if (number == 1 || number == 2) {
		return integralPoints[0];
	}
	else {
		return integralPoints[1];
	}

}

double UniversalElement::calculateShapeFunctionsMatrixValue(double ksi, double eta, int number) {
	if (number == 1) {
		return (0.25*(1 - ksi)*(1 - eta));
	}
	else if (number == 2) {
		return (0.25*(1 + ksi)*(1 - eta));
	}
	else if (number == 3) {
		return (0.25*(1 + ksi)*(1 + eta));
	}
	else {
		return (0.25*(1 - ksi)*(1 + eta));
	}
}

double UniversalElement::calculateKsiDerivativeMatrixValue(double eta, int number) {
	if (number == 1) {
		return (-0.25*(1 - eta));
	}
	else if (number == 2) {
		return (0.25*(1 - eta));
	}
	else if (number == 3) {
		return (0.25*(1 + eta));
	}
	else {
		return (-0.25*(1 + eta));
	}
}

double UniversalElement::calculateEtaDerivativeMatrixValue(double ksi, int number) {
	if (number == 1) {
		return (-0.25*(1 - ksi));
	}
	else if (number == 2) {
		return (-0.25*(1 + ksi));
	}
	else if (number == 3) {
		return (0.25*(1 + ksi));
	}
	else {
		return (0.25*(1 - ksi));
	}
}

void UniversalElement::showShapeFunctionsMatrix(){
	cout << "----------------------------" << endl;
	cout << "SHAPE FUNCTION MATRIX:\n" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << shapeFunctionsMatrix[i][j] << "\t";
		}
		cout << endl;
	}
}

void UniversalElement::showKsiDerivativeMatrix() {
	cout << "----------------------------" << endl;
	cout << "KSI DERIVATIVE MATRIX:\n" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << ksiDerivativeMatrix[i][j] << "\t";
		}
		cout << endl;
	}
}

void UniversalElement::showEtaDerivativeMatrix() {
	cout << "----------------------------" << endl;
	cout << "ETA DERIVATIVE MATRIX:\n" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << etaDerivativeMatrix[i][j] << "\t";
		}
		cout << endl;
	}
}