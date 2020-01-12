#include "stdafx.h"
#include <iostream>
using namespace std;

struct Matrix {
	int size;
	double **matrix;
	double determinant;

	void printMatrix();
	double calculateDeterminant2();
	double getMatrixValue(int a, int b);
	void setMatrixValue(int a, int b, double value);
	double getDeterminant();
	void multiplayMatrixByScalar(double scalar);

	Matrix();
	Matrix(int size);
};

Matrix::Matrix() {
	size = 0;
	determinant = 0;
}

Matrix::Matrix(int size) {
	this->size = size;
	

	matrix = new double*[size];
	
	for (int i = 0; i < size; i++) {
		matrix[i] = new double[size];
	}

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++){
			matrix[i][j] = j;
		}
	}

	if (size == 2) {
		determinant = calculateDeterminant2();
	}
	else
		determinant = 0;
}

void Matrix::printMatrix() {
	
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			cout << matrix[i][j] << "\t";
		}
		cout << endl;
	}
}

double Matrix::calculateDeterminant2() {
	determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
	return determinant;
}

double Matrix::getMatrixValue(int a, int b) {
	return matrix[a][b];
}

void Matrix::setMatrixValue(int a, int b, double value) {
	matrix[a][b] = value;
	calculateDeterminant2();
}

double Matrix::getDeterminant() {
	return determinant;
}

void Matrix::multiplayMatrixByScalar(double scalar) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			matrix[i][j] = matrix[i][j] * scalar;
		}
	}
	calculateDeterminant2();
}