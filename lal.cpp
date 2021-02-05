/*
 * lal.cpp
 *
 *  Created on: 14 Oca 2021
 *      Author: X550V
 */

#include <iostream>
#include "lal.h"
#include <math.h>
#include <cstdlib>
#include <typeinfo>
#include <time.h>

using namespace std;

Matrix::Matrix(double *entries, int rows, int cols) {
	this->entries = entries;
	this->rows = rows;
	this->cols = cols;
}

Matrix::Matrix(){
	double* entries = (double*) calloc(1,sizeof(double));
	this->entries = entries;
}


Matrix::~Matrix() {
	delete this->entries;
}

int Matrix::getRow() {
	return this->rows;
}
int Matrix::getCol() {
	return this->cols;
}
double* Matrix::getEntries() {
	return this->entries;
}

void Matrix::Print() {
	int counter = 0;

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; ++j) {
			cout << *(entries + counter++) << " ";
		}
		cout << endl;
	}
}

Matrix Matrix::operator+(const Matrix &m) {
	if (this->rows == m.rows && this->cols == m.cols) {
		int limit = this->cols * this->rows;
		double *entries = (double*) calloc(limit, sizeof(double));
		for (int i = 0; i < limit; i++) {
			entries[i] = (double) this->entries[i] + (double) m.entries[i];
		}
		Matrix sum(entries, m.rows, m.cols);
		sum.Print();
		return sum;
	} else {
		cout << "Operation is invalid. Incompatible lengths." << endl;
		return m;
	}
}

Matrix Matrix::operator-(const Matrix &m) {
	if (this->rows == m.rows && this->cols == m.cols) {
		int limit = this->cols * this->rows;
		double *entries = (double*) calloc(limit, sizeof(double));
		for (int i = 0; i < limit; i++) {
			entries[i] = (double) this->entries[i] - (double) m.entries[i];
		}
		Matrix sum(entries, m.rows, m.cols);
		return sum;
	} else {
		cout << "Operation is invalid. Incompatible lengths." << endl;
		return m;
	}
}

Matrix Matrix::operator*(const Matrix &m) {
	if (this->cols == m.rows) {

		double *entries = (double*) calloc(this->rows * m.cols, sizeof(double));
		int counter = 0;
		for (int i = 0; i < this->rows; i++) {
			for (int j = 0; j < m.cols; j++) {
				for (int k = 0; k < this->cols; k++) {
					entries[counter] += this->entries[i * this->cols + k]
							* m.entries[k * m.cols + j];
				}
				counter++;
			}

		}
		Matrix sum(entries, this->rows, m.cols);
		return sum;
	} else {
		cout << "Operation is invalid. Incompatible lengths." << endl;
		return m;
	}
}

Vector Matrix::operator *(Vector &v) {
	if (v.getRow() == this->cols) {
		double *entries = (double*) calloc(this->rows, sizeof(double));
		int counter = 0;
		double sum = 0.0;
		for (int i = 0; i < this->rows; i++) {
			entries[counter] = 0;
			for (int j = 0; j < this->cols; j++) {
				sum += (double) v.getEntries()[j]
						* (double) this->entries[j + i * this->cols];
			}
			*(entries + counter) = (double) sum;
			sum = 0;
			counter++;
		}
		Vector sumVector(entries, this->rows);
		return sumVector;
	} else {
		cout << "Operation is invalid. Incompatible lengths." << endl;
		return Vector(this->entries, this->rows);
	}
}

bool Matrix::operator ==(const Matrix &m) {
	bool isSame = true;
	if (this->rows == m.rows && this->cols == m.cols) {
		for (int i = 0; i < m.rows * m.cols; ++i) {
			if (this->entries[i] != m.entries[i]) {
				isSame = false;
			}
		}
	} else {
		isSame = false;
	}
	return isSame;
}


//-----------------------------------------------------------------------------//

Vector::Vector(double *entries, int rows) :
		Matrix(entries, rows, 1) {
	double l2Norm = 0;
	for (int i = 0; i < rows; ++i) {
		l2Norm += pow(entries[i], 2);
	}
	this->l2norm = sqrt(l2Norm);
}

Vector::~Vector() {
	delete this->entries;
}

int Vector::getRow() {
	return this->rows;
}

double* Vector::getEntries() {
	return this->entries;
}

double Vector::getNorm(){
	return this->l2norm;
}

void Vector::Print() {
	for (int i = 0; i < rows; i++) {
		cout << this->entries[i] << " ";
	}
	cout << endl;
}

Vector Vector::operator+(const Vector &v) {
	if (v.rows == this->rows) {
		double *entries = (double*) calloc(v.rows, sizeof(double));
		for (int i = 0; i < this->rows; ++i) {
			entries[i] = v.entries[i] + this->entries[i];
		}
		Vector sum(entries, v.rows);
		return sum;
	} else {
		cout << "Operation is invalid. Incompatible lengths." << endl;
		return v;
	}
}

Vector Vector::operator-(const Vector &v) {
	if (v.rows == this->rows) {
		double *entries = (double*) calloc(v.rows, sizeof(double));
		for (int i = 0; i < this->rows; ++i) {
			entries[i] = v.entries[i] - this->entries[i];
		}
		Vector sum(entries, v.rows);
		return sum;
	} else {
		cout << "Operation is invalid. Incompatible lengths." << endl;
		return v;
	}
}

Vector Vector::operator*(const Vector &v) {
	if (this->rows == v.rows) {
		double *entries = (double*) calloc(v.rows, sizeof(double));

		for (int i = 0; i < this->rows; i++) {
			entries[i] = this->entries[i] * v.entries[i];
		}
		Vector sum(entries, v.rows);
		return sum;

	} else {
		cout << "Operation is invalid. Incompatible lengths." << endl;
		return v;
	}
}

Vector Vector::operator *(Matrix &m) {
	if (this->rows == m.getCol()) {
		double *entries = (double*) calloc(m.getRow(), sizeof(double));
		int counter = 0;
		double sum = 0.0;
		for (int i = 0; i < m.getRow(); i++) {
			entries[counter] = 0;
			for (int j = 0; j < m.getCol(); j++) {
				sum += (double) this->entries[j]
						* (double) m.getEntries()[j + i * m.getCol()];
			}
			*(entries + counter) = (double) sum;
			sum = 0;
			counter++;
		}
		Vector sumVector(entries, m.getRow());
		return sumVector;
	} else {
		cout << "Operation is invalid. Incompatible lengths." << endl;
		return Vector(this->entries, this->rows);
	}
}

bool Vector::operator ==(const Vector &v) {
	bool isSame = true;
	if (this->rows == v.rows) {
		for (int i = 0; i < v.rows; ++i) {
			if (this->entries[i] != v.entries[i]) {
				isSame = false;
			}
		}
	} else {
		isSame = false;
	}
	return isSame;
}

//Pre-Increment
Vector Vector::operator ++() {
	for (int i = 0; i < this->rows; ++i) {
		this->entries[i]++;
	}
	cout << "pre";
	return *this;
}

//Post-Increment
Vector Vector::operator ++(int) {
	const Vector old(this->entries, this->rows);
	for (int i = 0; i < this->rows; ++i) {
		++*(old.entries + i);
	}
	cout << "post";
	return old;
}

void Matrix::Transpose() {
	if (typeid(*this).name() == typeid(Vector).name()) {
		double *entries = (double*) calloc(this->rows, sizeof(double));
		for (int i = 0; i < this->rows; i++) {
			entries[i] = this->entries[this->rows - i - 1];
		}
		this->entries = entries;
	} else {
		double *entries = (double*) calloc(this->rows * this->cols,
				sizeof(double));
		for (int i = 0; i < this->rows * this->cols; i++) {
			entries[i] = this->entries[this->rows * this->cols - i - 1];
		}
		this->entries = entries;
	}
	this->Print();
}

void Matrix::Norm() {
	if (typeid(*this).name() == typeid(Vector).name()) {
		cout << dynamic_cast<Vector *>(this)->getNorm();
	} else {
		double l2Norm = 0.0;
		for (int i = 0; i < this->rows * this->cols; ++i) {
			l2Norm += pow(entries[i], 2);
		}
		l2Norm = sqrt(l2Norm);
		cout << l2Norm ;
	}
}

Vector **CreateVectorArray(int m,int p){
	Vector **vectorArray = (Vector**) calloc(m, sizeof(Vector*));  
	for (int i = 0; i < m; i++) {           
		vectorArray[i] = (Vector*) calloc(p, sizeof(Vector));
	}
	return vectorArray;
}

void FillVectorArray(Vector** vectorArray, int m, int p){
	for(int i = 0; i < m; i++){
		double *entries = (double*) calloc(p, sizeof(double));
		for(int j = 0; j < p; j++){
			entries[j] = (double) rand() / (double) (RAND_MAX) + 1.0; 
		}
		Vector *v = new Vector(entries,p);
		vectorArray[i] = v;
		delete[]entries;
	}
}

Matrix vectorArray2Matrix(Vector** vectorArray, int m, int p){
	double* entries = (double*) calloc(m * p, sizeof(double));
	for(int i = 0; i < m; i++){
		for(int j = 0; j < p; j++){
			entries[i + j * m] = vectorArray[i]->getEntries()[j];
		}
	}
	Matrix newMatrix(entries, p, m);
	return newMatrix;
}
