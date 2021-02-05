/*
 * lal.h
 *
 *  Created on: 14 Oca 2021
 *      Author: X550V
 */

#ifndef LAL_H_
#define LAL_H_

class Vector;

class Matrix {

public:
	Matrix();
	Matrix(double*, int, int);
	virtual ~Matrix();
	void Print();
	int getRow();
	int getCol();
	double* getEntries();

	Matrix operator+(const Matrix &m);
	Matrix operator-(const Matrix &m);
	Matrix operator*(const Matrix &m);
	Vector operator*(Vector &m);
	bool operator==(const Matrix &m);

	virtual void Transpose();
	virtual void Norm();
	
protected:
	int rows;
	int cols;
	double *entries;
};

class Vector: public Matrix {
public:
	Vector(double*, int);
	~Vector();
	void Print();
	int getRow();
	double getNorm();
	double* getEntries();
	
	Vector operator+(const Vector &v);
	Vector operator-(const Vector &v);
	Vector operator*(const Vector &v);
	Vector operator*(Matrix &m);
	bool operator==(const Vector &v);
	Vector operator++();
	Vector operator++(int);

protected:
	double l2norm;
};

Vector **CreateVectorArray(int,int);
void FillVectorArray(Vector** vectorArray, int m, int p);
Matrix vectorArray2Matrix(Vector**,int,int);


#endif /* LAL_H_ */

