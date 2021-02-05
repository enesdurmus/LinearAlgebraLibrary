/*
 * test.cpp
 *
 *  Created on: 14 Oca 2021
 *      Author: X550V
 */

#include <iostream>
#include "lal.h"
#include <typeinfo>
#include <cstdlib>

using namespace std;

int main(int argc, char **argv) {

	double *entries = (double*) calloc(4, sizeof(double));
	for (int i = 0; i < 4; i++) {
		entries[i] = i + 1;
	}
	double *entries2 = (double*) calloc(8, sizeof(double));
	for (int i = 0; i < 8; i++) {
		entries2[i] = (double) i + 5;
	}
	Matrix m1(entries, 2, 2);
	m1.Print();

	Matrix m2(entries2, 2, 4);
	m2.Print();
	// bool a = m == m2;
	// cout << a;
	m2 = m1 + m2;
	m2.Print();
	//m2= m-m2;
//	 m2 = m * m2;
//	 m2.Print();

	Vector v(entries, 4);
//	v->Print();
	//++v;
	m2.Transpose();
	m1.Norm();
	//m = &v;
	//cout << typeid(*m).name();
	//Vector v2(entries, 4);
	// v2.Print();
	//a = v == v2;
	//cout << a;
	///v = v + v;
	// v = v- v;
//	v = v * m;
//	v.Print();
	cout << endl;
	
	int m = 2, p = 4;
	Vector** vectorArray = CreateVectorArray(m,p);
	FillVectorArray(vectorArray, m, p);
	for(int i = 0; i < m;i++){
		vectorArray[i]->Print();
	} 
	Matrix m3 = vectorArray2Matrix(vectorArray, m, p);
	m3.Print();
	
	Matrix dizi[5];
	dizi[0] = Matrix(entries, 2, 2);
	dizi[1] = Matrix(entries, 2, 2);
	dizi[2] = Vector(entries, 4);
	dizi[3] = Vector(entries, 4);
	dizi[4] = Vector(entries, 4);

	cout << typeid(dizi[1]).name();
	cout << typeid(dizi[4]).name();
	return 1;
}

