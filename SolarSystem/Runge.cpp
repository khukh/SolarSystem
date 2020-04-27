#pragma once

#include "iostream" 
#include "cmath"
#include "pch.h"
#include "iomanip"
#include "status.h"

/*using namespace std;*/

	//действи€ с векторами
std::vector <double> multAndSum(std::vector<double> xi, std::vector<double> k, double C) {
	std::vector <double> b(k.size());
	for (unsigned int i = 0; i <= k.size() - 1; i++) {
		b[i] = xi[i] + k[i] * C;
	}
	return(b);
}

void runge(status &sv, double h) {
	Vect<7> xi;
	xi = sv.getParam();
	const int size = xi.size;
	//sv.setParam(xi);
	//sv.nonIntegr();
	Vect<7> buf;
	Vect<7> k1;
	k1 = sv.rightPart();
	//static std::vector <double> a(xi.size());
	//static std::vector <double> null(xi.size());     //нулевой вектор, т.е. все элементы=0
	k1 = k1 * h;

	sv.setParam(k1 * 0.5 + xi);
	sv.nonIntegr();
	Vect<7> k2;
	k2 = sv.rightPart();
	k2 = k2 * h;

	sv.setParam(k2 * 0.5 + xi);
	sv.nonIntegr();
	Vect<7> k3;
	k3 = sv.rightPart();
	k3 = k3 * h;

	sv.setParam(k3 * 1 + xi);
	sv.nonIntegr();
	Vect<7> k4;
	k4 = sv.rightPart();
	k4 = k4 * h;
	Vect<7> result;

	result = xi + (k1 + (k2 * 2) + (k3 * 2) + k4) * (1.0 / 6.0);
	
	
	sv.setParam(result);
	sv.nonIntegr();
}


	
	//интегрирование методом Ёйлера
	void euler(status &sv, double h) {
	
		Vect<7> xi;
		xi = sv.getParam();

		Vect<7> f;
		f = sv.rightPart();
		sv.setParam(f * h + xi);		
		sv.nonIntegr();
		
	}

	//void setZnachEuler(vector <double> prir, double h) {
	//	parametr[0] = euler(parametr[0], prir[0], h);
	//	parametr[1] = euler(parametr[1], prir[1], h);
	//	parametr[2] = euler(parametr[2], prir[2], h);
	//	parametr[3] = euler(parametr[3], prir[3], h);
	//	parametr[4] = euler(parametr[4], prir[4], h);
	//	parametr[5] = euler(parametr[5], prir[5], h);
	//}