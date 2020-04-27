#pragma once
#include "pch.h"
class matrix {
public:
	matrix(double pitch, double yaw, double roll);
	void fillMatrix(double pitch, double yaw, double roll);
	double matr[3][3];
	//friend std::vector <double> operator*(matrix &A, std::vector <double> &b);
	//std::vector <double> operator * ( matrix &A, std::vector <double> &b);
	~matrix();
};

