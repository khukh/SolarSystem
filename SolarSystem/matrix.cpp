#include "pch.h"
#include "matrix.h"


matrix::matrix(double pitch, double yaw, double roll) { //из связанной в стартовую
	//matr.assign(3, std::vector<double>(3));

	matr[0][0] = cos(pitch)*cos(yaw);
	matr[0][1] = -sin(pitch)*cos(yaw)*cos(roll) + sin(yaw)*sin(roll);
	matr[0][2] = sin(pitch)*cos(yaw)*sin(roll)+sin(yaw)*cos(roll);

	matr[1][0] = sin(pitch);
	matr[1][1] = cos(pitch)*cos(roll);
	matr[1][2] = -cos(pitch)*sin(roll);

	matr[2][0] = -cos(pitch)*sin(yaw);
	matr[2][1] = sin(pitch)*sin(yaw)*cos(roll)+cos(yaw)*sin(roll);
	matr[2][2] = -sin(pitch)*sin(yaw)*sin(roll) + cos(yaw)*cos(roll);
}

void matrix::fillMatrix(double pitch, double yaw, double roll)
{
	matr[0][0] = cos(pitch)*cos(yaw);
	matr[1][0] = -sin(pitch)*cos(yaw)*cos(roll) + sin(yaw)*sin(roll);
	matr[2][0] = sin(pitch)*cos(yaw)*sin(roll) + sin(yaw)*cos(roll);

	matr[0][1] = sin(pitch);
	matr[1][1] = cos(pitch)*cos(roll);
	matr[2][1] = -cos(pitch)*sin(roll);

	matr[0][2] = -cos(pitch)*sin(yaw);
	matr[1][2] = sin(pitch)*sin(yaw)*cos(roll) + cos(yaw)*sin(roll);
	matr[2][2] = -sin(pitch)*sin(yaw)*sin(roll) + cos(yaw)*cos(roll);
}




matrix::~matrix() {}

//std::vector<double> operator*(matrix &A, std::vector<double> &b) {
//	std::vector<double> res(3);
//	
//	std::array <double, 3> arrdd= {1, 2, 3} ;
//	std::array <std::array <double, 3>, 3> arrdd2;
//
//	Vect <3> vv= { 1,2,3 };
//	double probe1[3];
//	double probe2[3][3];
//	for (int i = 0; i < 3; i++) {
//		for (int j = 0; j < 3; j++) {
//			double a = A.matr[i][j];
//			double r = b[j];
//
//			res[i] += a * r;
//
//			double p1 = probe1[i];
//			double p2 = probe2[i][j];
//			double p3 = res[i];
//			double p4 = arrdd[i];
//			double p5 = arrdd2[i][j];
//			double p6 = vv.vect[i];
//		}
//
//	}
//	return res;
//}
