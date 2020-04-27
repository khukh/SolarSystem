#pragma once
#ifndef CONSTANTS
#define CONSTANTS
const double PI = 3.14159265358979323846264338328;
const double toDeg = 180 / PI;
const double toRad = PI / 180;
//гравитационные постоянные
const double MU_SUN = 132712440018E+9;

//начальные условия положения планет
const Vect <7> initApophis = { -2.785344504736421E+03, -2.575857680062690E+04,  1.305309421257014E+03,-1.596306607427476E+11,  2.929153417513816E+10, -5.352392510420750E+09, 0 };
const double MU_APOPHIS = 2;


/////////////
const int N = 13;

const double R_P = 500; // м
const double R_A = R_P + 5 * N;

const double MU_APO = 2;//m^3/c^2

const double INC = PI / 2;
const double RAAN = 0;
const double W = 0;
const double U0 = 0;

const double V_X0 = 0;
const double V_Y0 = 0;
const double V_Z0 = sqrt(2 * MU_APO / R_P) * sqrt(R_A / (R_A + R_P));


const double X0 = R_P;
const double Y0 = 0;
const double Z0 = 0;

const double H = 10;

#endif // !CONSTANTS

