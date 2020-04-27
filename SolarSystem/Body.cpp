#include "pch.h"
#include "Body.h"

#include <fstream>
//#include "Constants.h"

Body::Body( double muBodyy, Vect<7> initCoord, Body *centralBody)
{
	centBody = centralBody;
	mu = (centralBody != NULL) ? (muBodyy + centralBody->mu):(muBodyy);
	
	muBody = muBodyy;
	parametr = initCoord;
	co = { parametr.vect[3], parametr.vect[4], parametr.vect[5] };
	ve = { parametr.vect[0], parametr.vect[1], parametr.vect[2] };
	if (centralBody != NULL) {
		DKE();
	}
	
}



//Body::Body(double muBodyy, Vect<7> initCoord)
//{
//	//centBody = centralBody;
//	//mu = muBodyy + centralBody->mu;
//	muBody = muBodyy;
//	parametr = initCoord;
//	co = { parametr.vect[3], parametr.vect[4], parametr.vect[5] };
//	ve = { parametr.vect[0], parametr.vect[1], parametr.vect[2] };
//	DKE();
//}

Body::~Body()
{
}

Vect<3> Body::getPosition(double t)
{
	if (centBody != NULL) {
		EDK(t);
	}
	
	return co;
}

double Body::getMu()
{
	return muBody;
}


void Body::DKE()
{


	double r = sqrt(parametr.vect[3] * parametr.vect[3] + parametr.vect[4] * parametr.vect[4] + parametr.vect[5] * parametr.vect[5]);
	

	double vSq = parametr.vect[0] * parametr.vect[0] + parametr.vect[1] * parametr.vect[1] + parametr.vect[2] * parametr.vect[2];//квадрат скороси
	h = vSq - 2 * mu / r;
	c = co * ve;
	cAbs = sqrt(c.vect[0] * c.vect[0] + c.vect[1] * c.vect[1] + c.vect[2] * c.vect[2]);
	f = ve * c + co * (-mu / r);
	fAbs = sqrt(mu*mu + cAbs * cAbs*h);
	/////////////////
	cepElem.vect[0] = acos(c.vect[2] / cAbs);//inc
	double cxy = sqrt(c.vect[0] * c.vect[0] + c.vect[1] * c.vect[1]);
	cepElem.vect[1] = atan2(c.vect[0] / cxy, -c.vect[1] / cxy); //raan
	p = cAbs * cAbs / mu;
	cepElem.vect[2] = -mu / h; //a
	cepElem.vect[3] = fAbs / mu; //ecc
	Vect<3> q = { cos(cepElem.vect[1]),sin(cepElem.vect[1]),0 };
	//w
	double cosW = (q.vect[0] * f.vect[0] + q.vect[1] * f.vect[1] + q.vect[2] * f.vect[2]) / (fAbs); //cos(w)
	double sinW = f.vect[2] / (fAbs * sin(cepElem.vect[0])); //sin(w)
	cepElem.vect[4] = atan2(sinW, cosW); //w

	//u
	double cosU = (q.vect[0] * co.vect[0] + q.vect[1] * co.vect[1] + q.vect[2] * co.vect[2]) / (r); //cos(u)
	double sinU = co.vect[2] / (r * sin(cepElem.vect[0])); //sin(u)
	u = atan2(sinU, cosU);
	u = (u < 0) ? (2 * PI + u) : u;

	//tetta
	double vR = (ve.vect[0] * co.vect[0] + ve.vect[1] * co.vect[1] + ve.vect[2] * co.vect[2]) / (r);
	double cosTetta = (p / r - 1) / cepElem.vect[3]; //cos(tetta)
	double sinTetta = vR * sqrt(p / mu) / cepElem.vect[3]; //sin(tetta)
	tetta = atan2(sinTetta, cosTetta);
	tetta = (tetta < 0) ? (2 * PI + tetta) : tetta;


	//////////////////////////////////


	if (h < 0) {
		Ea = 2 * atan2(sqrt((1 - cepElem.vect[3]) / (1 + cepElem.vect[3]))*tan(0.5*tetta), 1);
		//Ea = atan2(sqrt(1 - ecc * ecc) * sin(tetta) / (1 + ecc * cos(tetta)), ((cos(tetta) + ecc) / (1 + ecc * cos(tetta))));
		Ea = (Ea < 0) ? (2 * PI + Ea) : Ea;
		cepElem.vect[5] = parametr.vect[6] - sqrt(cepElem.vect[2] * cepElem.vect[2] * cepElem.vect[2])*(Ea - cepElem.vect[3] * sin(Ea)) / sqrt(mu);
	}
	else {
		double arg = sqrt((-1 + cepElem.vect[3]) / (1 + cepElem.vect[3]))*tan(0.5*tetta);
		Ea = 2 * atanh(arg);
		//Ea = atan2(sqrt(1 - ecc * ecc) * sin(tetta) / (1 + ecc * cos(tetta)), ((cos(tetta) + ecc) / (1 + ecc * cos(tetta))));
	//	Ea = (Ea < 0) ? (2 * PI + Ea) : Ea;
		cepElem.vect[5] = parametr.vect[6] - sqrt(abs(cepElem.vect[2] * cepElem.vect[2] * cepElem.vect[2]))*(-Ea + cepElem.vect[3] * sinh(Ea)) / sqrt(mu);
	}

}

void Body::EDK(double t)
{
	
	double nM = sqrt(mu) / abs(cepElem.vect[2]) / sqrt(abs(cepElem.vect[2]));
	double tetta1;
	if (cepElem.vect[3] < 1) {
		double Ma = nM * (t - cepElem.vect[5]);
		//double E0 = Ma + ecc * sin(Ma);

		double E1 = rootMA(Ma);
		tetta1 = 2 * atan2(sqrt((1 + cepElem.vect[3]) / (1 - cepElem.vect[3]))*tan(0.5*E1), 1);
	}
	else {
		double Ma = nM * (t - cepElem.vect[5]);
		double H1 = rootMAH(Ma);
		tetta1 = atan2(sqrt(cepElem.vect[3] * cepElem.vect[3] - 1)*sinh(H1) / (cepElem.vect[3] *cosh(H1) - 1), (cepElem.vect[3] - cosh(H1)) / (cepElem.vect[3] *cosh(H1) - 1));
	}


	//	tetta1 = (tetta1 < 0) ? (2 * PI + tetta1) : tetta1;
	//d_tetta = tetta1 - tetta;
	u = tetta1 + cepElem.vect[4];
	double r1 = cepElem.vect[2] * (1 - cepElem.vect[3] * cepElem.vect[3]) / (1 + cepElem.vect[3] * cos(tetta1));
	double r00 = cos(cepElem.vect[1]) * cos(u) - sin(u) * sin(cepElem.vect[1]) * cos(cepElem.vect[0]);
	double r01 = sin(cepElem.vect[1]) * cos(u) + sin(u) * cos(cepElem.vect[1]) * cos(cepElem.vect[0]);
	double r02 = sin(u) * sin(cepElem.vect[0]);
	co = { r1*r00, r1*r01, r1*r02 };

	double n00 = -cos(cepElem.vect[1]) * sin(u) - cos(u) * sin(cepElem.vect[1]) * cos(cepElem.vect[0]);
	double n01 = -sin(cepElem.vect[1]) * sin(u) + cos(u) * cos(cepElem.vect[1]) * cos(cepElem.vect[0]);
	double n02 = cos(u) * sin(cepElem.vect[0]);
	double vR = sqrt(mu / p) * cepElem.vect[3] * sin(tetta1);
	double vN = sqrt(mu / p) *(1 + cepElem.vect[3] * cos(tetta1));
	ve = { vR*r00 + vN * n00, vR*r01 + vN * n01, vR*r02 + vN * n02 };
	parametr = { vR*r00 + vN * n00, vR*r01 + vN * n01, vR*r02 + vN * n02, r1*r00, r1*r01, r1*r02, t };
}

double Body::rootMA(double Ma)
{
	double E1 = Ma + cepElem.vect[3] * sin(Ma);
	double E2 = E1;

	do {
		E1 = E2;
		E2 = E1 - (E1 - cepElem.vect[3] * sin(E1) - Ma) / (1 - cepElem.vect[3] * cos(E1));


	} while (abs(E2 - E1) > 1E-10);
	return E2;
}

double Body::rootMAH(double Ma)
{
	long double H1 = 2 * Ma;
	long double H2 = H1;

	do {
		H1 = H2;
		long double siH = (exp(H1) - exp(-H1)) / 2;
		long double coH = (exp(H1) + exp(-H1)) / 2;
		H2 = H1 - (-H1 + cepElem.vect[3] * siH - Ma) / (-1 + cepElem.vect[3] * coH);


	} while (abs(H2 - H1) > 1E-10);
	return H2;
}


void Body::printParam(std::ofstream &fout) {
	fout << parametr.vect[6] << '\t' << std::scientific;

	fout << parametr.vect[0] << '\t' << parametr.vect[1] << '\t' << parametr.vect[2] << '\t';
	fout << parametr.vect[3] << '\t' << parametr.vect[4] << '\t' << parametr.vect[5] << '\t';
	fout << h << '\t';
	fout << c.vect[0] << '\t' << c.vect[1] << '\t' << c.vect[2] << '\t';
	fout << cAbs << '\t';
	fout << f.vect[0] << '\t' << f.vect[1] << '\t' << f.vect[2] << '\t';
	fout << fAbs << '\t';

	//fout << k << '\t';
	//fout << inc * toDeg << '\t';
	//fout << raan * toDeg << '\t';
	//fout << p << '\t';
	//fout << a << '\t';
	//fout << w << '\t';
	//fout << tetta * toDeg << '\t';
	//fout << u * toDeg << '\t';
	//fout << v.vect[0]/*-parametr.vect[0] */ << '\t' << v.vect[1] /*- parametr.vect[1]*/ << '\t' << v.vect[2] /*- parametr.vect[2]*/ << '\t';
	//fout << coord.vect[0] /*- parametr.vect[3]*/ << '\t' << coord.vect[1] /*- parametr.vect[4]*/ << '\t' << coord.vect[2] /*- parametr.vect[5]*/ << '\t';
	//fout << tau << '\t';
	//fout << Ea * toDeg << '\t';
	//fout << nOrb << '\t';
	//fout << d_tetta << '\t';
	//fout << ecc << '\t';
	fout << '\n';

}

Body * Body::getCentralBody()
{
	return centBody;
}


