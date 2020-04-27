
#include "pch.h"

#define _USE_MATH_DEFINES
#include "status.h"

#include "cmath"
#include <fstream>
//#include "Constants.h"

status::status(Body *CB, std::vector <double> &coordinates) {
	
	//скорости
	parametr.vect[0] = coordinates[0];
	parametr.vect[1] = coordinates[1];
	parametr.vect[2] = coordinates[2];
	//кординаты
	parametr.vect[3] = coordinates[3];
	parametr.vect[4] = coordinates[4];
	parametr.vect[5] = coordinates[5];
	parametr.vect[6] = 0;
	nOrb = 0;
	centBody = CB;
	mu = CB->getMu();
	vInf = { 0,0,0 };
	bInf = { 0,0,0 };
}

status::~status() {}

//значения производных
Vect<7> status::rightPart() {
	Vect<7> prir;
	prir.vect[0] = -mu * parametr.vect[3] / r / r / r + dAlphaX;
	prir.vect[1] = -mu * parametr.vect[4] / r / r / r + dAlphaY;
	prir.vect[2] = -mu * parametr.vect[5] / r / r / r + dAlphaZ;

	prir.vect[3] = parametr.vect[0];
	prir.vect[4] = parametr.vect[1];
	prir.vect[5] = parametr.vect[2];

	prir.vect[6] = 1;
	return(prir);

}



void status::nonIntegr() {
	
	Vect <3> co = { parametr.vect[3], parametr.vect[4], parametr.vect[5] };
	Vect <3> ve = { parametr.vect[0], parametr.vect[1], parametr.vect[2] };

	DKE(co, ve);
	Vect <6> q = { inc, ecc, raan, w, a, tau };
	EDK(q);
	nOrb = int(parametr.vect[6] / (2*PI*sqrt(a*a*a / mu)));
	dAlphaX = dAlphaY = dAlphaZ = 0;
	for (int i = 0; i < inf.size(); i++) {
		Vect <3> position;// = inf[i]->getPosition(parametr.vect[6]);
		Vect <3> Rj;
		if (inf[i] == centBody->getCentralBody()) {
			position = centBody->getPosition(parametr.vect[6]); //(-r_j)
			Rj = position * (-1);
		}
		Vect <3> Rj0 = co + Rj*(-1);
		double absRj0 = sqrt(Rj0.vect[0] * Rj0.vect[0] + Rj0.vect[1] * Rj0.vect[1] + Rj0.vect[2] * Rj0.vect[2]);
		double absRj = sqrt(Rj.vect[0] * Rj.vect[0] + Rj.vect[1] * Rj.vect[1] + Rj.vect[2] * Rj.vect[2]);
		double q = (co.vect[0] * (co.vect[0] / absRj - 2 * Rj.vect[0] / absRj) + co.vect[1] * (co.vect[1] / absRj - 2 * Rj.vect[1] / absRj) + co.vect[2] * (co.vect[2] / absRj - 2 * Rj.vect[2] / absRj)) / absRj;
		double f = (3 + 3 * q + q * q)*q / (pow(1 + q, 1.5) + 1);
		dAlphaX -= inf[i]->getMu()*(co.vect[0] / absRj0 + f * Rj.vect[0] / absRj0) / absRj0 / absRj0;
		dAlphaY -= inf[i]->getMu()*(co.vect[1] / absRj0 + f * Rj.vect[1] / absRj0) / absRj0 / absRj0;
		dAlphaZ -= inf[i]->getMu()*(co.vect[2] / absRj0 + f * Rj.vect[2] / absRj0) / absRj0 / absRj0;
	}
}

//вывод параметров
void status::printParam(std::ofstream &fout) {
	fout << parametr.vect[6] << '\t' << std::scientific;
	
	fout << parametr.vect[0] << '\t' << parametr.vect[1] << '\t' << parametr.vect[2] << '\t';
	fout << parametr.vect[3] << '\t' << parametr.vect[4] << '\t' << parametr.vect[5] << '\t';
	fout << h << '\t';
	fout << c.vect[0] << '\t' << c.vect[1] << '\t' << c.vect[2] << '\t';
	fout << cAbs << '\t';
	fout << f.vect[0] << '\t' << f.vect[1] << '\t' << f.vect[2] << '\t';
	fout << fAbs << '\t';
	
	//fout << k << '\t';
	fout << inc * toDeg << '\t';
	fout << raan * toDeg << '\t';
	fout << p << '\t';
	fout << a << '\t';
	fout << w << '\t';
	fout << tetta*toDeg << '\t';
	fout << u*toDeg << '\t';
	fout << v.vect[0]/*-parametr.vect[0] */ << '\t' << v.vect[1] /*- parametr.vect[1]*/ << '\t' << v.vect[2] /*- parametr.vect[2]*/ << '\t';
	fout << coord.vect[0] /*- parametr.vect[3]*/ << '\t' << coord.vect[1] /*- parametr.vect[4]*/ << '\t' << coord.vect[2] /*- parametr.vect[5]*/ << '\t';
	fout << tau << '\t';
	fout << Ea * toDeg << '\t';
	fout << nOrb << '\t';
	fout << d_tetta << '\t';
	fout << ecc << '\t';
	fout << '\n';

}
void status::setParam(Vect<7> b) {
	parametr = b;
}

Vect<7> status::getParam() {
	Vect<7> get = parametr;
	return(get);
}

status& status::operator=(const status& right) {
	//проверка на самоприсваивание
	if (this == &right) {
		return *this;
	}
	parametr = right.parametr;
	m = right.m;
	r = right.r;
	nOrb = 0;
	nonIntegr();

	return *this;
}

void status::DKE(Vect<3> co, Vect<3> ve)
{

	r = sqrt(parametr.vect[3] * parametr.vect[3] + parametr.vect[4] * parametr.vect[4] + parametr.vect[5] * parametr.vect[5]);
	//if (nOrb == 0) {
		/////////////////
	double vSq = parametr.vect[0] * parametr.vect[0] + parametr.vect[1] * parametr.vect[1] + parametr.vect[2] * parametr.vect[2];//квадрат скороси
	h = vSq - 2 * mu / r;
	c = co * ve;
	cAbs = sqrt(c.vect[0] * c.vect[0] + c.vect[1] * c.vect[1] + c.vect[2] * c.vect[2]);
	f = ve * c + co * (-mu / r);
	fAbs = sqrt(mu*mu + cAbs * cAbs*h);
	/////////////////
	inc = acos(c.vect[2] / cAbs);
	double cxy = sqrt(c.vect[0] * c.vect[0] + c.vect[1] * c.vect[1]);
	raan = atan2(c.vect[0] / cxy, -c.vect[1] / cxy);
	p = cAbs * cAbs / mu;
	a = -mu / h; ///
	ecc = fAbs / mu;
	Vect<3> q = { cos(raan),sin(raan),0 };
	//w
	double cosW = (q.vect[0] * f.vect[0] + q.vect[1] * f.vect[1] + q.vect[2] * f.vect[2]) / (fAbs); //cos(w)
	double sinW = f.vect[2] / (fAbs * sin(inc)); //sin(w)
	w = atan2(sinW, cosW);

	//u
	double cosU = (q.vect[0] * co.vect[0] + q.vect[1] * co.vect[1] + q.vect[2] * co.vect[2]) / (r); //cos(u)
	double sinU = co.vect[2] / (r * sin(inc)); //sin(u)
	u = atan2(sinU, cosU);
	u = (u < 0) ? (2 * PI + u) : u;

	//tetta
	double vR = (ve.vect[0] * co.vect[0] + ve.vect[1] * co.vect[1] + ve.vect[2] * co.vect[2]) / (r);
	double cosTetta = (p / r - 1) / ecc; //cos(tetta)
	double sinTetta = vR * sqrt(p / mu) / ecc; //sin(tetta)
	tetta = atan2(sinTetta, cosTetta);
	tetta = (tetta < 0) ? (2 * PI + tetta) : tetta;

	if (h > 0) {
		Vect<3> alpha = c * f;
		
		vInf = (f * (-mu) + alpha*sqrt(h))*(sqrt(h)/fAbs/fAbs);

		bInf = (f * (cAbs*cAbs) + alpha * (mu / sqrt(h)))*(1/fAbs/fAbs);
	}

	//////////////////////////////////
	

	if (h < 0) {
		Ea = 2 * atan2(sqrt((1 - ecc) / (1 + ecc))*tan(0.5*tetta), 1);
		//Ea = atan2(sqrt(1 - ecc * ecc) * sin(tetta) / (1 + ecc * cos(tetta)), ((cos(tetta) + ecc) / (1 + ecc * cos(tetta))));
		Ea = (Ea < 0) ? (2 * PI + Ea) : Ea;
		tau = parametr.vect[6] - sqrt(a*a*a)*(Ea - ecc * sin(Ea)) / sqrt(mu);
	}
	else {
		double arg = sqrt((-1 + ecc) / (1 + ecc))*tan(0.5*tetta);
		Ea = 2 * atanh(arg);
		//Ea = atan2(sqrt(1 - ecc * ecc) * sin(tetta) / (1 + ecc * cos(tetta)), ((cos(tetta) + ecc) / (1 + ecc * cos(tetta))));
	//	Ea = (Ea < 0) ? (2 * PI + Ea) : Ea;
		tau = parametr.vect[6] - sqrt(abs(a*a*a))*(-Ea + ecc * sinh(Ea)) / sqrt(mu);
	}
	

	//////////////////////////
}

void status::EDK(Vect<6> q)
{
	if (parametr.vect[6] > 922989.00) {
		double gfr = 0;
	}

	double nM = sqrt(mu) / abs(a) / sqrt(abs(a));
	double tetta1;
	if (ecc < 1) {
		double Ma = nM * (parametr.vect[6] - tau);
		//double E0 = Ma + ecc * sin(Ma);

		double E1 = rootMA(Ma);
		tetta1 = 2 * atan2(sqrt((1 + ecc) / (1 - ecc))*tan(0.5*E1), 1);
	}
	else {
		double Ma = nM * (parametr.vect[6] - tau);
		double H1 = rootMAH(Ma);
		tetta1 = atan2(sqrt(ecc*ecc - 1)*sinh(H1) / (ecc*cosh(H1) - 1), (ecc - cosh(H1)) / (ecc*cosh(H1) - 1));
	}
	
	
	//	tetta1 = (tetta1 < 0) ? (2 * PI + tetta1) : tetta1;
	d_tetta = tetta1 - tetta;
	double r1 = a * (1 - ecc * ecc) / (1 + ecc * cos(tetta1));
	double r00 = cos(raan) * cos(u) - sin(u) * sin(raan) * cos(inc);
	double r01 = sin(raan) * cos(u) + sin(u) * cos(raan) * cos(inc);
	double r02 = sin(u) * sin(inc);
	coord = { r1*r00, r1*r01, r1*r02 };

	double n00 = -cos(raan) * sin(u) - cos(u) * sin(raan) * cos(inc);
	double n01 = -sin(raan) * sin(u) + cos(u) * cos(raan) * cos(inc);
	double n02 = cos(u) * sin(inc);
	double vR = sqrt(mu / p) * ecc * sin(tetta1);
	double vN = sqrt(mu / p) *(1 + ecc * cos(tetta1));
	v = { vR*r00 + vN * n00, vR*r01 + vN * n01, vR*r02 + vN * n02 };
}

double status::rootMA(double Ma)
{
	double E1 = Ma +ecc*sin(Ma);
	double E2 = E1;
	
	do {
		E1 = E2;
		E2 = E1 - (E1 - ecc * sin(E1) - Ma) / (1 - ecc * cos(E1));
		

	} while (abs(E2 - E1) > 1E-10);
	return E2;
}

double status::rootMAH(double Ma)
{
	long double H1 =2 * Ma;
	long double H2 = H1;

	do {
		H1 = H2;
		long double siH = (exp(H1)-exp(-H1))/2;
		long double coH = (exp(H1) + exp(-H1)) / 2;
		H2 = H1 - (-H1 + ecc * siH - Ma) / (-1 + ecc * coH);


	} while (abs(H2 - H1) > 1E-10);
	return H2;
}

void status::addInfluence(Body *a)
{
	inf.push_back(a);
	double mu = inf[inf.size() - 1]->getMu();
}


