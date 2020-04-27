#pragma once

#include "pch.h"

#include "iostream" 
#include <fstream>




class status {
public:
	status(Body *CB, std::vector <double> &a);
	~status();

	virtual void nonIntegr();		//пересчет неинтегрируемых параметров
	virtual Vect<7> rightPart();	//значени€ производных	
	void addInfluence(Body * a);
	void printParam( std::ofstream &fout);	//вывод параметров
	void setParam(Vect<7> b);
	Vect<7> getParam();
	double r;
	double d_tetta;
	int nOrb;
	status& operator=(const status& right);
protected:
	Body* centBody;
	std::vector <Body*> inf;
	Vect<7> parametr;
	/*	parametr [0] = Vx
		parametr [1] = Vy
		parametr [2] = Vz
		parametr [3] = x
		parametr [4] = y
		parametr [5] = z
		parametr [6] = t
		*/
	void DKE(Vect<3> co, Vect<3> ve);
	void EDK(Vect <6> q);
	double rootMA(double Ma);
	double rootMAH(double Ma);



	double m;//масса

	////////////////////»нтегралы энергии и их проекции/////////////////////////////////////

	double h; //интеграл энергии
	Vect<3> c; //интеграл площадей
	double cAbs; // модуль интеграла площадей
	Vect<3> f; // интеграл Ћапласса
	double fAbs; //модуль интеграла Ћапласса

	////////////////////Ёлементы орбиты ( еплеровы элементы)/////////////////////////////////////
	double inc;
	double raan;
	double p, a;
	double Ea;
	double ecc;
	double w, tetta, u;
	double tau;
	////////////////////ƒекартовы элементы/////////////////////////////////////
	Vect<3> v;
	double vAbs;
	Vect<3> coord;
	double coordAbs;
	double mu;
	Vect<3> vInf, bInf;

	double dAlphaX, dAlphaY, dAlphaZ;
};

