#pragma once
class Body
{
public:
	Body( double muBody, Vect <7> initCoord, Body *centralBody = NULL );
	//Body(double muBodyy, Vect<7> initCoord);
	~Body();
	Vect <3> getPosition(double t);
	double getMu();
	void printParam(std::ofstream & fout);
	Body *getCentralBody();
	
private:
	Body* centBody;
	void DKE();
	void EDK(double t);
	double rootMA(double Ma);
	double rootMAH(double Ma);


	double muBody;
	double mu;
	Vect<6> cepElem;/*
					cepElem[0] - inclination
					cepElem[1] - raan
					cepElem[2] - a
					cepElem[3] - ecc
					cepElem[4] - w
					cepElem[5] - tau
					*/
	Vect<7> parametr;/*
					 parametr[0] - Vx
					 parametr[1] - Vy
					 parametr[2] - Vz
					 parametr[3] - X
					 parametr[4] - Y
					 parametr[5] - Z
					 parametr[6] - t
					 */

	
	//double m;//�����

	////////////////////��������� ������� � �� ��������/////////////////////////////////////

	double h; //�������� �������
	Vect<3> c; //�������� ��������
	double cAbs; // ������ ��������� ��������
	Vect<3> f; // �������� ��������
	double fAbs; //������ ��������� ��������

	////////////////////�������� ������ (��������� ��������)/////////////////////////////////////
	double p;
	double Ea;
	double tetta, u;
	////////////////////��������� ��������/////////////////////////////////////
	double vAbs;
	double coordAbs;

	Vect <3> co, ve;

	//CentralBody CB;

};

