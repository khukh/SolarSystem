#pragma once
template <int SIZE>
class Vect
{
public:

	double vect[SIZE];
	int size = SIZE;

	friend Vect<SIZE> operator + (const Vect<SIZE>& a, const Vect<SIZE>& b) {
		Vect<SIZE> c;
		for (int i = 0; i < SIZE; i++) {
			c.vect[i] = a.vect[i] + b.vect[i];
		}
		return c;
	}

	friend Vect<SIZE> operator * (const Vect<SIZE>& a, const Vect<SIZE>& b) {
		Vect<SIZE> c;
		for (int i = 0; i < SIZE; i++) {
			c.vect[i] = a.vect[(i + 1) % 3] * b.vect[(i + 2) % 3] - b.vect[(i + 1) % 3] * a.vect[(i + 2) % 3];
		}
		return c;
	}

	friend Vect<SIZE> operator * (const Vect<SIZE>& a, double b) {
		Vect<SIZE> c;
		for (int i = 0; i < SIZE; i++) {
			c.vect[i] = a.vect[i] * b;
		}
		return c;
	}
	Vect<SIZE>& operator=(const Vect<SIZE>& right) {
		//проверка на самоприсваивание
		if (this == &right) {
			return *this;
		}
		for (int i = 0; i < SIZE; i++) {
			vect[i] = right.vect[i];
		}

		return *this;
	}

	friend Vect<SIZE> operator*(matrix &A, Vect<SIZE> &b) {
		Vect<3> res;


		for (int i = 0; i < 3; i++) {
			double r1 = 0;
			for (int j = 0; j < 3; j++) {
				double a = A.matr[i][j];
				double r = b.vect[j];

				r1 += a * r;

			}
			res.vect[i] = r1;

		}
		return res;
	}
};