#pragma once

#include <iostream>
#include "nr3.h"
#include "utilities.h"

using namespace std;
using namespace util;

#pragma region headers

template<class T>
Doub trapez(T &func, Doub a, Doub b, Doub acc);

template <class T>
double crossSection(T &func, T &GAMMA, double &t, double &y, double &z);

#pragma endregion 

#pragma region Functions
const double a = 1;
const double b = 1;
const double v = 1;
const double k = 1;
const double w = 1;
const double rho = 1;

struct gamma {
	gamma() {};
	double operator()(double x, double y, double z) { return pow((pow(x, 2) / (pow(a, 2) * pow(z, 2))) + (pow(y, 2) / (pow(b, 2) * pow(z, 2))) + 1, -0.5); }
};
struct GAMMA {
	GAMMA() {};
	double operator()(double t, double y, double z) { return pow(((pow(v, 2) * pow(t, 2)) / (pow(a, 2) * pow(z, 2))) + (pow(y, 2) / (pow(b, 2) * pow(z, 2))) + 1, -0.5); }
};
struct func_gamma {
	func_gamma() {};
	double operator()(double gamma) { return (1 - (pow(1 - (gamma), 2) / (3 / 2 - sqrt(2)))) * cos(k * pow(1 - (gamma), 2)); }
};
struct func_GAMMA {
	func_GAMMA() {};
	double operator()(double GAMMA) { return (1 - (pow(1 - GAMMA, 2) / (3 / 2 - sqrt(2)))) * cos(k * pow(1 - GAMMA, 2)); }
};
struct alpha {
	alpha() {};
	double operator()(double y, double z) { return (a*z / v) * sqrt(1 - (pow(y, 2) / (pow(b, 2)*pow(z, 2)))); }
};
#pragma endregion 

int main()
{
	Doub acc = pow(10, -3);

	GAMMA G;
	func_GAMMA fG;
	

	system("pause");
	return 0;
}

#pragma region trapez
template<class T>
Doub trapez(T &func, Doub a, Doub b, Doub acc)
{
	Doub N = 1;
	Doub h = 0;
	Doub sum = 0;
	Doub res = 0;
	Doub lastRes = 0;

	while(true)
	{
		//calculate h
		h = (b - a) / N;

		//sum
		for (int i = 1; i < N; i++)
			sum += func(a + i*h);

		//calculate rest
		res = h * (0.5*func(a) + 0.5*func(b) + sum);

		//finish check
		if (abs(lastRes - res) <= acc)
		{
			return res;
		}

		//reset values for next iteration
		lastRes = res;
		sum = 0;
		res = 0;
		N = N * 2;

	}
}

#pragma endregion 

#pragma region crossSection
template <class T>
double crossSection(T &func, T &GAMMA, double &t, double &y, double &z) {
	if (abs(y) > b*z)
		return 0;
	alpha a;
	double a_integralBorder = a(y, z);
	double b_integralBorder = -a(y, z);

	// Call trapetz on: (w/pow(z,2) * pow(GAMMA,3) * func(GAMMA) , a_integralBorder, b_integralborder)
	Doub res = trapez(w / pow(z, 2) * pow(GAMMA, 3) * func(GAMMA), a_integralBorder, b_integralBorder);

	res = res / rho;
	return res;
}
#pragma endregion 