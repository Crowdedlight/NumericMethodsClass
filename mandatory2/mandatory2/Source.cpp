//#pragma once

#include <iostream>
#include "nr3.h"
#include "utilities.h"
#include <math.h>

using namespace std;
using namespace util;

#pragma region headers

double gamma(double &x, double &y, double &z);
double GAMMA(double t, double &y, double &z);

template <class T>
double crossSection(T &func, double &y, double &z, double acc);
template <class T>
Doub trapez_mod(T &func, Doub a, Doub b, Doub acc, Doub y, Doub z);

#pragma endregion

#pragma region Functions
const double a = 0.1;
const double b = 0.4;
const double v = 1;
const double k = 1.2;
const double w = 0.186;
const double rho = 1200;

double gamma(double &x, double &y, double &z) {
	return pow((pow(x, 2) / (pow(a, 2) * pow(z, 2))) + (pow(y, 2) / (pow(b, 2) * pow(z, 2))) + 1, -0.5);
}
double GAMMA(double t, double &y, double &z) {
	return pow(((pow(v, 2) * pow(t, 2)) / (pow(a, 2) * pow(z, 2))) + (pow(y, 2) / (pow(b, 2) * pow(z, 2))) + 1, -0.5);
}
struct func_gamma {
	func_gamma() {};
	double operator()(double gamma) { return (1 - (pow(1 - (gamma), 2) / (3 / 2 - sqrt(2)))) * cos(k * pow(1 - (gamma), 2)); }
};
struct func_GAMMA {
	func_GAMMA() {};
	double operator()(double GAMMA) { return (1 - (pow(1 - GAMMA, 2) / (3 / 2 - sqrt(2)))) * cos(k * pow(1 - GAMMA, 2)); }
};
struct Alpha {
	Alpha() {};
	double operator()(double y, double z) { return (a*z / v) * sqrt(1 - (pow(y, 2) / (pow(b, 2)*pow(z, 2)))); }
};
#pragma endregion

#pragma region Print functions
void printTable(vector<vector<Doub>> colums, vector<string> headers)
{
	//print headers
	for (int i = 0; i < headers.size(); i++)
	{
		if (i == 0)
			cout << setw(5) << headers[i] << setw(3) << "|";
		else
			cout << setw(15) << headers[i] << setw(3) << "|";
	}
	cout << endl;
	for (int j = 0; j < (8 + (18 * (headers.size() - 1))); j++) { cout << "-"; };
	cout << endl;

	//print values
	//cout << scientific;
	//loop though all rows
	for (int i = 0; i < colums[0].size(); i++)
	{
		//loop though all colums vectors
		for (int j = 0; j < colums.size(); j++)
		{
			//index
			if (j == 0)
				cout << setw(5) << int(colums[j][i]) << setw(3) << "|";
			else
				cout << setw(15) << colums[j][i] << setw(3) << "|";
		}
		cout << endl;
	}

	//End spacing 
	cout << endl << endl << endl;
}
#pragma endregion 

int main()
{
	Doub acc = pow(10, -3);
	Doub Y = 0.1;
	Doub Z = 0.5;

	func_GAMMA fG;

	auto res22 = crossSection(fG, Y, Z, acc);	
	cout << endl << "final result = " << res22 << endl;

	system("pause");
	return 0;
}

#pragma region trapez

template<class T>
Doub trapez_mod(T &func, Doub a_i, Doub b_i, Doub acc, Doub y, Doub z)
{
	Alpha alpha;
	Doub N = round(2 * alpha(y, z) / 0.02);
	Doub h = 0;
	Doub sum = 0;
	Doub res = 0;
	Doub lastErr = 0;
	Doub iter = 0;

	vector<double> vIter;
	vector<double> vRes;
	vector<double> vAlphaK;
	vector<double> vErr;
	vector<double> vAr;
	vector<string> vHeader = { "n", "Result", "alphaK", "Err", "Ar" };

	while (true)
	{
		//calculate h
		Alpha al;
		h = 2 * al(y, z) / N;

		//sum
		for (int i = 1; i < N; i++)
			sum += w / pow(z, 2) * pow(GAMMA((a_i + i * h), y, z), 3) * func(GAMMA((a_i + i * h), y, z));

		//calculate rest
		res = h * (w / pow(z, 2) * pow(GAMMA(a_i, y, z), 3) * func(a_i) / 2 + w / pow(z, 2) * pow(GAMMA(b_i, y, z), 3) * func(b_i) / 2 + sum);
		res = res / rho;

		//Save values
		vIter.push_back(N);
		vRes.push_back(res);

		Doub alphaK, error;
		if (iter > 1)
		{
			alphaK = (vRes[iter - 2] - vRes[iter - 1]) / (vRes[iter - 1] - vRes[iter]);
			error = abs((alphaK*vRes[iter] - vRes[iter - 1]) / (alphaK - 1) - vRes[iter]);
		}
		else
		{
			alphaK = 0;
			error = 0;
		}

		vAlphaK.push_back(alphaK);
		vErr.push_back(error);

		//finish check
		if (abs(lastErr - error) <= acc && iter > 1)
		{
			//Ricardson exterpolartion
			int aK = round(vAlphaK.back());
			vAr.push_back(0);
			for (auto i = 1; i < vRes.size(); i++)
			{
				auto Ar = (aK * vRes[i] - vRes[i - 1]) / (aK - 1);
				vAr.push_back(Ar);
			}

			vector<vector<Doub>> print;
			print.push_back(vIter);
			print.push_back(vRes);
			print.push_back(vAlphaK);
			print.push_back(vErr);
			print.push_back(vAr);
			printTable(print, vHeader);

			return res;
		}

		//reset values for next iteration
		lastErr = error;
		sum = 0;
		res = 0;
		N = N * 2;
		iter++;
	}
}

#pragma endregion

template <class T>
double crossSection(T &func, double &y, double &z, double acc) {
	if (abs(y) > b*z)
		return 0;
	Alpha al;
	double a_integralBorder = -al(y, z);
	double b_integralBorder = al(y, z);

	//double res = integralTrapz_mod(func, a_integralBorder, b_integralBorder, N, y, z, h);
	double res = trapez_mod(func, a_integralBorder, b_integralBorder, acc, y, z);
	return res;
}
