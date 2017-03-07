#include <iostream>
#include "nr3.h"
#include "svd.h"
#include "utilities.h"

using namespace std;
using namespace util;

int main()
{
	//Data load
	VecDoub xFilip(82); VecDoub yFilip(82);
	ifstream Filip("FilipData.dat");
	for (int i = 0; i < 82; i++) {
		Filip >> yFilip[i];
		Filip >> xFilip[i];
	}

	VecDoub xPont(40); VecDoub yPont(40);
	ifstream Pont("PontiusData.dat");
	for (int i = 0; i < 40; i++) {
		Pont >> yPont[i];
		Pont >> xPont[i];
	}

#pragma region SetupMatrix
	MatDoub A(xPont.size(), 3);
	for (auto i = 0; i < xPont.size(); i++)
	{
		for (auto j = 0; j < 3; j++)
		{
			A[i][j] = pow(xPont[i], j);
		}
	}

	MatDoub A2(xFilip.size(), 11);
	for (auto i = 0; i < xFilip.size(); i++)
	{
		for (auto j = 0; j < 11; j++)
		{
			A2[i][j] = pow(xFilip[i], j);
		}
	}
	printMatrix(A);
#pragma endregion

	//pontius
	SVD svdP(A);

	MatDoub w1 = diag(svdP.w);
	for (int i = 0; i < w1.nrows(); i++)
	{
		w1[i][i] = 1 / w1[i][i];
	}

	auto x_1 = svdP.v * w1 * svdP.u.transpose() * yPont;

	x_1.print();

	//fillip
	SVD svdF(A2);
	MatDoub w2 = diag(svdF.w);
	for (int i = 0; i < w2.nrows(); i++)
	{
		w2[i][i] = 1 / w2[i][i];
	}

	auto x_2 = svdF.v * w2 * svdF.u.transpose() * yFilip;

	x_2.print();

	system("pause");
	return 0;
}