
#include <iostream>
#include <fstream>
#include "nr3.h"
#include "utilities.h"
#include "ludcmp.h"
#include "cholesky.h"

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
			A[i][j] = pow(xPont[i],j);
		}
	}

	MatDoub A2(xFilip.size(), 11);
	for (auto i = 0; i < xFilip.size(); i++)
	{
		for (auto j = 0; j < 11; j++)
		{
			A2[i][j] = pow(abs(xFilip[i]), j);
		}
	}
	//printMatrix(B);
#pragma endregion 

#pragma region calculate

	//(A^t * A) * a = (A^t * b)
	auto Atrans = Transpose(A);
	auto Asymm = Atrans * A;
	auto rhs = Atrans * yPont;
	LUdcmp lucomp(Asymm);

	VecDoub x(3);
	lucomp.solve(rhs, x);
	printMatrix(x);

	Cholesky csky(Asymm);
	VecDoub xC(3);
	csky.solve(rhs, xC);
	printMatrix(xC);

#pragma endregion 

	system("pause");
	return 0;
}