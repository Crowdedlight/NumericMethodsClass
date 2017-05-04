#pragma once

#include "nr3.h"
#include "utilities.h"
#include <iostream>

using namespace std;
using namespace util;

struct rhs {
	Doub operator() (const Doub i, const Doub j, const Doub h) {
		Doub Xi = i*h;
		Doub Yj = j*h;
		return (1 + Xi + Yj);
	}
};

MatDoub buildAMatrix(int n, rhs &ro, Doub func_Rand, MatDoub &A, VecDoub &B)
{
	Doub h = 1 / n;

	for (auto i = 1; i <= n-1; i++)
	{
		for (auto j = 1; j <= n-1; j++)
		{
			int alpha = i - 1 + (n - 1)*(j - 1);
			A[alpha][alpha] = -4;
			B[alpha] = pow(h, 2) * ro(i,j,h);

			if (i > 1 && i <= n)
			{
				A[alpha - 1][alpha] = 1; 
			}

			if(i >= 1 && i < n-1)
			{
				A[alpha + 1][alpha] = 1;
			}

			if (j > 1 && j <= n)
			{
				A[alpha][alpha - (n-1)] = 1;
			}

			if (j >= 1 && j < n-1)
			{
				A[alpha][alpha + (n-1)] = 1;
			}

		}
	}
	return A;
}

int main()
{
	rhs f1;
	int n = 4;
	int Asize = pow(n - 1, 2);
	MatDoub A(Asize, Asize, 0.0);
	VecDoub B(Asize, 0.0);
	buildAMatrix(n, f1, 0, A, B);

	A.print();
	cout << endl << endl;
	B.print();


	system("pause");
	return 0;
}