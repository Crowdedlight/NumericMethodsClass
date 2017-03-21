#include <iostream>
#include "nr3.h"
#include "utilities.h"
#include "svd.h"

using namespace std;
using namespace util;

int main()
{
	
#pragma region SetupMatrix

	int N = 500;
	MatDoub A(2 * N, 4);	
	VecDoub z(2 * N);

	ifstream d1("d1");
	for (int i = 0; i < 2*N; i++)
	{
		double theta_1, theta_2, x, y;
		d1 >> theta_1;
		d1 >> theta_2;
		d1 >> x;
		d1 >> y;

		if (!(i % 2))
		{
			A[i][0] = 1;
			A[i][1] = 0;
			A[i][2] = 1 * cos(theta_1);
			A[i][3] = 1 * cos(theta_1 + theta_2);
			z[i] = x;
		}
		else
		{
			A[i][0] = 0;
			A[i][1] = 1;
			A[i][2] = 1 * sin(theta_1);
			A[i][3] = 1 * sin(theta_1 + theta_2);
			z[i] = y;
		}

		//rhs
		//cout << A[i][0] << "," << A[i][1] << "," << A[i][2] << "," << A[i][3] << endl;
	}


#pragma endregion 

#pragma region computation

	SVD svd(A);


	//Psudoinverse as we have more equations than unknowns
	printMatrix(svd.w);
	MatDoub w1 = diag(svd.w);
	for (int i = 0; i < w1.nrows(); i++)
	{
		w1[i][i] = 1 / w1[i][i];
	}

	auto x_1 = svd.v * w1 * svd.u.transpose() * z;

	//x_1.print();

#pragma endregion 

	system("pause");
	return 0;
}