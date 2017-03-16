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

    MatDoub A2(2 * N, 4);
    VecDoub z2(2 * N);

	ifstream d1("d1");
    ifstream d2("d2");
	for (int i = 0; i < 2*N; i++)
	{
		double theta_1, theta_2, x, y;

        //dataset 1
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

        //dataset 2
        d2 >> theta_1;
        d2 >> theta_2;
        d2 >> x;
        d2 >> y;

        if (!(i % 2))
        {
            A2[i][0] = 1;
            A2[i][1] = 0;
            A2[i][2] = 1 * cos(theta_1);
            A2[i][3] = 1 * cos(theta_1 + theta_2);
            z2[i] = x;
        }
        else
        {
            A2[i][0] = 0;
            A2[i][1] = 1;
            A2[i][2] = 1 * sin(theta_1);
            A2[i][3] = 1 * sin(theta_1 + theta_2);
            z2[i] = y;
        }
		//cout << A[i][0] << "," << A[i][1] << "," << A[i][2] << "," << A[i][3] << endl;
	}

#pragma endregion 

#pragma region computation

	SVD svd1(A);
    SVD svd2(A2);


	//Psudoinverse as we have more equations than unknowns
	MatDoub w1 = diag(svd1.w);
    MatDoub w2 = diag(svd2.w);
	for (int i = 0; i < w1.nrows(); i++)
	{
		w1[i][i] = 1 / w1[i][i];
        w2[i][i] = 1 / w2[i][i];
	}

	auto x_1 = svd1.v * w1 * svd1.u.transpose() * z;
    auto x_2 = svd2.v * w2 * svd2.u.transpose() * z2;

    cout << "-------------------------------" << endl;
    cout << "-------- Paremeters d1 --------" << endl;
    cout << "-------------------------------" << endl;
	x_1.print();
    cout << "-------------------------------" << endl;
    cout << "-------- Paremeters d2 --------" << endl;
    cout << "-------------------------------" << endl;
    x_2.print();

#pragma endregion 

    cout << endl << endl;
	system("pause");
	return 0;
}