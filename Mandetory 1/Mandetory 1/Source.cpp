#include <iostream>
#include "nr3.h"
#include "utilities.h"
#include "svd.h"

using namespace std;
using namespace util;

#pragma region decleareFunc
    double varians(int j, SVD& svd);
#pragma endregion

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

    //svd1.w.print();


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

    //cout << "-------------------------------" << endl;
    //cout << "-------- Paremeters d1 --------" << endl;
    //cout << "-------------------------------" << endl;
	//x_1.print();
    //cout << "-------------------------------" << endl;
    //cout << "-------- Paremeters d2 --------" << endl;
    //cout << "-------------------------------" << endl;
    //x_2.print();

#pragma endregion 

#pragma region Residual Error
    //calculate residual error ||aq-z||

    //dataset 1
    auto rErr = (A*x_1 - z).length();
    //dataset 2
    auto rErr2 = (A2*x_2 - z2).length();


    //cout << "-------------------------------" << endl;
    //cout << "------ Residual error d1 ------" << endl;
    //cout << "-------------------------------" << endl;
    //cout << rErr << endl;
    //
    //cout << "-------------------------------" << endl;
    //cout << "------ Residual error d2 ------" << endl;
    //cout << "-------------------------------" << endl;
    //cout << rErr2 << endl;

#pragma endregion

#pragma region Varians

    //dataset 1
    auto sigma1_x = varians(0, svd1);
    auto sigma1_y = varians(1, svd1);
    auto sigma1_a = varians(2, svd1);
    auto sigma1_b = varians(3, svd1);

    auto sigma2_x = varians(0, svd2);
    auto sigma2_y = varians(1, svd2);
    auto sigma2_a = varians(2, svd2);
    auto sigma2_b = varians(3, svd2);

    //cout << "-------------------------------" << endl;
    //cout << "---------- Varians d1 ---------" << endl;
    //cout << "-------------------------------" << endl;
    //cout << sigma1_x << ", " << sigma1_y << ", " << sigma1_a << ", " << sigma1_b << endl;
    //
    //cout << "-------------------------------" << endl;
    //cout << "---------- Varians d2 ---------" << endl;
    //cout << "-------------------------------" << endl;
    //cout << sigma2_x << ", " << sigma2_y << ", " << sigma2_a << ", " << sigma2_b << endl;

#pragma endregion


    cout << "Results for dataset 1:" << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "q: ";    x_1.print();
    cout << endl;
    cout << "Residul Error: " << rErr << endl;
    cout << endl;
    cout << "Singular Values: "; svd1.w.print();
    cout << endl;
    cout << "Varians (x0,y0,a,b): " << sigma1_x << "    " << sigma1_y << "    " << sigma1_a << "    " << sigma1_b << endl;
    cout << endl << endl;

    cout << "Results for dataset 2:" << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "q: ";    x_2.print();
    cout << endl;
    cout << "Residul Error: " << rErr2 << endl;
    cout << endl;
    cout << "Singular Values: "; svd2.w.print();
    cout << endl;
    cout << "Varians (x0,y0,a,b): " << sigma2_x << "    " << sigma2_y << "    " << sigma2_a << "    " << sigma2_b << endl;
    cout << endl << endl;


    cout << endl << endl;
	system("pause");
	return 0;
}


double varians(int j, SVD& svd)
{
    double tempRes = 0;

    for (int i = 0; i < svd.w.size(); i++)
    {
        tempRes += pow((svd.v[j][i] / svd.w[i]), 2);
    }

    return tempRes;
}