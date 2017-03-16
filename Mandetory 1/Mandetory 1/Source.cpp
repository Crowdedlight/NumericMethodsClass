#include <iostream>
#include "nr3.h"
#include "utilities.h"
#include "svd.h"
#include <math.h>
#include <vector>

using namespace std;
using namespace util;

#define PI 3.14159265

#pragma region decleareFunc
    vector<double> varians(SVD& svd);
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
	for (int i = 0; i < 2*N; i+= 2)
	{
		double theta_1, theta_2, x, y;

        //dataset 1
		d1 >> theta_1;
		d1 >> theta_2;
		d1 >> x;
		d1 >> y;

        //x part
        A[i][0] = 1;
        A[i][1] = 0;
        A[i][2] = 1 * cos(theta_1);
        A[i][3] = 1 * cos(theta_1 + theta_2);
        z[i] = x;

        //y part
        A[i+1][0] = 0;
        A[i+1][1] = 1;
        A[i+1][2] = 1 * sin(theta_1);
        A[i+1][3] = 1 * sin(theta_1 + theta_2);
        z[i+1] = y;
	}

    for (int i = 0; i < 2*N; i+=2)
    {
        double theta_1, theta_2, x, y;

        //dataset 2
        d2 >> theta_1;
        d2 >> theta_2;
        d2 >> x;
        d2 >> y;

        //x part
        A2[i][0] = 1;
        A2[i][1] = 0;
        A2[i][2] = 1 * cos(theta_1);
        A2[i][3] = 1 * cos(theta_1 + theta_2);
        z2[i] = x;

        //y part
        A2[i+1][0] = 0;
        A2[i+1][1] = 1;
        A2[i+1][2] = 1 * sin(theta_1);
        A2[i+1][3] = 1 * sin(theta_1 + theta_2);
        z2[i+1] = y;
    }

#pragma endregion 

#pragma region computation

	SVD svd1(A);
    SVD svd2(A2);

    cout << "V dataset 1:" << endl;
    cout << "------------------------------------------------------" << endl;
    svd1.v.print();
    cout << endl << endl;

    cout << "W dataset 1:" << endl;
    cout << "------------------------------------------------------" << endl;
    diag(svd1.w).print();
    cout << endl << endl;

    cout << "V dataset 2:" << endl;
    cout << "------------------------------------------------------" << endl;
    svd2.v.print();
    cout << endl << endl;

    cout << "W dataset 2:" << endl;
    cout << "------------------------------------------------------" << endl;
    diag(svd2.w).print();
    cout << endl << endl;

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

#pragma endregion 

#pragma region Residual Error
    //calculate residual error ||aq-z||

    //dataset 1
    auto rErr = (A*x_1 - z).length();
    //dataset 2
    auto rErr2 = (A2*x_2 - z2).length();

#pragma endregion

#pragma region Varians

    //dataset 1
    auto sigma1 = varians(svd1);
    //dataset 2
    auto sigma2 = varians(svd2);


#pragma endregion

#pragma region printout
    cout << "Results for dataset 1:" << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "q: ";    x_1.print();
    cout << endl;
    cout << "Residul Error: " << rErr << endl;
    cout << endl;
    cout << "Singular Values: "; svd1.w.print();
    cout << endl;
    cout << "Varians (x0,y0,a,b): " << sigma1[0] << "    " << sigma1[1] << "    " << sigma1[2] << "    " << sigma1[3] << endl;
    cout << endl << endl;

    cout << "Results for dataset 2:" << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "q: ";    x_2.print();
    cout << endl;
    cout << "Residul Error: " << rErr2 << endl;
    cout << endl;
    cout << "Singular Values: "; svd2.w.print();
    cout << endl;
    cout << "Varians (x0,y0,a,b): " << sigma2[0] << "    " << sigma2[1] << "    " << sigma2[2] << "    " << sigma2[3] << endl;
    cout << endl << endl;

#pragma endregion

    cout << endl << endl;
	system("pause");
	return 0;
}


vector<double> varians(SVD& svd)
{
    vector<double> Res;
    double tempRes = 0;

    for (int j = 0; j < svd.w.size(); j++)
    {        
        tempRes = 0;
        for (int i = 0; i < svd.v.nrows(); i++)
        {
            tempRes += pow((svd.v[i][j] / svd.w[i]), 2);
        }
        Res.push_back(tempRes);
    }
    return Res;
}