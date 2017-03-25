#pragma once
#include "nr3.h"
#include <iostream>
#include "utilities.h"
#include <math.h>

using namespace std;
using namespace util;

#define PI 3.141592653589793

#pragma region functors

struct x_sub_cosx {
	x_sub_cosx(){}
	double operator()(double x) { return x - cos(x); }
};

#pragma endregion 

#pragma region rootFunctions

///Print table based on colums vectors
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
	for (int j = 0; j < (8+(18*(headers.size()-1))); j++) { cout << "-"; };
	cout << endl;

	//print values
	cout << scientific;
	//loop though all rows
	for (int i = 0; i < colums[0].size(); i++)
	{
		//loop though all colums vectors
		for( int j = 0; j < colums.size(); j++)
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

///Bisection
template <class T>
Doub rtbis(T &func, const Doub x1, const Doub x2, const Doub xacc) {
	const Int JMAX = 50;
	Doub dx, xmid, rtb;
	Doub f = func(x1);
	Doub fmid = func(x2);

	//table vectors colums
	vector<string> headers = { "k", "Xk", "Dk", "C" };
	vector<Doub> indexes;
	vector<Doub> xk;
	vector<Doub> dk;
	vector<Doub> c;
	vector<vector<Doub>> tableColums;

	if (f*fmid >= 0.0) throw("Root must be bracketed for bisection in rtbis");
	rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
	for (Int j = 0; j<JMAX; j++) {
		fmid = func(xmid = rtb + (dx *= 0.5));
		if (fmid <= 0.0) rtb = xmid;


		if (abs(dx) < xacc || fmid == 0.0) 
		{
			//print vetor
			tableColums.push_back(indexes);
			tableColums.push_back(xk);
			tableColums.push_back(dk);
			tableColums.push_back(c);
			printTable(tableColums, headers);
			return rtb;
		}

		//Table values
		indexes.push_back(j);
		//X_k
		xk.push_back(rtb);
		//D_k
		if (j == 0) dk.push_back(0);
		else dk.push_back(xk[j] - xk[j - 1]);
        //C := xi - x / x-1 - xi-2
        if (j < 2)
        {
            c.push_back(0);
        }
        else
        {
            Doub cCal = (xk[j] - xk[j - 1]) / (xk[j - 1] - xk[j - 2]);
            c.push_back((cCal));
        }

	}
	throw("Too many bisections in rtbis");
}

///False Position
template <class T>
Doub rtflsp(T &func, const Doub x1, const Doub x2, const Doub xacc) {
    const Int MAXIT = 30;
    Doub xl, xh, del;
    Doub fl = func(x1);
    Doub fh = func(x2);
    if (fl*fh > 0.0) throw("Root must be bracketed in rtflsp");
    if (fl < 0.0) {
        xl = x1;
        xh = x2;
    }
    else {
        xl = x2;
        xh = x1;
        SWAP(fl, fh);
    }
    Doub dx = xh - xl;

    //table vectors colums
    vector<string> headers = { "k", "Xk", "Dk", "C" };
    vector<Doub> indexes;
    vector<Doub> xk;
    vector<Doub> dk;
    vector<Doub> c;
    vector<vector<Doub>> tableColums;


    for (Int j = 0; j < MAXIT; j++) {
        Doub rtf = xl + dx*fl / (fl - fh);
        Doub f = func(rtf);
        if (f < 0.0) {
            del = xl - rtf;
            xl = rtf;
            fl = f;
        }
        else {
            del = xh - rtf;
            xh = rtf;
            fh = f;
        }
        dx = xh - xl;

        //Table values
        indexes.push_back(j);
        //X_k
        xk.push_back(rtf);
        //D_k
        if (j == 0) dk.push_back(0);
        else dk.push_back(xk[j] - xk[j-1]);
        //C := xi - x / x-1 - xi-2
        if (j < 2)
        {
            c.push_back(0);
        }
        else
        {
            Doub cCal = (xk[j] - xk[j - 1]) / (xk[j - 1] - xk[j - 2]);
            c.push_back((cCal));
        }


        //End check
        if (abs(del) < xacc || f == 0.0)
        {
            //print vetor
            tableColums.push_back(indexes);
            tableColums.push_back(xk);
            tableColums.push_back(dk);
            tableColums.push_back(c);
            printTable(tableColums, headers);
            return rtf;
        }
    }
    throw("Maximum number of iterations exceeded in rtflsp");
}

#pragma endregion 




int main()
{
	x_sub_cosx f_1;
	double x1 = 0;
	double x2 = PI / 2;
	double acc1 = pow(10, -8);

    //Bisection
    cout << "----------------------" << endl;
    cout << "------ BISECTION -----" << endl;
    cout << "----------------------" << endl;
	auto result = rtbis(f_1, x1, x2, acc1);

    //False position - Regular falsi
    double acc2 = pow(10, -16);
    cout << "----------------------" << endl;
    cout << "--- Regular Falsi ----" << endl;
    cout << "----------------------" << endl;
    auto res2 = rtflsp(f_1, x1, x2, acc2);


	system("pause");
	return 0;
}
