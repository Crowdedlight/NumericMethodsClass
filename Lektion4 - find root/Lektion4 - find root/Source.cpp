#pragma once
#include "nr3.h"
#include <iostream>
#include "utilities.h"
#include <math.h>
//#include "roots.h"

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
		if (j == 0) dk.push_back(abs(x1-x2));
		else dk.push_back(dx);	
		//C := xi+1 - x / x - xi-1
		c.push_back((dx*0.5)/dx);

	}
	throw("Too many bisections in rtbis");
}

#pragma endregion 


int main()
{
	x_sub_cosx f_1;
	double x1 = 0;
	double x2 = PI / 2;
	double acc1 = pow(10, -8);

	auto result = rtbis(f_1, x1, x2, acc1);

	system("pause");
	return 0;
}
