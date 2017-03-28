#pragma once
#include <iostream>
#include "nr3.h"
#include "roots_multidim.h"
#include "qrdcmp.h"
#include "ludcmp.h"
#include "utilities.h"

using namespace util;
using namespace std;

Doub k = 2.5;
Doub v = 120;
Doub w = 4.0;
Doub alpha = 2*pow(10, -7);
Doub d;
vector<Doub> n = { 5.0, 2.0, 1.0, 0.5, 0.2, 0.1 };
int iter = 0;

VecDoub vecFunc(VecDoub_I x)
{
	VecDoub f(8);
	f[0] = x[6] * cosh(x[4]/x[6])-1-x[2];				// L0
	f[1] = 2*x[6]*sinh(x[4]/x[6])-x[1];					// L 
	f[2] = 2*x[4]+2*k*cos(x[5])- d;						// p
	f[3] = x[2]+k*sin(x[5])-n[iter];					// phi
	f[4] = sinh(x[4]/x[6])-tan(x[3]);					// x
	f[5] = (1 + (v/(w*x[0])))*tan(x[3])-tan(x[5]);		// theta
	f[6] = x[0]*(1+alpha*x[7])-x[1];					// a
	f[7] = (w*x[0])/(2*sin(x[3])) - x[7];				// h

	return f;
};

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
	for (int j = 0; j < (8 + (18 * (headers.size() - 1))); j++) { cout << "-"; };
	cout << endl;

	//print values
	cout << scientific;
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

///Print table based on colums vectors
void printTableVector(VecDoub data, vector<string> headers)
{
	cout << "-------------------------------------------" << endl;
	for(auto n = 0; n < headers.size(); n++)
	{
		cout << setw(10) << headers[n] << setw(3) << "|" << setw(10) << data[n] << endl;
	}

	cout << "-------------------------------------------" << endl;
	
	//End spacing 
	cout << endl << endl << endl;
}

int main()
{

#pragma region setup vectors
	
	//d guess 1
	d = 30;

	//Start Guess
	VecDoub_IO x(8);
	bool check = false;

	x[0] = 30;		// L0
	x[1] = 30;		// L 
	x[2] = 0.1;		// p
	x[3] = 3.14/6;	// phi
	x[4] = 15;		// x
	x[5] = 3.14/3;	// theta
	x[6] = 40;		// a
	x[7] = 5;		// h

#pragma endregion

	vector<string> headers = { "L0", "L", "p", "phi", "x", "theta", "a", "h" };

	for(iter = 0; iter < n.size(); iter++)
	{
		newt(x, check, vecFunc);

		printTableVector(x, headers);
	}
	system("pause");
	return 0;
}