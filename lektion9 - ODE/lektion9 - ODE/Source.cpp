#pragma once
#include "nr3.h"
#include <iostream>
#include "utilities.h"

using namespace std;
using namespace util;

#pragma region HEADERS
template<class T>
void Euler(VecDoub_IO &x, T &func, int x_max, int n);
#pragma endregion 

#pragma region Print functions
void printTable(vector<vector<Doub>> colums, vector<string> headers)
{
	//print headers
	for (int i = 0; i < headers.size(); i++)
	{
		if (i == 0)
			cout << setw(10) << headers[i] << setw(3) << "|";
		else
			cout << setw(15) << headers[i] << setw(3) << "|";
	}
	cout << endl;
	for (int j = 0; j < (13 + (18 * (headers.size() - 1))); j++) { cout << "-"; };
	cout << endl;

	//print values
	//cout << scientific;
	//loop though all rows
	for (int i = 0; i < colums[0].size(); i++)
	{
		//loop though all colums vectors
		for (int j = 0; j < colums.size(); j++)
		{
			//index
			if (j == 0)
				cout << setw(10) << int(colums[j][i]) << setw(3) << "|";
			else
				cout << setw(15) << colums[j][i] << setw(3) << "|";
		}
		cout << endl;
	}

	//End spacing 
	cout << endl << endl << endl;
}
#pragma endregion 

VecDoub vecFunc(VecDoub_I x)
{
	VecDoub f(2);
	f[0] = x[0]*x[1];			//y1
	f[1] = -pow(x[0],2);		//y2

	return f;
};

int main()
{
	
#pragma region ODE manually implemented

	//guess variables
	//VecDoub_IO x(2);
	//x[0] = 1;			// y1
	//x[1] = 1;			// y2

	//Euler(x, vecFunc, 14, 20);
	//cout << x[0] << endl;
	//cout << x[1] << endl;

#pragma endregion 



	system("pause");
	return 0;
}


#pragma region METHODS

template<class T>
void Euler(VecDoub_IO &x, T &func, int x_max, int n)
{
	Doub h = Doub(x_max) / n;
	VecDoub_IO res(x.size());

	vector<double> vN;
	vector<double> vRes1;
	vector<double> vRes2;
	vector<double> vConstant;
	//vector<double> vOrder;
	//vector<double> vError;
	vector<string> vHeader = { "N", "y1", "y2", "Constant" };

	//initial guess
	vN.push_back(0);
	vRes1.push_back(1);
	vRes2.push_back(1);
	vConstant.push_back(pow(x[0],2) + pow(x[1],2));
	int iter = 0;

	while(true)
	{
		h = Doub(x_max) / n;
		
		//calculate res for this n
		for(int i = 1; i <= x_max; i++)
		{
			res = x + (func(x) * h);
		}
		x = res;
		
		//vectors
		vN.push_back(n);
		vRes1.push_back(res[0]);
		vRes2.push_back(res[1]);

		//constant
		vConstant.push_back(pow(x[0],2) + pow(x[1],2));

		//finish check 
		if (iter == 10)
		{
			//richardson error on y(20)


			vector<vector<double>> printVec;
			printVec.push_back(vN);
			printVec.push_back(vRes1);
			printVec.push_back(vRes2);
			printVec.push_back(vConstant);
			//printVec.push_back(vError);
			//printVec.push_back(vOrder);

			printTable(printVec, vHeader);
			break;
		}

		//update n
		n = n*2;
		iter++;
	}
}

#pragma endregion 