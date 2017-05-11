#pragma once
#include "nr3.h"
#include "utilities.h"
#include <iostream>

using namespace std;
using namespace util;


#define D   0.001   // kg/m
#define M   0.058   // kg
#define G   9.81


#pragma region Print functions
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
	//cout << scientific;
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
#pragma endregion 


template<class T>
void midpoint(T &func, Doub a, Doub b, Doub acc, VecDoub &yi)
{
	//vector for print
	vector<Doub> indexV;
	vector<Doub> nV;
	vector<Doub> resultV1;
	vector<Doub> resultV2;
	vector<Doub> resultV3;
	vector<Doub> resultV4;
	vector<Doub> alphaKV;
	vector<Doub> errorV;
	vector<string> headers = { "N" , "x1'", "x2'", "y1'", "y2'", "AlphaK", "Error" };

	bool running = true;
	Doub N = 2;
	VecDoub initialConditions = yi;
	
	Doub lastErr = 0;
	VecDoub xmid(yi.size());
	VecDoub dxdy(yi.size());
	int iter = 0;
	Doub h;
	MatDoub res(yi.size(), 20, 0.0);
	VecDoub temp(yi.size(), 0.0);
	VecDoub temp2(yi.size(), 0.0);

	//insert initial Values in Res
	for (int i = 0; i < yi.size(); i++)
		res[i][iter] = yi[i];

	//repeat calcs until integral stabilizes
	while (running)
	{
		h = (b - a) / N;

		for (int i = 0; i <= N; i++)
		{
			func(yi, dxdy);
			for (int j = 0; j < yi.size(); j++)
				xmid[j] = yi[j] + 0.5 * dxdy[j] * h;

			func(xmid, dxdy);
			for (int j = 0; j < yi.size(); j++)
			{
				res[j][iter+1] = res[j][iter] + dxdy[j] * h;
				yi[j] = res[j][iter + 1];
			}
		}

		nV.push_back(N);
		resultV1.push_back(res[0][iter]);
		resultV2.push_back(res[1][iter]);
		resultV3.push_back(res[2][iter]);
		resultV4.push_back(res[3][iter]);

		//print saving
		Doub alphaK, error;
		if (iter > 1)
		{
			for (int i = 0; i < temp.size(); i++)
			{
				temp[i] = res[i][iter - 2] - res[i][iter - 1];
				temp2[i] = res[i][iter - 1] - res[i][iter];
			}
			alphaK = temp.length() / temp2.length();
			//order will always be 2 if midpoint is implemented correct, so alphaK = 4
			for (int i = 0; i < temp.size(); i++)
			{
				temp[i] = abs(res[i][iter] - res[i][iter-1]);
			}
			error = (temp.length()) / static_cast<Doub>(4 - 1);
		}
		else
		{
			alphaK = 0;
			error = 0;
		}
		
		alphaKV.push_back(alphaK);
		errorV.push_back(error);

		//finish check
		if (abs(lastErr - error) <= acc && iter > 2)
		{
			//print stuff
			vector<vector<Doub>> printV;
			printV.push_back(nV);
			printV.push_back(resultV1);
			printV.push_back(resultV2);
			printV.push_back(resultV3);
			printV.push_back(resultV4);
			printV.push_back(alphaKV);
			printV.push_back(errorV);
			printTable(printV, headers);
			return;
		}

		//reset values for next iteration
		lastErr = error;
		N = N * 2;
		iter++;

		//reset initial values
		yi = initialConditions;
	}
}

struct rhs {
	void operator() (VecDoub_I &y, VecDoub_O &dydx) {
		dydx[0] = y[1];
		dydx[1] = -(D / M)*y[1] * sqrt(pow(y[1], 2) + pow(y[3], 2));
		dydx[2] = y[3];
		dydx[3] = -G - (D / M)*y[3] * sqrt(pow(y[1], 2) + pow(y[3], 2));
	}
};

int main()
{
	
	Doub acc = pow(10, -3);
	rhs rhs;
	VecDoub yi(4);
	yi[0] = 0;
	yi[1] = 1;
	yi[2] = 0;
	yi[3] = 3;
	midpoint(rhs, 0.0, 1.0, acc, yi);

	yi.print();



	system("pause");
	return 0;
}