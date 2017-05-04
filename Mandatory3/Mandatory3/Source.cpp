#pragma once
#include "nr3.h"
#include "utilities.h"
#include <iostream>

using namespace std;
using namespace util;


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
Doub midpoint(T &func, Doub a, Doub b, Doub acc)
{
	//vector for print
	vector<Doub> indexV;
	vector<Doub> nV;
	vector<Doub> resultV;
	vector<Doub> alphaKV;
	vector<Doub> errorV;
	vector<string> headers = { "N" , "Result", "AlphaK", "Error" };

	bool running = true;
	Doub N = 1;
	Doub lastErr = 0;
	Doub res = 0;
	int iter = 0;
	Doub h;

	//repeat calcs until integral stabilizes
	while (running)
	{
		h = (b - a) / N;
		res = 0;
		for (int i = 1; i <= N; i++)
		{
			//sum of f(xi-1/2)
			Doub xi = (a + i*h) - (0.5*h);
			res += func(xi);
		}

		//final res
		res = res*h;

		nV.push_back(N);
		resultV.push_back(res);

		//print saving
		Doub alphaK, error;
		if (iter > 1)
		{
			alphaK = (resultV[iter - 2] - resultV[iter - 1]) / (resultV[iter - 1] - resultV[iter]);
			error = abs((alphaK*resultV[iter] - resultV[iter - 1]) / (alphaK - 1) - resultV[iter]);
		}
		else
		{
			alphaK = 0;
			error = 0;
		}
		
		alphaKV.push_back(alphaK);
		errorV.push_back(error);

		//finish check
		if (abs(lastErr - error) <= acc && iter > 1)
		{
			//print stuff
			vector<vector<Doub>> printV;
			printV.push_back(nV);
			printV.push_back(resultV);
			printV.push_back(alphaKV);
			printV.push_back(errorV);
			printTable(printV, headers);
			return res;
		}

		//reset values for next iteration
		lastErr = error;
		res = 0;
		N = N * 2;
		iter++;
	}
}



int main()
{
	
	Doub acc = pow(1, -3);
	



	system("pause");
	return 0;
}