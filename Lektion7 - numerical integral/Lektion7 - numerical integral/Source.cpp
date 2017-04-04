#pragma once

#include "nr3.h"
#include <iostream>
#include <utilities.h>

using namespace std;
using namespace util;

template<class T>
Doub midpoint(T &func, Doub a, Doub b, Doub acc);
template<class T>
Doub trapez(T &func, Doub a, Doub b, int N);
template<class T>
Doub simpson(T &func, Doub a, Doub b, int N);

#pragma region functions
struct func_one
{
	func_one() {}
	Doub operator()(Doub x) { return cos(pow(x, 2)) * exp(-x); }
};

struct func_two
{
	func_two() {}
	Doub operator()(Doub x) { return sqrt(x) * cos(pow(x, 2)) * exp(-x); }
};

struct func_three
{
	func_three() {}
	Doub operator()(Doub x) { return 1000 * exp(-1/x) * exp(-1/(1-x)); }
};

struct func_four
{
	func_four() {}
	Doub operator()(Doub x) { return 1/sqrt(x) * cos(pow(x, 2)) * exp(-x); }
};

#pragma endregion 

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

int main()
{
	Doub acc = pow(10, -8);

	func_one f1;
	func_two f2;
	func_three f3;
	func_four f4;

	Doub mid1 = midpoint(f1, 0, 1, acc);
	//Doub trap1 = trapez(f1, 0, 1, 100);
	//Doub simpson1 = simpson(f1, 0, 1, 100);


	system("pause");
	return 0;
}

template<class T>
Doub midpoint(T &func, Doub a, Doub b, Doub acc)
{
	//vector for print
	vector<Doub> indexV;
	vector<Doub> nV;
	vector<Doub> resultV;
	vector<Doub> alphaKV;
	vector<Doub> errorV;
	vector<string> headers = { "Iter", "N" , "A", "AlphaK", "Error" };


	bool running = true;
	Doub N = 1;
	Doub lastRes = 0;
	Doub res = 0;
	int iter = 0;
	Doub h = 0;
	//repeat calcs until integral stabilizes
	while(running)
	{		
		h = (b - a) / N;
		res = 0;
		for(int i = 1; i <= N; i++)
		{
			//sum of f(xi-1/2)
			Doub xi = (a + i*h) - (0.5*h);
			res += func(xi);
		}

		//final res
		res = res*h;

		indexV.push_back(iter);
		nV.push_back(N);
		resultV.push_back(res);
		
		//print saving
		if (iter > 1)
		{
			Doub alphaK = (resultV[iter - 2] - resultV[iter - 1]) / (resultV[iter - 1] - resultV[iter]);
			alphaKV.push_back(alphaK);
		}
		else
			alphaKV.push_back(0);

		//finish check
		if(abs(lastRes - res) < acc)
		{
			running = false;
			//calculate errors
			errorV.push_back(0); // first entry can't be calculated
			for(int j = 1; j < resultV.size(); j++)
			{
				Doub err = (resultV[j] - resultV[j - 1]) / (alphaKV.back() - 1);
				errorV.push_back(err);
			}

			//print stuff
			vector<vector<Doub>> printV;
			printV.push_back(indexV);
			printV.push_back(nV);
			printV.push_back(resultV);
			printV.push_back(alphaKV);
			printV.push_back(errorV);
			printTable(printV, headers);
			return res;
		}

		//reset values for next iteration
		lastRes = res;
		res = 0;
		N = N * 2;
		iter++;
	}
	return res*h;
}

template<class T>
Doub trapez(T &func, Doub a, Doub b, int N)
{
	Doub h = (b - a) / N;
	Doub sum = 0;
	Doub res = 0;

	//sum
	for (int i = 1; i < N; i++)
		sum += func(a+i*h);

	//calculate rest
	res = h * (0.5*func(a) + 0.5*func(b) + sum);

	return res;
}

template<class T>
Doub simpson(T &func, Doub a, Doub b, int N)
{
	Doub h = (b - a) / N;
	Doub res = 0;
	for (int i = 1; i < N; i++)
	{
		//sum of f(xi-1/2)
		Doub xi = (a + i*h) - (0.5*h);
		res += func(xi);
	}
	return res;
}