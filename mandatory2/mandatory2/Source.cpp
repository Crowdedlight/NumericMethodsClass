#pragma once

#include <iostream>
#include "nr3.h"
#include "utilities.h"

using namespace std;
using namespace util;

template<class T>
Doub trapez(T &func, Doub a, Doub b, Doub acc);

struct func_test
{
	func_test() {}
	Doub operator()(Doub x) { return ; }
};

int main()
{
	Doub acc = pow(10, -3);

	func_test f1;

	cout << trapez(f1, 0, 4, acc) << endl;;


	system("pause");
	return 0;
}


template<class T>
Doub trapez(T &func, Doub a, Doub b, Doub acc)
{
	Doub N = 1;
	Doub h = 0;
	Doub sum = 0;
	Doub res = 0;
	Doub lastRes = 0;

	while(true)
	{
		//calculate h
		h = (b - a) / N;

		//sum
		for (int i = 1; i < N; i++)
			sum += func(a + i*h);

		//calculate rest
		res = h * (0.5*func(a) + 0.5*func(b) + sum);

		//finish check
		if (abs(lastRes - res) <= acc)
		{
			return res;
		}

		//reset values for next iteration
		lastRes = res;
		sum = 0;
		res = 0;
		N = N * 2;

	}
}