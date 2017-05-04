#include "nr3.h"
#include "utilities.h"
#include <iostream>
#include "stepper.h"
#include "odeint.h"
#include "stepperdopr853.h"
#include "shoot.h"
#include "shootf.h"
#include "ludcmp.h"
#include "roots.h"
#include "qrdcmp.h"
#include "roots_multidim.h"

using namespace std;
using namespace util;

struct rhs {
	void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
		dydx[0] = y[1];
		dydx[1] = -cos(y[0]) * sin(y[1]);
	}
};

struct load
{
	VecDoub operator() (const Doub x, VecDoub_I &y)
	{
		VecDoub ystart(2);
		ystart[0] = 0;
		ystart[1] = y[0];
		return ystart;
	}
};

struct score
{
	Doub t;
	score(Doub target) :t(target) {}
	VecDoub operator() (const Doub x, VecDoub_I &y)
	{
		VecDoub error(1);
		error[0] = y[0] - t;
		return error;
	}
};

int main()
{
	Int nvar = 2;
	Doub x1 = 0.0, x2 = 10.0;
	load l;
	rhs d;
	score s(3);
	VecDoub v(1);
	v[0] = 0; //Need to init otherwise it gets init wrong
	bool check = false;

	Shoot<load, rhs, score> shoot(nvar, x1, x2, l, d, s);
	newt(v, check, shoot);

	v.print();

	system("pause");
	return 0;
}