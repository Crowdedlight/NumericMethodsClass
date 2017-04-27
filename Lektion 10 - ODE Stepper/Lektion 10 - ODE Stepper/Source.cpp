#pragma once

#include "nr3.h"
#include <iostream>
#include "utilities.h"
#include "rk4.h"
#include "stepperdopr5.h"
#include "stepperross.h"
#include "ludcmp.h"
#include "stepper.h"
#include "odeint.h"

using namespace std;
using namespace util;

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
				cout << setw(10) << right << (colums[j][i]) << setw(3) << "|";
			else
				cout << setw(15) << colums[j][i] << setw(3) << "|";
		}
		cout << endl;
	}

	//End spacing 
	cout << endl << endl << endl;
}
#pragma endregion 


struct rhs {
	rhs(){}
	void operator() (const Doub x, VecDoub &y, VecDoub_O &dydx) {
		dydx[0] = y[0]*y[1];
		dydx[1] = -pow(y[0], 2);
	}
	void jacobian(const Doub x, VecDoub_I &y, VecDoub_O &dfdx, MatDoub_O &dfdy) {
		Int n = y.size();
		for (Int i = 0; i<n; i++) dfdx[i] = 0.0;
		dfdy[0][0] = y[1];
		dfdy[0][1] = y[0];
		dfdy[1][0] = -2*y[0];
		dfdy[1][1] = 0;
	}
};

int main()
{
	//**********************************//
	//********** STEPPERDOPR ***********//
	//**********************************//

	const Int nvar = 2;
	const Doub atol = pow(10, -6), rtol = atol, h1 = 0.01, hmin = 0.0, x1 = 0.0, x2 = 10.0;
	VecDoub ystart(nvar);
	ystart[0] = 1.0;
	ystart[1] = 1.0;
	Output out(20); //Dense output at 40 points plus x1.
		rhs d; //Declare d as a rhs object.
		Odeint<StepperDopr5<rhs> > ode(ystart, x1, x2, atol, rtol, h1, hmin, out, d);
	ode.integrate();

	//**********************************//
	//********** STEPPERROSS ***********//
	//**********************************//
	ystart[0] = 1.0;
	ystart[1] = 1.0;
	Output out2(20); //Dense output at 20 points plus x1.
	Odeint<StepperRoss<rhs> > ode2(ystart, x1, x2, atol, rtol, h1, hmin, out2, d);
	ode2.integrate();


	//put in vector for table
	vector<double> vStep;
	vector<double> vX;
	vector<double> vY;
	vector<double> vConstant;
	vector<double> vConstantRoss;
	vector<double> vStepRoss;
	vector<double> vRossX;
	vector<double> vRossY;
	vector<string> vHeader = { "N", "Xres", "Yres", "Constant", "N_ross", "Xres Ross", "Yres Ross", "Constant Ross" };

	for (Int i = 0; i<out.count; i++)
	{
		vStep.push_back(out.xsave[i]);
		vX.push_back(out.ysave[0][i]);
		vY.push_back(out.ysave[1][i]);
		vConstant.push_back(pow(out.ysave[0][i], 2) + pow(out.ysave[1][i], 2));

		vStepRoss.push_back(out2.xsave[i]);
		vRossX.push_back(out2.ysave[0][i]);
		vRossY.push_back(out2.ysave[1][i]);
		vConstantRoss.push_back(pow(out2.ysave[0][i], 2) + pow(out2.ysave[1][i], 2));
	}

	vector<vector<double>> printV = { vStep, vX, vY, vConstant, vStepRoss, vRossX, vRossY, vConstantRoss };
	printTable(printV, vHeader);

	system("pause");
	return 0;
}
