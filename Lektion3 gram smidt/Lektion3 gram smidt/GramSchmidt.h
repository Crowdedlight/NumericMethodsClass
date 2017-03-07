#pragma once
#include "nr3.h"
#include "utilities.h"
#include <numeric>
#include <iostream>
#include <math.h>

using namespace std;
class GramSchmidt
{
public:
    GramSchmidt();    
    
    vector<VecDoub> calculate(vector<VecDoub> input);
    void print();
    ~GramSchmidt();

private:
    VecDoub sumCalc(VecDoub &input, vector<VecDoub> &eVec, int i);
    vector<VecDoub> eVectors;
};

