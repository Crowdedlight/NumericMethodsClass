#include "GramSchmidt.h"


//Tester
GramSchmidt::GramSchmidt()
{
}


GramSchmidt::~GramSchmidt()
{
}

vector<VecDoub> GramSchmidt::calculate(vector<VecDoub> input)
{    
    vector<VecDoub> eVec(input.size);
    for (auto vec : input)
    {
        //e1 := x1/||x1||
        eVec[0] = (vec / vec.norm);

        for (int i = 1; i < input.size(); i++)
        {
            //ei := xi - sum
            eVec[i] = input[i] - sumCalc(input[i], eVec, i);
            //ei := ei / ||ei||
            eVec[i] = eVec[i] / eVec[i].norm;
        }
    }
    this->eVectors = eVec;
    return eVec;
}

void GramSchmidt::print()
{
    //use utilities to print vectors
}

VecDoub GramSchmidt::sumCalc(VecDoub &input, vector<VecDoub> &eVec, int i)
{
    VecDoub sum;
    for (int j = 0; j < i - 1; j++)
    {
        sum = (input * eVec[j]) * eVec[j];
    }
    return sum;
}


