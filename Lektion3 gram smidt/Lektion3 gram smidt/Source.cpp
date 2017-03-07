#include <iostream>
#include "nr3.h"
#include "GramSchmidt.h"
using namespace std;

int main()
{	
    //init class
    GramSchmidt myGS();
    
    //crate input arrays
    Doub x1_a[] = { 3, 0, 4, 12 };
    Doub x2_a[] = { 1, 1, 0, 0 };
    Doub x3_a[] = { -1, 1, 2, 0 };
    VecDoub x1(4, x1_a), x2(4, x2_a), x3(4, x3_a);

    vector<VecDoub> input;
    input.push_back(x1);
    input.push_back(x2);
    input.push_back(x3);

    

	system("pause");
	return 0;
}