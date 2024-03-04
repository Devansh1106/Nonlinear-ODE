#include<cmath>
#include<iostream>
#include "Master.hpp"
using namespace std;
double model_prob_rhs(double x){return 0;}

int main (int argc, char** argv)
{
	BoundaryCondition bc_mp;
	Input mp(1.0,0.0,1.0,&model_prob_rhs);
	Grid grid_1(0.0,1.0,500);
	int max_iter = 10;
	bc_mp.SetDirichletLhs(0.0);	//original = -5,4
    bc_mp.SetDirichletRhs(1.0);
    BvpOde bvpode(&mp, &bc_mp, &grid_1);
	Solving(argc, argv, &mp, &bc_mp, &grid_1, &bvpode, max_iter);
	return 0;
}