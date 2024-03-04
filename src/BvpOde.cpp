#include <iostream>
#include <cassert>
#include <fstream>
#include <cstdio>
#include "Master.hpp"


BvpOde::BvpOde(Input* Ode, BoundaryCondition* tbc_mp, Grid* grid_1)
{
	mOde = Ode;
	mgrid_1 = grid_1;
	n1 = mgrid_1->Nodes;
	mbc_mp = tbc_mp;
	assert(mbc_mp != NULL);
	char* r = NULL;
	mu = new Vector(n1, r);
	mb = new Vector(n1);
	mA = new Matrix(n1);
	assert(mA != NULL);
}

BvpOde::~BvpOde()
{
	delete mA;
	delete mb;
	delete mu;
}

void BvpOde::ApplyBcVector()
{
	bool right_bc_applied = false;
	bool left_bc_applied = false;

	if (mbc_mp->IsNeumannLhs == true)
	{
		(*mb)[0] = mbc_mp->mLhsBcValue;
		left_bc_applied = true;
	}
																	//Working fine
	if (mbc_mp->IsNeumannRhs == true)
	{
		(*mb)[n1-1] = mbc_mp->mRhsBcValue;
		right_bc_applied = true;
	}

	if (mbc_mp->IsDirichletLhs == true)
	{
		(*mb)[0] = mbc_mp->mLhsBcValue;
		left_bc_applied = true;
	}

	if (mbc_mp->IsDirichletRhs == true)
	{
		(*mb)[n1-1] = mbc_mp->mRhsBcValue;
		right_bc_applied = true;
	}
	assert(right_bc_applied && left_bc_applied);
}

void BvpOde::ApplyBcMatrix()
{
	assert(mA != NULL);
	bool right_bc_applied_1 = false;
	bool left_bc_applied_1 = false;
	if (mbc_mp->IsDirichletLhs == true)
	{
		(*mA)(0,0) = 1;
		left_bc_applied_1 = true;
	}

	if (mbc_mp->IsNeumannRhs == true)
	{
		double h = mgrid_1->Mesh[n1] - mgrid_1->Mesh[n1-1];
		(*mA)(n1-1, n1-2) = -1/h;
		(*mA)(n1-1, n1-1) = 1/h;
		right_bc_applied_1 = true;
	}

	if (mbc_mp->IsNeumannLhs == true)
	{
		double h = mgrid_1->Mesh[1] - mgrid_1->Mesh[0];
		(*mA)(0,0) = -1/h;
		(*mA)(0,1) = 1/h;
		left_bc_applied_1 = true;
	}

	if (mbc_mp->IsDirichletRhs == true)
	{
		(*mA)(n1-1,n1-1) = 1;
		right_bc_applied_1 =  true;
	}


	assert(left_bc_applied_1);
	assert(right_bc_applied_1);
}

Vector* BvpOde::InitialGuess()
{
	assert(mu != NULL);
	double num;
	int i = 0;
	FILE *file_1 = fopen("solution_values.txt", "r");
	if (file_1 != NULL)
	{																		//Working fine
		while(fscanf(file_1, "%lf", &num)!=EOF and i < n1)
	  	{
	  		(*mu)(i) = num;
	  		i++;
	  	}
	}
	else {std::cout<<"Error opening file InitialGuess.txt! ";}
	fclose(file_1);
	return mu;
}

Matrix* BvpOde::MatrixGen(Vector& mpu, int iter)
{
	mmu = &mpu;
	double A = mOde->mCoeffUxx;
	double B = mOde->mCoeffUx;
	double C = mOde->mCoeffU;
	for (int i = 1; i < n1-1; i++)
	{
		double xm = mgrid_1->Mesh[i-1];
		double x = mgrid_1->Mesh[i];
		double xp = mgrid_1->Mesh[i+1];
		double u = (*mmu)(i);
		double alpha = 2.0*A/(xp-xm)/(xp-x) + B/(xp-xm);
		double beta = -2.0*A/(xp-xm)/(xp-x) - 2.0*A/(xp-xm)/(x-xm) + 2*C*u;
		double gamma = 2.0*A/(xp-xm)/(x-xm) - B/(xp-xm);
		(*mA)(i, i-1) = gamma;
		(*mA)(i, i) = beta;
		(*mA)(i,i+1) = alpha;
	}
	if(iter == 4){
		std::ofstream f2("matrix.txt");
	    for(int i = 0; i < n1; i++)
		{
			for(int j = 0; j < n1; j++)
			{
				f2<<(*mA)(i,j)<<" ";
			}
			f2<<"\n";
		}
	}
	return mA;
}

Vector* BvpOde::VectorGen(Vector& mpu, int iter)
{
	mmu = &mpu; 
	double A = mOde->mCoeffUxx;
	double B = mOde->mCoeffUx;
	double C = mOde->mCoeffU;
	for (int i = 1; i < n1-1; i++)
	{
		double xm = mgrid_1->Mesh[i-1];
		double x = mgrid_1->Mesh[i];
		double xp = mgrid_1->Mesh[i+1];
		double up = (*mmu)(i+1);
		double u = (*mmu)(i);
		double um = (*mmu)(i-1);
		double rhs = mOde->mRhsFunc(x);
		double alpha = 2.0*A/(xp-xm)/(xp-x) + B/(xp-xm);
		double beta = -2.0*A/(xp-xm)/(xp-x) - 2.0*A/(xp-xm)/(x-xm) + C*u;
		double gamma = 2.0*A/(xp-xm)/(x-xm) - B/(xp-xm);
		(*mb)[i] = (-1.0) * (alpha*up + beta*u + gamma*um - rhs);
	}
	// if (iter == 4){
	std::ofstream f21("rhs_vec.txt");
    for(int i = 0; i < n1; i++)
	{
		f21<<(*mb)[i]<<"\n";
	}
	f21.close();
	// }
	return mb;
}