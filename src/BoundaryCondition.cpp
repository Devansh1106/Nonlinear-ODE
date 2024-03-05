#include <iostream>
#include <cassert>
#include "../include/Master.hpp"
BoundaryCondition::BoundaryCondition()
{
	IsNeumannLhs = false;
	IsNeumannRhs = false;
	IsDirichletLhs = false;
	IsDirichletRhs = false;
}

void BoundaryCondition::SetNeumannLhs(double Lhsbcvalue)
{
	assert(IsDirichletLhs == false);
	IsNeumannLhs = true;
	mLhsBcValue = Lhsbcvalue;
}

void BoundaryCondition::SetNeumannRhs(double Rhsbcvalue)
{
	assert(IsDirichletRhs == false);
	IsNeumannRhs = true;
	mRhsBcValue = Rhsbcvalue;
}

void BoundaryCondition::SetDirichletLhs(double Lhsbcvalue)
{
	assert(IsNeumannLhs == false);
	IsDirichletLhs = true;
	mLhsBcValue = Lhsbcvalue;
}

void BoundaryCondition::SetDirichletRhs(double Rhsbcvalue)
{
	assert(IsNeumannRhs == false);
	IsDirichletRhs = true;
	mRhsBcValue = Rhsbcvalue;
}