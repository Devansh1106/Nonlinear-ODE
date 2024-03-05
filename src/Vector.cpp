#include <iostream>
#include "../include/Master.hpp"

Vector::Vector(int anumNodes)
{
	mvNumNodes = anumNodes;
	mData = new double [mvNumNodes];
	for (int i = 0; i < mvNumNodes; i++)
	{
		mData[i] = 0.0;
	}
}

Vector::Vector(int anumNodes, char* t)
{
	muNumNodes = anumNodes;
	muData = new double [muNumNodes];
	for (int i = 0; i < muNumNodes; i++)
	{
		muData[i] = 0.0;
	}
}

Vector::Vector(const Vector& othervector)
{
	mData = new double [othervector.mvNumNodes];
	for (int i = 0; i < othervector.mvNumNodes; i++)
	{
		mData[i] = othervector.mData[i];
	}
}

Vector::Vector(const Vector& othervector, char* t)
{
	muData = new double [othervector.muNumNodes];
	for (int i = 0; i < othervector.muNumNodes; i++)
	{
		muData[i] = othervector.muData[i];
	}
}

Vector::~Vector()
{
	delete[] mData;
	delete[] muData;
}

double& Vector::operator[](int i)
{   
    return mData[i];
}

double& Vector::operator()(int i)
{   
    return muData[i];
}