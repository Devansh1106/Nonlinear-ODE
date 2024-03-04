#include <iostream>
#include <cassert>
#include <fstream>
#include "Master.hpp"

Matrix::Matrix(const Matrix& othermatrix)
{
	mData = new double* [othermatrix.mpNumNodes];
	for(int i = 0; i < othermatrix.mpNumNodes; i++)
	{
		mData[i] = new double [othermatrix.mpNumNodes];
	}
	for(int i = 0; i < othermatrix.mpNumNodes; i++)
	{
		for(int j = 0; j < othermatrix.mpNumNodes; j++)
		{
			mData[i][j] = othermatrix.mData[i][j];
		}
	}
}

Matrix::Matrix(int vnumNodes)
{
	assert(vnumNodes > 0);
	// std::ofstream f2("is11.txt");
	mpNumNodes = vnumNodes;
	mData = new double* [mpNumNodes];
	for(int i = 0; i < mpNumNodes; i++)
	{
		mData[i] = new double [mpNumNodes];
	}
	for(int i = 0; i < mpNumNodes; i++)
	{
		for(int j = 0; j < mpNumNodes; j++)
		{
			mData[i][j] = 0.0;
		}
	}
	// for(int i = 0; i < mpNumNodes; i++)
	// {
	// 	for(int j = 0; j < mpNumNodes; j++)
	// 	{
	// 		f2<<mData[i][j]<<" ";
	// 	}
	// 	f2<<"\n";
	// }
	// f2.close();
}

Matrix::~Matrix()
{
    for (int i = 0; i < mpNumNodes; i++)
    {
        delete[] mData[i];
    }
    delete[] mData;
}

double& Matrix::operator()(int i, int j)
{   
    return mData[i][j];
}