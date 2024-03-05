#include <iostream>
#include "../include/Master.hpp"
Input::Input(double coeffUxx,double coeffUx,double coeffU,double (*righthandside)(double))
{
	mCoeffUxx = coeffUxx;
	mCoeffUx = coeffUx;
	mCoeffU = coeffU;
	mRhsFunc = righthandside;
}
