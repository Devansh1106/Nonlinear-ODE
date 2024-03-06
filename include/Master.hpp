#ifndef MASTERHEADERDEF
#define MASTERHEADERDEF
#include <iostream>
#include <vector>

class Input {
public:
	friend class BvpOde;
	Input(double coeffUxx,double coeffUx,double coeffU,double (*righthandside)(double));
	double rhsfunc(double x);
	double mCoeffUxx;
	double mCoeffUx;
	double mCoeffU;
	double (*mRhsFunc) (double x);
};

class Matrix {
private:
	double** mData;
	int mpNumNodes;
public:
	friend class BvpOde;
	Matrix(const Matrix& othermatrix);
	Matrix(int vnumNodes);
	~Matrix();
	double& operator() (int i, int j);
};

class Vector {
private:
	int mvNumNodes;
	int muNumNodes;
	double* mData;
	double* muData;
public:
	friend class BvpOde;
	Vector(int anumNodes);
	Vector(int anumNodes, char* t);
	Vector(const Vector& othervector);
	Vector(const Vector& othervector, char* t);
	void Vector_delete();
	~Vector();
	double& operator[](int i);
	double& operator()(int i);
};

class BoundaryCondition {
public:
	bool IsNeumannLhs;
	bool IsNeumannRhs;
	bool IsDirichletLhs;
	bool IsDirichletRhs;
	double mLhsBcValue;
	double mRhsBcValue;
	friend class BvpOde;
	BoundaryCondition();
	void SetNeumannLhs(double Lhsbcvalue);
	void SetNeumannRhs(double Rhsbcvalue);
	void SetDirichletRhs(double Rhsbcvalue);
	void SetDirichletLhs(double Lhsbcvalue);
};

class Grid {
public:
	std::vector<double> Mesh;
	Grid(double xMinimum, double xMaximum, int mNodes);
	double xMin , xMax;
	int Nodes;
};

class BvpOde {
private:
	Matrix* mA;
	Vector* mb;
	Vector* mu;
	Vector* mmu;
	int n1;
	Input* mOde;
	BoundaryCondition* mbc_mp;
	Grid* mgrid_1;
	
public:

	friend class LinearSystem;
	BvpOde(Input* Ode, BoundaryCondition* tbc_mp, Grid* grid_1);
	~BvpOde();
	void ApplyBcMatrix();
	void ApplyBcVector();
	Matrix* MatrixGen(Vector& mpu);
	Vector* VectorGen(Vector& mpu);
	Vector* InitialGuess();
};

#ifndef __cplusplus
extern "C"
{
#endif
	void Solving(int argc, char** argv, Input* Ode, BoundaryCondition* tbc_mp, Grid* grid_1, BvpOde* bvpode, int max_iter);
#ifndef __cplusplus
}
#endif
#endif