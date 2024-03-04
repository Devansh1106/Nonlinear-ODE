#include <math.h>
#include <fstream>
#include <complex>
#include <assert.h>
#include <cstdlib>
#include "Master.hpp"
#include "mpi.h"
#include "dmumps_c.h"

extern "C++"
{
	Vector* U;
	Matrix* M_final;
	Vector* V_final;
}

MUMPS_INT n;
MUMPS_INT8 len_a = 0;	//len_a = len_irn = len_jcn
DMUMPS_STRUC_C id;

void Solving(int argc, char** argv, Input* Ode, BoundaryCondition* tbc_mp, Grid* grid_1, BvpOde* bvpode, int max_iter)
{
	int N = grid_1->Nodes;
	n = N;
	int* irn;
	int* jcn;
	double* rhs;
	double* a;
	double time_taken,time_taken_2;
	double* sum_rhs;
	double* sum_sol;
	double saved_err_rhs;
	double err_rhs, err_sol;
	static int flag_2 = 0;
	static int flag_3 = 0;
	static int flag_1 = 0;
	int myid, ierr;
	int error =0;

	ierr = MPI_Init(&argc, &argv);		//MPI initialization
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	if(myid == 0){time_taken = MPI_Wtime();}		//record the start time

	// if (myid == 0){
	// 	std::ofstream a1("a1.txt");
	// 	// for(size_t f = 0; f < N; f++){
	// 	// 	a1<<rhs[f]<<"\n";
	// 	// }
	// 	a1.close();
	// 	// std::cout<<"rank 7";
	// }

	sum_rhs = (double*)malloc(1 * sizeof(double));
	// *sum_rhs = 0;

	if (myid == 0)
	{
		sum_sol = (double*)malloc(1 * sizeof(double));
		*sum_sol = 0;		
	}

	if (tbc_mp->IsNeumannLhs)
	{
		len_a = (N-2)*3 + 2;
	}

	if (tbc_mp->IsDirichletLhs)
	{
		len_a = (N-2)*3 + 1;
	}

	if (tbc_mp->IsNeumannRhs)
	{
		len_a = len_a + 2;
	}

	if (tbc_mp->IsDirichletRhs)
	{
		len_a = len_a + 1;
	}

	/*****Memory Allocation*****/
	irn = (int*)malloc(len_a * sizeof(int));
	jcn = (int*)malloc(len_a * sizeof(int));
	rhs = (double*)malloc(N * sizeof(double));
	a = (double*)malloc(len_a * sizeof(double));

	// if (myid == 1){
	// 	std::ofstream a1("a1.txt");
	// 	// for(size_t f = 0; f < N; f++){
	// 	a[2] = 0;
	// 	a1<<a[2]<<"\n";
	// 	// }
	// 	a1.close();
	// 	// std::cout<<"rank 7";
	// }
	
	int iter = 0;
	while (iter < max_iter)
	{
		*sum_rhs = 0;
		err_sol = 0;
		if (myid == 0)
		{
			if (flag_3 == 0)
			{
				U = bvpode->InitialGuess();
				// std::ofstream f0("che.txt");
				bvpode->ApplyBcMatrix();
				bvpode->ApplyBcVector();
				flag_3 = 1;
			}
			// std::ofstream f121("update.txt");
    		// for(int i = 0; i < N; i++)
			// {
			// 	f121<<(*U)(i)<<"\n";
			// }
			// f121.close();

			V_final = bvpode->VectorGen(*U, iter);
			// double summation = 0;
			// for(int i = 0; i < N; i++){
			// 	double y = (*V_final)[i];
			// 	y = std::pow(y,2.0);
			// 	summation = summation + y;
			// }
			// summation = std::sqrt(summation);
			// if(summation < 1e-8){
			// 	break;
			// }
			// printf("\nSummation - %lf\n",summation);
			M_final = bvpode->MatrixGen(*U, iter);
			// if (iter == 4){
			// 	std::ofstream mat("mat4.txt");
			// 	std::ofstream vec("vec4.txt");
			// 	for (size_t i = 0; i < N; i++){
			// 		for (size_t j = 0; j < N; j++)
			// 		{
			// 			mat<<(*M_final)(i,j)<<" ";
			// 		}
			// 		mat<<"\n";
			// 		vec<<(*V_final)[i]<<"\n";
			// 	}
			// mat.close();
			// vec.close();
			// }

			int k = 0;
			for (int i = 0; i < N; i++)
			{
				rhs[i] = (*V_final)[i];
				for (int j = 0; j < N; j++)
				{
					if((*M_final)(i,j) != 0)
					{
						if(flag_2 == 0)
						{
							irn[k] = i+1;
							jcn[k] = j+1;							
						}
						a[k] = (*M_final)(i,j);
						k++;
					}
				}
			}
			if (flag_2 == 0)
			{
				MPI_Bcast(irn, len_a, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(jcn, len_a, MPI_INT, 0, MPI_COMM_WORLD);
				flag_2 = 1;			
			}
			MPI_Bcast(a, len_a, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(rhs, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}
		else
		{
			if(flag_2 == 0)
			{
				MPI_Bcast(irn, len_a, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(jcn, len_a, MPI_INT, 0, MPI_COMM_WORLD);				
				flag_2 = 1;
			}
			MPI_Bcast(a, len_a, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(rhs, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}

		// if (myid == 3){
		// 	std::ofstream a1("a1.txt");
		// 	for(size_t f = 0; f < len_a; f++){
		// 		a1<<a[f]<<"\n";
		// 	}
		// 	a1.close();
		// 	// std::cout<<"rank 7";
		// }

		MPI_Barrier(MPI_COMM_WORLD);
		// if (myid == 7){
		// 	// std::ofstream a1("a1.txt");
		// 	// for(size_t f = 0; f < N; f++){
		// 	// 	a1<<rhs[f]<<"\n";
		// 	// }
		// 	// a1.close();
		// 	std::cout<<"rank 7";
		// }
		id.par = 1; id.sym = 0;
		id.job = -1;
		dmumps_c(&id);
		if(myid == 0)
		{	
			id.n = n; id.nnz = len_a; id.irn = irn; id.jcn = jcn;
			id.a = a; id.rhs = rhs;
		}

		// if (iter == 4){
		// 	std::ofstream mat("mat4.txt");
		// 	std::ofstream vec("vec4.txt");
		// 	for (size_t i = 0; i < len_a; i++){
		// 		// for (size_t j = 0; j < N; j++)
		// 		// {
		// 		// 	mat<<(*M_final)(i,j)<<" ";
		// 		// }
		// 		// mat<<"\n";
		// 		// vec<<(*V_final)[i]<<"\n";
		// 		vec<<id.jcn[i]<<"\n";
		// 	}
		// 	mat.close();
		// 	vec.close();
		// }
		// if (myid == 3){
		// 	std::cout<<"rank 3";
		// }

		#define ICNTL(I) icntl[(I)-1] 
		id.ICNTL(1) = -1; id.ICNTL(2) = -1; id.ICNTL(3) = -1; id.ICNTL(4) = 0; //Supressing error msgs

		// Call the MUMPS package (analyze , factorization and solve)
		id.job = 6;
		dmumps_c(&id);

		if (id.infog[0] < 0)
		{
			printf("(PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",myid, id.infog[0], id.infog[1]);
			error = 1;
		}
		//Terminate instance
		id.job = -2;
		dmumps_c(&id);
		// char buffer[50];
		// sprintf(buffer, "mumps_sol%d", iter);
		// std::ofstream f001(buffer);
		// for (size_t c = 0; c < N; c++){
		// 	f001<<(*U)(c)<<"\t"<<rhs[c]<<"\n";
		// }
		// f001.close();
		if (myid == 0)
		{
			std::ofstream f09("anjs.txt");
			double err = 0;
			double old_rhs[N];
			
			for (int i=0; i < N; i++)
			{
				// if(i == 0)
				// 	err = abs(rhs[i]);
				// else
				// 	err = abs(rhs[i]);
				err = pow(rhs[i],2);
				// err = pow((*V_final)[i],2);
				*sum_sol += err;
				(*U)(i) += rhs[i]; 
				f09<<(*V_final)[i]<<" ";
				f09<<*sum_sol<<"\n";
			}
			for(int z = 0; z < N; z++){
				old_rhs[z] = rhs[z];
			}

			f09.close();
			std::ofstream f1("sum.txt");
			*sum_sol = *sum_sol/(double)(N);
			*sum_rhs = sqrt(*sum_sol);

			f1<<"Sum "<<(*sum_rhs);
			f1.close();
		}
		MPI_Bcast(sum_rhs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if((*sum_rhs) < 1e-8)
		{
			break;
		}
		else
		{
			iter++;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	if (iter < max_iter)
	{
		if(myid == 0)
		{
			FILE *file_sol = fopen("Solution.txt","w");
	    	if(file_sol != NULL)
	    	{
	    		fprintf(file_sol,"%d\n" ,iter);
	    		for(int i=0; i < N; i++)
	    		{
	    			double M = grid_1->Mesh[i];
	    			fprintf(file_sol, "%lf %lf\n",M, (*U)(i));
	    		}
	    		printf("Solution is in Solution.txt file!");
	    	}
	    	else
	    	{
	    		printf("Error opening solution file!");
	    	}
	    	fclose(file_sol);
	    }
	}
	else
	{
		if(myid == 0)
		{
			printf("Newton-Raphson method fails! Exceeded max iterations. ");
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	free(irn);
	free(jcn);
	free(a);
	free(rhs);
	free(sum_rhs);
	if (myid == 0){free(sum_sol);
	time_taken_2 = MPI_Wtime();
	printf("\n%f",time_taken_2 - time_taken);}
	ierr = MPI_Finalize();
}
