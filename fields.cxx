#include <iostream>
#include <cstdlib>
#include <cmath>
#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "utils.hpp"
#include "fields.hpp"
#include "solver.hpp"
using namespace std;
////////////////////////////////////////////////////////////////////
/////  N: row; M: Column; x: initial guess
///////////////////////////////////////////////////////////////////
//void grd( double **A, double *b, double *x, unsigned int N, unsigned int M);
//void CG(double **A, double *b, double *x, unsigned int N, unsigned int M);
//void Jacobi(double **A, double *b, double *x, unsigned int N, unsigned int M);
//void Gauss_seidel(double **A, double *b, double *x,unsigned int N,unsigned int M);

void allocate_variables(const Param &param, Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;

    var.volume = new double_vec(e);
    var.temperature = new double_vec(n);

    var.shpdx = new shapefn(e);
    var.shpdz = new shapefn(e);

    var.mat = new MatProps(param, var);
}



void update_temperature(const Param &param, const Variables &var,
                        double_vec &temperature, double_vec &tdot)
{

    // To be completed
	double *b, *x,tol=1e-5,sum1;
	unsigned int N=200, M=200,n,m;
	double **A = new double *[N];
	for(unsigned int i=0; i<N;i++){
		A[i] = new double[M];
	}
	b = new double[N];
	x = new double[N];

	double upper,lower, rk ;
	for(n=0;n<N;n++)
	{
		for(m=0;m<M;m++)
		{
			if(n==m)
			{
				A[n][m]=-2.0;
			}
			else if(n==m+1)
			{
				A[n][m]=1.0;
			}
			else if(n==m-1)
			{
				A[n][m]=1.0;
			}
			else
				A[n][m]=0.0;
//			cout << A[n][m] <<" ";
		}
//		cout << "\n";
	}
	for(n=0; n<N; n++)
	{
		b[n]=0.5;
		x[n]=0.5;
	}
////select the solver, if you do not use other solvers, please comment it.
//// each time could only use one solver.
//    grd( A, b, x, N, M);
	for(unsigned int step=0; step<1000000; step++)
	{	
		double sum=0.0;
		sum1 = CG(A, b, x, N, M);
//		sum1 = grd( A, b, x, N, M);
//		sum1 = Jacobi(A, b, x, N, M);
//		sum1 = Gauss_seidel(A, b, x, N, M);
		cout<<sqrt(sum1)<<"\n";
		if( sqrt(sum1) <= tol)
		{
//			cout<<"the input b is : "<<"\n";
//			for(unsigned int i=0; i<N; i++)
//			{
//				cout<<b[i]<<" ";
//			}
//			cout<<"\n";
//			cout<<"the estimated b value: "<<"\n";
//			for(unsigned int i=0; i<N; i++)
//			{
//				cout<< b1[i]<<" ";
//			}
//			cout<<"\n";
			std::cout <<"n step to converge(CG method):  " << step <<"\n" ;
			std::cout <<"the tolerance is : " << tol << "\n" ;
			break;
		}
	}

//    Jordan(A, b, x, N, M);
//    Gauss_seidel(A, b, x, N, M);

	for(unsigned int i=0; i<N; i++)
	{
    	delete [] A[i];
    }
    	delete [] A;
    	delete [] b;
    	delete [] x;
}
