#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iterator>
#include <vector>
#include "constants.hpp"
#include "parameters.hpp"
#include "matprops.hpp"
#include "geometry.hpp"
#include "utils.hpp"
#include "fields.hpp"
#include "solver.hpp"
#include "ic.hpp"
using namespace std;

void allocate_variables(const Param &param, Variables& var)
{
    const int n = var.nnode;
    const int e = var.nelem;

    var.volume = new double_vec(e);
    var.temperature = new double_vec(n);
//    var.temperature = new double(n);
    var.shpdx = new shapefn(e);
    var.shpdz = new shapefn(e);

    var.mat = new MatProps(param, var);
}

//Variables  var;

//*var.temperature = initial_temperature(param,var,*var.temperature); 

 void update_temperature(const Param &param, const Variables &var,
                         double_vec &temperature, double_vec &tdot)
 {


 	cout << "element "<< var.nelem <<"\n";
 	cout << "this is the node"<< var.nnode <<"\n";
// //////////////////////////////////////////////////////////////////
// //
// //    Define variable and allocate memory
// //      	time_step: time step 
// ///		pho: density
// //		c: heat capacity
// //		k: heat diffusion coefficient
// //		f: heat production rate
// //////////////////////////////////////////////////////////////////
	double *F,*T,local_T,*f;
	double time_step=1.0,pho=1.29,c=1.004,k_diff_inside=1.9e-5,f1=2000,tol=1e-5,sum1;
	const int node = var.nnode;
	const int element = var.nelem;
/////////////////////////////////////////////////////////////////
//      Global stiffness matrix
///////////////////////////////////////////////////////////////
	double **K = new double *[var.nnode];
	for(int n=0; n<var.nnode;n++){
		K[n] = new double[var.nnode];
	}
	double **M = new double *[var.nnode];
	for(int n=0; n<var.nnode;n++){
		M[n] = new double[var.nnode];
	}
//	double **A = new double *[var.nnode];
//	for(int n=0; n<var.nnode;n++){
//		A[n] = new double[var.nnode];
//	}
/////////////////////////////////////////////////////////////////
///   heat production rate and temperature
////////////////////////////////////////////////////////////////
	F = new double[var.nnode];
	f = new double[3];
	T = new double[var.nnode];
	for(int i=0; i<var.nnode; i++)
	{
		T[i]=20.0;
//		cout<<T[i]<<"\n";
	}
	for(int i=0; i<var.nnode;i++)
	{
		cout <<"before"<< temperature[i] <<"\n";
	}
	for(int i=0; i<var.nnode; i++)
	{
		temperature.at(i)=T[i];
	}
	for(int i=0; i<var.nnode;i++)
	{
		cout <<"after"<< temperature[i] <<"\n";
	}
//	std::vector<double> vec( std::begin(T), std::end(T));
//	std::vector<double> temperature(std::begin(T),std::end(T));
	// for(int i=0; i<var.nnode; i++)
	// {
	// 	var.temperature[i]=T[i];
	// 	cout<<"the temperature is"<<var.temperature[i]<<"\n";
	// }

///////////////////////////////////////////////////////////////////
//     For local stiffness matrix
//////////////////////////////////////////////////////////////////
	double **k = new double *[3];
	for(unsigned int i=0; i<3;i++){
                 k[i] = new double[3];
        }
	double **m = new double *[3];
	for(unsigned int i=0; i<3;i++){
		m[i] = new double[3];
	}
	
	double upper,lower, rk ;

/////////////////////////////////////////////////////////////////////////////
//
//            Time dependent heat diffusion problem
/////////////////////////////////////////////////////////////////////////////
	
	double ke[3][3]= {
		{2, -1, -1},
		{-1, 1, 0},
		{-1, 0, 1}
	};

	for(int e=0; e<var.nelem; ++e)
	{
		const int  *conn = (*var.connectivity)[e];
		const array_t& coord = *var.coord;
		double f[3]={0.0,0.0,0.0};
		int n0 = conn[0];
		int n1 = conn[1];
		int n2 = conn[2];
		const double *a1 = coord[n0];
		const double *b1 = coord[n1];
		const double *c1 = coord[n2];
//		cout <<"the global node number "<< n0 <<"\n";
		for(int i=0; i<3; i++)
		{
			for(int j=0; j<3; j++)
			{
				m[i][j]=0.5*pho*c*phi(i)*phi(j)*2*((*var.volume)[e]);
//				cout<<m[i][j]<<"\n";
				M[conn[i]][conn[j]] += m[i][j];
			}
		}
		for(int i=0; i<3; i++)
		{
			if(((a1[0]+b1[0]+c1[0])/3)-100<=1 && ((a1[1]+b1[1]+c1[1])/3)-100<1)
			{			
				f[i]= 0.5*f1*phi(i)*2*((*var.volume)[e]);
			}
			else
			{
				f[i]=0;
			}
//			cout <<"the local f is "<<f[i]<<"\n";
			F[conn[i]] += f[i];
		}
// 		////////// If it is inside the domain 

		if((*var.bcflag)[n0] == 0 && (*var.bcflag)[n1]==0 && (*var.bcflag)[n2]==0)
		{		
			for(int i=0; i<3; i++)
			{
				for(int j=0; j<3; i++)
				{
//					k[i][j]=0.5*pho*c*phi(i)*phi(j)*2*((*var.volume)[e]) + time_step * 0.5*k_diff*ke[i][j]*2*((*var.volume)[e]);
					k[i][j]=0.5*k_diff_inside*ke[i][j]*2*((*var.volume)[e]);
//					b[i]=b[i] + (0.5)*pho*c*phi(i)*phi(j)*2*((*var.volume)[e]) * local[j];
					K[conn[i]][conn[j]] += k[i][j];
//					A[i][j]=M[i][j]+time_step*K[i][j];
				}
				
			}

		}
// 		////////////////if two nodes are on the boundary( 0,1 corresponded global nodes on boundary)
		else if((*var.bcflag)[n0] != 0 && (*var.bcflag)[n1] != 0 && (*var.bcflag)[n2] == 0)
 		{
			
			k[0][0]=1;
			k[1][1]=1;
			k[2][2]= 0.5*k_diff_inside*ke[2][2]*2*((*var.volume)[e]);
			k[0][1]=0;k[0][2]=0;
			k[1][0]=0;k[1][2]=0;
			k[2][0]=0;k[2][1]=0;
			
			for(int i=0; i<3; i++)
			{	
				for(int j=0; j<3; j++)
				{
					K[conn[i]][conn[j]] += k[i][j];
				}
 			}
			
 		}
// 		////////////////if two nodes are on the boundary( 0, 2 nodes on boundary)
		else if((*var.bcflag)[n0] !=0 &&(*var.bcflag)[n1]==0 && (*var.bcflag)[n2] != 0 )
		{
			k[0][0]=1;
			k[2][2]=1;
			k[1][1]= 0.5*k_diff_inside*ke[1][1]*2*((*var.volume)[e]);
			k[0][1]=0;k[0][2]=0;
			k[1][0]=0;k[1][2]=0;
			k[2][0]=0;k[2][1]=0;
			for(int i=0; i<3; i++)
			{
				for(int j=0; j<3; i++)
				{
					K[conn[i]][conn[j]] += k[i][j];
				}
			}
		}
// 		////////////////if two nodes are on the boundary( 1, 2 nodes on boundary)
		else if((*var.bcflag)[n0] == 0 && (*var.bcflag)[n1] != 0 && (*var.bcflag)[n2] != 0)
		{
			k[2][2]=1;
			k[1][1]=1;
			k[0][0]= 0.5*k_diff_inside*ke[0][0]*2*((*var.volume)[e]);
			k[0][1]=0;k[0][2]=0;
			k[1][0]=0;k[1][2]=0;
			k[2][0]=0;k[2][1]=0;
			for(int i=0; i<3; i++)
			{
				for(int j=0; j<3; j++)
				{
					K[conn[i]][conn[j]] += k[i][j];
				}
			}
		}
// 		/////  node 0 on the boundary , but nodes 1, 2 are inside the domain
		else if ((*var.bcflag)[n0] != 0 && (*var.bcflag)[n1] == 0 && (*var.bcflag)[n2] == 0)
		{
			cout<<"node 0 corresoned global node on the boundary";
			for(int i=1; i<3; i++)
			{
				for(int j=1; j<3; j++)
				{
					k[i][j]=0.5*k_diff_inside*ke[i][j]*2*((*var.volume)[e]);
				}
			}
			k[0][0]=1; k[0][1]=0;k[0][2]=0;
			k[1][0]=0;k[2][0]=0;

			for(int i=0; i<3; i++)
			{
				for(int j=0; j<3; i++)
				{
					K[conn[i]][conn[j]] += k[i][j];
				}
			}

		}
// 		/////  node 1 on the boundary , but nodes 0, 2 are inside the domain
		else if((*var.bcflag)[n0]== 0 &&(*var.bcflag)[n1] != 0 && (*var.bcflag)[n2] == 0)
		{	
			k[0][0]=0.5*k_diff_inside*ke[0][0]*2*((*var.volume)[e]);
			k[0][1]=0; k[0][2]=0;
			k[1][0]=0; k[1][1]=1; k[1][2]=0;
			k[2][0]=0.5*k_diff_inside*ke[2][0]*2*((*var.volume)[e]);
			k[2][1]=0;
			k[2][2]=0.5*k_diff_inside*ke[2][2]*2*((*var.volume)[e]);
			for(int i=0; i<3; i++)
			{
				for(int j=0; j<3; i++)
				{
					K[conn[i]][conn[j]] += k[i][j];
				}
			}
		}
// 		/////  node 2 on the boundary , but nodes 0, 1 are inside the domain
		else if ((*var.bcflag)[n0]== 0 &&(*var.bcflag)[n1] == 0 && (*var.bcflag)[n2] != 0)
		{
			for(int i=0; i<2; i++)
			{
				for(int j=0; j<2; j++)
				{	
					k[i][j]=0.5*k_diff_inside*ke[i][j]*2*((*var.volume)[e]);
				}
			}	
			k[0][2]=0;k[1][2]=0;k[2][0]=0;k[2][1]=0;k[2][2]=1;
			for(int i=0; i<3; i++)
			{
				for(int j=0; j<3; i++)
				{
					K[conn[i]][conn[j]] += k[i][j];
				}
			}
		}
		else if ((*var.bcflag)[n0] != 0 &&(*var.bcflag)[n1] != 0 && (*var.bcflag)[n2] != 0)
		{
			k[0][0]=1;k[1][1]=1;k[2][2]=1;
			k[0][1]=0;k[0][2]=0;k[1][0]=0;k[1][2]=0;k[2][0]=0;k[2][1]=0;
			for(int i=0; i<3; i++)
			{
				for(int j=0; j<3; j++)
				{
					K[conn[i]][conn[j]] += k[i][j];
				}
			}
		}

 	}
// /////////////////////////////////////////////////////////////////////////////
// ////select the solver, if you do not use other solvers, please comment it.
// //// each time could only use one solver.
// // 
// ////////////////////////////////////////////////////////////////////////////
	for(unsigned int step=0; step<1000000; step++)
	{	
		double sum=0.0;
//		sum1 = CG(A, F, T, node, node);
//		sum1 = grd( A, b, x, N, M);
//		sum1 = Jacobi(A, b, x, N, M);
//		sum1 = Gauss_seidel(A, b, x, N, M);
//		cout<<sqrt(sum1)<<"\n";
//		if( sqrt(sum1) <= tol)
//		{
//			std::cout <<"n step to converge(CG method):  " << step <<"\n" ;
//			std::cout <<"the tolerance is : " << tol << "\n" ;
//			break;
//		}
	}


	for(int i=0; i<var.nnode; i++)
	{
    	delete [] K[i];
		delete [] M[i];
    }
    delete [] K;
	delete [] M;

	for(unsigned int i=0; i<3;i++)
	{
		delete [] k[i];
		delete [] m[i];
//		delete [] ke[i];
	}
	delete [] k;
	delete [] m;
//	delete [] ke;
    delete [] F;
	delete [] T;
	delete [] f;
//	delete [] local;
}
