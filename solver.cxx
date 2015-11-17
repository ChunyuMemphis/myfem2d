#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace std;
void grd(double **A, double *b, double *x, int N, int M)
{
	double dot[N];
	double r[N];
	double upper=0.0, lower=0.0, rk;
//	for(unsigned int i = 0; i < N; i++)
//	{
//		x[i] = 1.0;
//	}
	for(unsigned int step=0; step<500; step++)
	{
		
//		double r_tran[1][2]={0.0,0.0};
//		double r_tran_A[1][2]={0.0,0.0};
        double upper=0.0;
        double lower=0.0;
        double sum=0.0;
		// zero the dot array
		for(unsigned int i=0; i < N; i++)
		{
			dot[i]=0.0;
		}
	    // compute A * x0
		for(unsigned int i = 0; i < N; i++)
	    {
	    	for(unsigned int j = 0; j < M; j++)
			{
				dot[i] += A[i][j] * x[j];
//			    std::cout << dot[i] <<" ";
		    }
//		    std::cout << "\n";
		}

		// compute residual
	    for(unsigned int i = 0; i < N; i++)
	    {
	    	r[i] = b[i] - dot[i];
//		    std::cerr << r[i] << "\n";
	    }

	    // r^T*r
        for(unsigned int i = 0; i < N; i++)
	    {
			upper += r[i] * r[i];
	    }

	    // r^T*A*r
		for(unsigned int i = 0; i < N; i++)
	    {
	    	for(unsigned int j = 0; j < M; j++)
			{
				lower +=  r[i] * A[i][j] * r[j];
		    }
		}	    
				
		// new residual
        if(lower != 0)
        {
    	  rk=upper/lower;
        }
        else
        {
        	std::cerr << "R^T*A*r is ZERO!! Aborting!!" << "\n";
//        	exit(-1);
        }

        // update the solution
        for(unsigned int i = 0; i < N; i++)
        {
    	  x[i] += rk * r[i];
        }

        // computing the norm of the residual.
        for(int i=0;i<N;i++)
        {
        	sum = sum+abs(r[i]);
        }
        // computing the norm of the residual.
        if( sum <= 10e-5)
        {
    	std::cout <<"n step to converge:(gradient)  " << step ; 
    	std::cout <<"\n";
    	break;
        }
    }
}

double CG(double **A, double *b, double *x, int N, int M)
{
	double dot[N];
	double r[N],p0[N];
	double A_p0[N],beta,rk1;

	// zero the dot array
	for(unsigned int i=0; i < N; i++)
	{
		dot[i]=0.0;
	}

	 // compute A * x0
	for(unsigned int i = 0; i < N; i++)
	{
	    for(unsigned int j = 0; j < M; j++)
		{
			dot[i] += A[i][j] * x[j];
		}
//		cout<<dot[i]<<" ";
//		std::cout << "\n";
	}
	// compute residual
	for(unsigned int i = 0; i < N; i++)
	{
	    r[i] = b[i] - dot[i];
//		std::cerr << r[i]<<" ";
	}
//	cout<<"\n";
	for(int i=0; i<N; i++)
    {
    	p0[i]=r[i];
//    	cout<<p0[i]<<" ";
    }   
	for(unsigned int step=0; step<500; step++)
	{
		
//		double r_tran[1][2]={0.0,0.0};
//		double r_tran_A[1][2]={0.0,0.0};
        double upper1=0.0;
        double upper2=0.0;
        double lower1=0.0;
        double lower2=0.0;
        double sum=0.0;
        for(unsigned int i = 0; i < N; i++)
	    {
			upper1 += p0[i] * r[i];
//			cout<<upper1<<"\n";
		}
//		

	    // r^T*A*r
		for(unsigned int i = 0; i < N; i++)
	    {
	    	for(unsigned int j = 0; j < M; j++)
			{
				lower1 +=  p0[i] * A[i][j] * p0[j];
		    }
		}	    
//		cout<<lower1<<"\n";		
		// new residual
        if(lower1 != 0)
        {
    	  rk1=upper1/lower1;
//    	  cout<<rk1<<"\n";
        }
        else
        {
        	std::cerr << "R^T*A*r is ZERO!! Aborting!!" << "\n";
//        	exit(-1);
        }

        // update the solution
        for(unsigned int i = 0; i < N; i++)
        {
    	  x[i] += rk1 * p0[i];
//    	  cout<<x[i]<<" ";
        }
//        cout<<"\n";
        for(unsigned int i=0; i < N; i++)
		{
			A_p0[i]=0.0;
		}
        //update the residual 
        for(unsigned int i = 0; i < N; i++)
	    {
	    	for(unsigned int j = 0; j < M; j++)
			{
				A_p0[i] += A[i][j] * p0[j];
		    }
		}
		for(int i=0;i<N;i++)
        {
        	r[i]=r[i]-rk1*A_p0[i];
//        	cout<<r[i]<<" ";
        }

		for(unsigned int i = 0; i < N; i++)
	    {
			upper2 += p0[i] * r[i];
		}
//		cout<<upper2<<"\n";
	    for(unsigned int i = 0; i < N; i++)
	    {
	    	for(unsigned int j = 0; j < M; j++)
			{
				lower2 +=  p0[i] * A[i][j] * p0[j];
		    }
		}
//		cout << lower2;	  
		if(lower2 != 0)
        {
    	  beta=upper2/lower2;
        }  
		//update p0;
		for(int i=0;i<N;i++)
        {
        	p0[i]=r[i]-beta*p0[i];
        }
        // computing the norm of the residual.
        for(int i=0;i<N;i++)
        {
        	sum = sum+abs(r[i]);
        }
//        cout<<sum<<"\n";
        if( sum <= 10e-5)
        {
    	  std::cout <<"n step to converge(CG method):  " << step ; 
    	  std::cout <<"\n";
    	  break;
        }
    }
}
double Jordan(double **A, double *b, double *x, int N, int M)
{
	for(int step=0; step<500; step++)
	{
		int i,j;
		double sum1=0.0,dot[N],r[N];
		for(unsigned int i=0; i < N; i++)
		{
			dot[i]=0.0;
		}
		
		for(i=0; i<N; i++)
		{
			double sum=0.0;
			for( j=0; j<M; j++)
			{
				if(j!= i)
				{
					sum = sum + A[i][j]*x[j];
				}
			}
			if(A[i][i]!=0)
			{	
				x[i]=(1./A[i][i])*(b[i]-sum);
			}
//			cout<<x[i]<<" ";
		}
		for(unsigned int i = 0; i < N; i++)
	    {
	    	for(unsigned int j = 0; j < M; j++)
			{
				dot[i] += A[i][j] * x[j];
//			    std::cout << dot[i] <<" ";
		    }
//		    std::cout << "\n";
		}

		// compute residual
	    for(unsigned int i = 0; i < N; i++)
	    {
	    	r[i] = b[i] - dot[i];
//		    std::cerr << r[i] << "\n";
	    }
//		cout<<"\n";
	    for(int i=0;i<N;i++)
        {
        	sum1 = sum1+abs(r[i]);
        }
//        cout << sum1 <<"\n";
        // computing the norm of the residual.
        if( sum1 <= 10e-5)
        {
    	std::cout <<"n step to converge:(Jacobi)  " << step ; 
    	std::cout <<"\n";
    	break;
        }
	}
}

double Gauss_seidel(double **A, double *b, double *x, int N, int M)
{	
//	double x1[200][N];
	for(int step=0; step<500; step++)
	{	
		double sum1=0.0,dot[N],r[N];
		for(unsigned int i=0; i < N; i++)
		{
			dot[i]=0.0;
		}
		for(int i=0; i<N; i++)
		{
			double sum_before = 0.0;
			double sum_after = 0.0;
			for(int j=0; j<M; j++)
			{
				if(j<i)
				{
					sum_before=sum_before + A[i][j]*x[j];
				}
				if(j>i)
				{
					sum_after=sum_after + A[i][j]*x[j];
				}
			}
			x[i]=(1./A[i][i])*(b[i]-sum_before-sum_after);
//			cout<<x[i]<<" ";
		}
//		cout<<"\n";
		for(unsigned int i = 0; i < N; i++)
	    {
	    	for(unsigned int j = 0; j < M; j++)
			{
				dot[i] += A[i][j] * x[j];
//			    std::cout << dot[i] <<" ";
		    }
//		    std::cout << "\n";
		}

		// compute residual
	    for(unsigned int i = 0; i < N; i++)
	    {
	    	r[i] = b[i] - dot[i];
//		    std::cerr << r[i] << "\n";
	    }
//		cout<<"\n";
	    for(int i=0;i<N;i++)
        {
        	sum1 = sum1+abs(r[i]);
        }
//        cout << sum1 <<"\n";
        // computing the norm of the residual.
        if( sum1 <= 10e-5)
        {
    	std::cout <<"n step to converge:(Gauss_seidel)  " << step ; 
    	std::cout <<"\n";
    	break;
    	}
	}
}



int main()
{	
	double *b, *x;
    int N=10, M=10,n,m;
//	A = (double **) malloc(N*sizeof(double *));
//	b = (double *) malloc(N*sizeof(double *));
//	x = (double *) malloc(N*sizeof(double *));
    double **A = new double *[N];
	for(int i=0; i<N;i++){
		A[i] = new double[M];
	}
	b = new double[N];
	x = new double[N];
//    double A[N][M],b[N], x[N];
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
//// each time only use "ONE SOLVER"
//	grd( A, b, x, N, M);
//	CG( A, b, x, N, M);
//	Jordan(A, b, x, N, M);
	Gauss_seidel(A,b,x,N,M);
	for(int i=0; i<N; i++){
    	delete [] A[i];
    }
    delete [] A;
    delete [] b;
    delete [] x;
}

