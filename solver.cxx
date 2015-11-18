#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace std;
void grd(double **A, double *b, double *x, unsigned int N, unsigned int M)
{
	double dot[N];
	double r[N];
	double upper=0.0, lower=0.0, rk, tol=10e-5;
//	for(unsigned int i = 0; i < N; i++)
//	{
//		x[i] = 1.0;
//	}
	for(unsigned int step=0; step<100000; step++)
	{
		
//		double r_tran[1][2]={0.0,0.0};
//		double r_tran_A[1][2]={0.0,0.0};
        double upper=0.0;
        double lower=0.0;
        double sum=0.0;
	double b1[N];
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

		// compute resid
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
//        		exit(-1);
        	}

        // update the solution
        	for(unsigned int i = 0; i < N; i++)
        	{
    	  		x[i] += rk * r[i];
        	}

        // computing the norm of the residual.
        	for(unsigned int i=0;i<N;i++)
        	{
        		sum = sum+pow(r[i],2);
        	}
		cout<<"the residual is : " << sqrt(sum) << "\n";
        // computing the norm of the residual.
        	if( sqrt(sum) <= tol)
        	{	
			cout<<"the  input b is : "<<"\n";
			for(unsigned int i=0; i<N; i++)
			{
				cout<<b[i]<<" ";
			}
			cout<<"\n";
			cout<<"the estimated b value: "<<"\n";
			for(unsigned int i=0; i<N; i++)
			{
				cout<< dot[i]<<" ";
			}
			cout<<"\n";
    			std::cout <<"n step to converge:(gradient)  " << step <<"\n"; 
    			std::cout <<"the tolerance is : " << tol << "\n";
    			break;
        	}
	}
}

void CG(double **A, double *b, double *x, unsigned int N,unsigned int M)
{
	double dot[N];
	double r[N],p0[N];
	double A_p0[N],beta,rk1;
	double tol=10e-5;
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
	for(unsigned int i=0; i<N; i++)
	{
    		p0[i]=r[i];
//    		cout<<p0[i]<<" ";
    	}   
	for(unsigned int step=0; step<1000000; step++)
	{
		
//		double r_tran[1][2]={0.0,0.0};
//		double r_tran_A[1][2]={0.0,0.0};
        double upper1=0.0;
        double upper2=0.0;
        double lower1=0.0;
        double lower2=0.0;
        double sum=0.0, b1[N];
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
//    	  		cout<<rk1<<"\n";
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
//    	  		cout<<x[i]<<" ";
        	}
//        cout<<"\n";
		for(unsigned int i=0; i<N; i++)
		{
			b1[i]=0.0;
		}
		for(unsigned int i=0; i<N; i++)
		{
			for(unsigned j=0; j<M; j++)
			{
				b1[i] += A[i][j]*x[j];
			}
		}
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
		for(unsigned int i=0;i<N;i++)
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
		for(unsigned int i=0;i<N;i++)
        	{
        		p0[i]=r[i]-beta*p0[i];
        	}
        // computing the norm of the residual.
        	for(unsigned int i=0;i<N;i++)
        	{
        		sum = sum+pow(r[i],2);
        	}
        	cout<<sqrt(sum)<<"\n";
        	if( sqrt(sum) <= tol)
        	{
			cout<<"the input b is : "<<"\n";
                        for(unsigned int i=0; i<N; i++)
                        {   
                                cout<<b[i]<<" ";
                        }   
                        cout<<"\n";
                        cout<<"the estimated b value: "<<"\n";
                        for(unsigned int i=0; i<N; i++)
                        {   
                                cout<< b1[i]<<" ";
                        }   
                        cout<<"\n";			
    	  		std::cout <<"n step to converge(CG method):  " << step <<"\n" ; 
    	  		std::cout <<"the tolerance is : " << tol << "\n" ;
    	  		break;
        	}	
    	}
}
void Jacobi(double **A, double *b, double *x, unsigned int N,unsigned int M)
{
	double tol=10e-5,x1[N];
	for(unsigned int step=0; step<100000; step++)
	{
		double sum1=0.0,dot[N],r[N];
		for(unsigned int i=0; i < N; i++)
		{
			dot[i]=0.0;
		}
		
		for(unsigned int i=0; i<N; i++)
		{
			double sum=0.0;
			for(unsigned int j = 0; j < M; j++)
			{
				if(j!= i)
				{
					sum = sum + A[i][j]*x[j];
				}
			}
			if(A[i][i]!=0)
			{	
				x1[i]=(1./A[i][i])*(b[i]-sum);
			}
//			cout<<x[i]<<" ";
		}
		for(unsigned int i = 0; i< N; i++)
		{
		x[i]=x1[i];
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
	    	for(unsigned int i=0;i<N;i++)
        	{
        		sum1 = sum1+pow(r[i],2);
        	}
	        cout <<"the residual is :"<< sqrt(sum1) <<"\n";
        // computing the norm of the residual.
        	if( sqrt(sum1) <= tol)
        	{
			cout<<"the input b is : "<<"\n";
                        for(unsigned int i=0; i<N; i++)
                        {   
                                cout<<b[i]<<" ";
                        }   
                        cout<<"\n";
                        cout<<"the estimated b value: "<<"\n";
                        for(unsigned int i=0; i<N; i++)
                        {   
                                cout<< dot[i]<<" ";
                        }   
                        cout<<"\n";
    			std::cout <<"n step to converge:(Jacobi)  " << step<<"\n" ; 
    			std::cout <<"the tolerance is:"<< tol <<"\n";
    			break;
        	}
	}
}

void Gauss_seidel(double **A, double *b, double *x, unsigned int N,unsigned int M)
{	
//	double x1[200][N];
	double tol=10e-5;
	for(unsigned int step=0; step<100000; step++)
	{	
		double sum1=0.0,dot[N],r[N];
		for(unsigned int i=0; i < N; i++)
		{
			dot[i]=0.0;
		}
		for(unsigned int i=0; i<N; i++)
		{
			double sum_before = 0.0;
			double sum_after = 0.0;
			for(unsigned int j=0; j<M; j++)
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
//		    	std::cerr << r[i] << "\n";
	    	}
//		cout<<"\n";
	    	for(unsigned int i=0;i<N;i++)
        	{
        		sum1 = sum1+pow(r[i],2);
        	}
	        cout << "the residual is: "<< sqrt(sum1) <<"\n";
        // computing the norm of the residual.
        	if( sqrt(sum1) <= tol)
        	{
    			cout<<"the input b is : "<<"\n";
                        for(unsigned int i=0; i<N; i++)
                        {   
                                cout<<b[i]<<" ";
                        }   
                        cout<<"\n";
                        cout<<"the estimated b value: "<<"\n";
                        for(unsigned int i=0; i<N; i++)
                        {   
                                cout<< dot[i]<<" ";
                        }   
                        cout<<"\n";
			std::cout <<"n step to converge:(Gauss_seidel)  " << step<<"\n" ; 
    			std::cout <<"the tolerance is :"<< tol <<"\n";
    			break;
    		}
	}
}
