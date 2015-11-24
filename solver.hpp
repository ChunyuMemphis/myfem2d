#ifndef DYNEARTHSOL3D_solver_HPP
#define DYNEARTHSOL3D_solver_HPP

double grd( double **A, double *b, double *x, unsigned int N,unsigned int M);
double CG(double **A, double *b, double *x,unsigned int N,unsigned int M);
double Jacobi(double **A, double *b, double *x,unsigned int N,unsigned int M);
double Gauss_seidel(double **A, double *b, double *x,unsigned int N,unsigned int M);

#endif
