#ifndef DYNEARTHSOL3D_solver_HPP
#define DYNEARTHSOL3D_solver_HPP

void grd( double **A, double *b, double *x, unsigned int N,unsigned int M);
void CG(double **A, double *b, double *x,unsigned int N,unsigned int M);
void Jacobi(double **A, double *b, double *x,unsigned int N,unsigned int M);
void Gauss_seidel(double **A, double *b, double *x,unsigned int N,unsigned int M);

#endif
