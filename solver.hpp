fndef DYNEARTHSOL3D_solver_HPP
#define DYNEARTHSOL3D_solver_HPP

void grd( double **A, double *b, double *x, int N, int M);
void CG(double **A, double *b, double *x, int N, int M);
void Jordan(double **A, double *b, double *x, int N, int M);
void Gauss_seidel(double **A, double *b, double *x, int N, int M);

