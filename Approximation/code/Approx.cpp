#ifndef _APPROX_CPP_
#define _APPROX_CPP_

#include "Approx.h"

double Approx::one(double _x)
{
  return 1;
}

double Approx::x(double _x)
{
  return _x;
}

double Approx::x2(double _x)
{
  return pow(_x,2);
}


double Approx::innerProd(double (*u)(double), double (*v)(double), double (*p)(double))
{
  double res = 0;
  int N = xdata.size();
  for(int i = 1; i <= N; i++)
  {
    res += u(xdata[i-1])*v(xdata[i-1])*p(xdata[i-1]);
  }

  return res;
}

double Approx::innerProd(double (*v)(double), double (*p)(double))
{
  double res = 0;
  int N = xdata.size();
  for(int i = 1; i <= N; i++)
  {
    res += (ydata[i-1]) * (v(xdata[i-1])) * (p(xdata[i-1]));
  }

  return res;
}

Approx::Approx(std::vector<double> _xdata, std::vector<double> _ydata)
{
  xdata = _xdata;
  ydata = _ydata;
}

void Approx::solveNormal()
{  
  std::vector<double (*)(double)> list{one, x, x2};
  int index = 0;
  double G[3*3];
  for(int i = 0; i < 3; i++)
  {
    c[i] = innerProd(list[i], one);
    for(int j = 0; j < 3; j++)
    {
      G[index] = innerProd(list[i], list[j], one);
      Gram[index] = innerProd(list[i], list[j], one);
      index++;
    }
  }
  
  lapack_int info, m, n, lda, ldb, nrhs;
  m = 3;
  n = 3;
  nrhs = 1;
  lda = 3;
  ldb = 1;
  
  int ipiv[3];
  // Solve Normal Equation
  info = LAPACKE_dgesv(LAPACK_ROW_MAJOR,n,nrhs,G,lda,ipiv,c,ldb);

  std::cout << "a_0 = " << c[0] << std::endl;
  std::cout << "a_1 = " << c[1] << std::endl;
  std::cout << "a_2 = " << c[2] << std::endl;

}

int Approx::solveQR()
{
  int M = xdata.size();
  double A[M*3];
  int j = 0;
  // Matrix for LS problem
  for(int i = 0; i < M*3; i+=3)
  {
    A[i] = pow(xdata[j],2);
    A[i+1] = xdata[j];
    A[i+2] = 1;
    j++;
  }

  lapack_int m, n, lda, ldb, nrhs, info;
  double tau[3];
  m = M;
  n = 3;
  lda = 3;
  nrhs = 1;
  ldb = 1;

  // Rewrites A
  info = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR,m,n,A,lda,tau);

  double R[3*3];
  // Extract the matrix R1 from A
  for(int i = 0; i < 9; i++)
  {
    R1[i] = A[i];
    R[i] = A[i];
    if(i == 3 || i == 6 || i == 7)
    {
      R1[i] = 0;
      R[i] = 0;
    }
  }

  // The RHS of the DLS problem
  double b[M];
  for(int i = 0; i < M; i++)
    b[i] = ydata[i];

  // Rewrites A = Q1 (m*n)
  info = LAPACKE_dorgqr(LAPACK_ROW_MAJOR,m,n,3,A,lda,tau);

  // The array for cblas_dgemv parameter 
  double c1[M];
  for(int i = 0; i < M; i++)
    c1[i] = 0.0;
  // do Q1**T * b (Rewrites c1)
  cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1.0, A, lda, b, 1, 0.0, c1, 1);

  double c2[3];
  // Extract the first three elements from c1
  for(int i = 0; i < 3; i++)
  {
    c2[i] = c1[i];
  }

  // Solve the equation R1 * a = c2
  int ipiv[3];
  info = LAPACKE_dgesv(LAPACK_ROW_MAJOR,n,nrhs,R,lda,ipiv,c2,ldb);

  
  std::cout << "a_0 = " << c2[2] << std::endl;
  std::cout << "a_1 = " << c2[1] << std::endl;
  std::cout << "a_2 = " << c2[0] << std::endl;

  return 0;
}

double Approx::condG()
{
  lapack_int info;
  double S[3]; // singular values
  double U[3]; // Left singular vector
  double VT[3]; // Right singular vector (transposed)
  info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR,'N',3,3,Gram,3,S,U,3,VT,3);
  // 2-norm of matrix = the largest singular value
  double norm = S[0];
  // 2-norm of the inverse = reciprocal of the smallest singular value
  double norminv = 1 / S[2]; 

  return norm*norminv;
}

double Approx::condR1()
{
  lapack_int info;
  double S[3]; // singular values
  double U[3]; // left singular vectors
  double VT[3]; // right singular vectors (transposed)
  info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR,'N',3,3,R1,3,S,U,3,VT,3);
  // 2-norm of matrix = the largest singular value
  double norm = S[0];
  // 2-norm of the inverse = reciprocal of the smallest singular value
  double norminv = 1 / S[2];

  return norm*norminv;
}

#endif
