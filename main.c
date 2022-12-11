#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double *alocaVetor(int tam) { return (double *)calloc(tam, sizeof(double)); }

double p(double x)
{
  return -2/x;
}

double q(double x)
{
  return 2/(x*x);
}

double r(double x)
{
  return sin(log(x))/(x*x);
}

double *DiferencaFinita(double a, int N, double h, double alpha, double beta){
  double *D, *U, *L, *R, *y, *l, *u, *z, t=a+h, h2;
  int i;

  h2= h*h;

  D = alocaVetor(N + 1);
  U = alocaVetor(N + 1);
  L = alocaVetor(N + 1);
  u = alocaVetor(N + 1);
  l = alocaVetor(N + 1);
  z = alocaVetor(N + 1);
  y = alocaVetor(N + 2);
  R = alocaVetor(N + 1);

  D[1]= 2+h2*q(t);
  U[1]= -1 + (h/2)*p(t);
  R[1]= -h2*r(t) + (1-(h/2)*p(t))*alpha;

  for(i=2; i<N; i++){
    t+=h;
    D[i]= 2+h2*q(t);
    U[i]= -1+(h/2)*p(t);
    L[i]= -1 - (h/2)*p(t);
    R[i]= -h2*r(t);
  }

  t+=h;
  D[N]= 2+h2*q(t);
  L[N]= -1 - (h/2)*p(t);
  R[N]= -h2*r(t) + (1- (h/2)*p(t))*beta;

  l[1]= D[1];
  u[1]= U[1]/D[1];
  z[1]= R[1]/l[1];

  for(i=2; i<N; i++){
    l[i]= D[i] - L[i]*u[i-1];
    u[i]= U[i]/l[i];
    z[i]= (R[i] - L[i]*z[i-1])/l[i];
  }
  l[N]= D[N] - L[N]*u[N-1];
  z[N]= (R[N] - L[N]*z[N-1])/l[N];
  y[N]= z[N];
  y[0]= alpha;
  y[N+1] = beta;

  for(i=N-1; i>0; i--){
    y[i]= z[i]-u[i]*y[i+1];
  return y;
  }

}

int main(int argc, char **argv) {
  double *x, a, b, h, t, alpha, beta;
  int N1, N2, N3, i;

  a= atof(argv[1]);
  b= atof(argv[2]);
  N1= atoi(argv[3]);
  N2= atoi(argv[4]);
  N3= atoi(argv[5]);
  alpha= atof(argv[6]);
  beta= atof(argv[7]);
  
  h= (b-a)/(N1-1);

  x= DiferencaFinita(a, N1, h, alpha, beta);
  
  for(i=0; i<N1+1; i++){
    printf("%lf %lf\n", a+i*h, x[i]);
  }

  h= (b-a)/(N2-1);

  x= DiferencaFinita(a, N2, h, alpha, beta);
  
  for(i=0; i<N2+1; i++){
    printf("%lf %lf\n", a+i*h, x[i]);
  }

  h= (b-a)/(N3-1);

  x= DiferencaFinita(a, N3, h, alpha, beta);
  
  for(i=0; i<N3+1; i++){
    printf("%lf %lf\n", a+i*h, x[i]);
  }
  

  return 0;
}