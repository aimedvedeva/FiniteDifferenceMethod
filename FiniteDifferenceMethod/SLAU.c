#include <stdlib.h>
#include "SLAU.h"

double *passSLAU(double *b, double *c, double *d, double *r, int matrSize) // метод прогонки для 3-ёх диагональной матрицы
{
  int i = 0;
  double *x = NULL, *lambda = NULL, *delta = NULL;

  if ((x = malloc(sizeof(double)*matrSize)) == NULL)
    return NULL;
  if ((lambda = malloc(sizeof(double)*matrSize)) == NULL)
  {
    free(x);
    return NULL;
  }
  if ((delta = malloc(sizeof(double)*matrSize)) == NULL)
  {
    free(x);
    free(lambda);
    return NULL;
  }

  delta[0]  = -d[0] / c[0];
  lambda[0] =  r[0] / c[0];
  for (i = 1; i < matrSize; i++)
  {
    delta[i]  = -d[i] / (c[i] + b[i] * delta[i - 1]);
    lambda[i] = (r[i] - b[i] * lambda[i - 1]) / (c[i] + b[i] * delta[i - 1]);
  }

  x[matrSize - 1] = lambda[matrSize - 1];
  for (i = matrSize - 2; i >= 0; i--)
  {
    x[i] = delta[i] * x[i + 1] + lambda[i];
  }

  free(lambda);
  free(delta);
  return x;
}