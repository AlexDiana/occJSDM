#include <Accelerate/Accelerate.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9*t.tv_nsec; }
int main(int argc, char** argv){
  int sizes[4][3] = {{100,25,2000},{300,65,200},{900,185,50},{900,305,20}};
  for (int k = 0; k < 4; ++k) {
    int n = sizes[k][0], q = sizes[k][1], reps = sizes[k][2];
    double* Z = malloc(sizeof(double)*n*q); double* Q = calloc(q*q, sizeof(double));
    for (long i = 0; i < (long)n*q; ++i) Z[i] = (double)rand()/RAND_MAX - 0.5;
    double t0 = now();
    for (int r = 0; r < reps; ++r) cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, q, q, n, 1.0, Z, n, Z, n, 0.0, Q, q);
    double t1 = now();
    for (int r = 0; r < reps; ++r) cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, q, n, 1.0, Z, n, 0.0, Q, q);
    double t2 = now();
    printf("Accelerate n=%d q=%d: dgemm %.3f ms, dsyrk %.3f ms\n", n, q, 1000*(t1-t0)/reps, 1000*(t2-t1)/reps);
    free(Z); free(Q);
  }
  return 0;
}
