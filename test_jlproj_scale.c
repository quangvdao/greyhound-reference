#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "aesctr.h"
#include "data.h"
#include "jlproj.h"
#include "malloc.h"
#include "poly.h"
#include "polx.h"

#define BENCH_REPETITIONS 3

static double wall_time(void) {
  struct timespec t;

  clock_gettime(CLOCK_MONOTONIC,&t);
  return (double)t.tv_sec+1e-9*t.tv_nsec;
}

static double median3(double a[BENCH_REPETITIONS]) {
  double t;

  if(a[0] > a[1]) { t=a[0]; a[0]=a[1]; a[1]=t; }
  if(a[1] > a[2]) { t=a[1]; a[1]=a[2]; a[2]=t; }
  if(a[0] > a[1]) { t=a[0]; a[0]=a[1]; a[1]=t; }
  return a[1];
}

int main(int argc, char **argv) {
  size_t i,len,matrix_bytes;
  double start;
  double derive_dense[BENCH_REPETITIONS],derive_ternary[BENCH_REPETITIONS];
  double project_dense[BENCH_REPETITIONS],project_ternary[BENCH_REPETITIONS];
  double collapse_dense[BENCH_REPETITIONS],collapse_ternary[BENCH_REPETITIONS];
  volatile uint64_t sink = 0;
  __attribute__((aligned(16)))
  uint8_t seed[16] = {1};
  __attribute__((aligned(64)))
  uint8_t challenge[256*QBYTES] = {2};
  aes128ctr_ctx aesctx;
  poly *witness;
  polx *collapsed;
  uint8_t *mat1,*mat2;
  __attribute__((aligned(64)))
  int32_t projection[256];

  if(argc != 2) {
    fprintf(stderr,"Usage: %s number-of-64-coefficient-polynomials\n",argv[0]);
    return 1;
  }
  len = strtoull(argv[1],NULL,0);
  if(len == 0 || len > SIZE_MAX/JL_MATRIX_POLY_BYTES) return 2;
  matrix_bytes = len*JL_MATRIX_POLY_BYTES;

  witness = _aligned_alloc(64,len*sizeof(poly));
  collapsed = _aligned_alloc(64,len*sizeof(polx));
  mat1 = _aligned_alloc(64,matrix_bytes);
  mat2 = _aligned_alloc(64,matrix_bytes);
  polyvec_ternary(witness,len,seed,0);

  for(i=0;i<BENCH_REPETITIONS;i++) {
    aes128ctr_init(&aesctx,seed,i);
    start = wall_time();
    aes128ctr_squeezeblocks(mat1,matrix_bytes/AES128CTR_BLOCKBYTES,&aesctx);
    derive_dense[i] = wall_time()-start;
    sink ^= mat1[i];

    aes128ctr_init(&aesctx,seed,i);
    start = wall_time();
    aes128ctr_squeezeblocks(mat1,matrix_bytes/AES128CTR_BLOCKBYTES,&aesctx);
    aes128ctr_squeezeblocks(mat2,matrix_bytes/AES128CTR_BLOCKBYTES,&aesctx);
    derive_ternary[i] = wall_time()-start;
    sink ^= mat1[i] ^ mat2[i];

    memset(projection,0,sizeof(projection));
    start = wall_time();
    polyvec_jlproj_add(projection,witness,len,mat1);
    project_dense[i] = wall_time()-start;
    sink ^= jlproj_normsq(projection);

    memset(projection,0,sizeof(projection));
    start = wall_time();
    polyvec_jlproj_add_ternary(projection,witness,len,mat1,mat2);
    project_ternary[i] = wall_time()-start;
    sink ^= jlproj_normsq(projection);

    start = wall_time();
    polxvec_jlproj_collapsmat(collapsed,mat1,len,challenge);
    collapse_dense[i] = wall_time()-start;
    sink ^= (uint16_t)collapsed[i % len].vec[0].vec->c[0];

    start = wall_time();
    polxvec_jlproj_collapsmat_ternary(collapsed,mat1,mat2,len,challenge);
    collapse_ternary[i] = wall_time()-start;
    sink ^= (uint16_t)collapsed[i % len].vec[0].vec->c[0];
  }

  const double dd = median3(derive_dense);
  const double dt = median3(derive_ternary);
  const double pd = median3(project_dense);
  const double pt = median3(project_ternary);
  const double cd = median3(collapse_dense);
  const double ct = median3(collapse_ternary);
  printf("polynomials=%zu scalar_dimension=%zu dense_matrix_bytes=%zu ternary_matrix_bytes=%zu\n",
         len,len*N,matrix_bytes,2*matrix_bytes);
  printf("derive_dense=%.9f derive_ternary=%.9f ratio=%.4f\n",dd,dt,dt/dd);
  printf("project_dense=%.9f project_ternary=%.9f ratio=%.4f\n",pd,pt,pt/pd);
  printf("collapse_dense=%.9f collapse_ternary=%.9f ratio=%.4f\n",cd,ct,ct/cd);
  if(sink == UINT64_MAX) fprintf(stderr,"benchmark sink: %llu\n",(unsigned long long)sink);

  free(witness);
  free(collapsed);
  free(mat1);
  free(mat2);
  return 0;
}
