#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "cpucycles.h"
#include "data.h"
#include "randombytes.h"
#include "jlproj.h"
#include "poly.h"
#include "polz.h"

int main(void) {
  size_t i;
  int64_t c;
  volatile uint64_t benchmark_sink = 0;
  unsigned long long t[21], overhead;
  __attribute__((aligned(16)))
  uint8_t seed[16];
  uint64_t nonce = 0;
  __attribute__((aligned(64)))
  uint8_t mat[JL_MATRIX_POLY_BYTES];
  __attribute__((aligned(64)))
  uint8_t mat2[JL_MATRIX_POLY_BYTES];
  __attribute__((aligned(64)))
  uint8_t buf[256*QBYTES];
  __attribute__((aligned(64)))
  int32_t p[256],p1[256],p2[256];
  poly s;
  polx r,r1,r2,tref,tmp;
  zz x,y;

  overhead = cpucycles_overhead();

  randombytes(seed,16);
  randombytes(mat,sizeof(mat));
  randombytes(mat2,sizeof(mat2));
  randombytes(buf,sizeof(buf));
  polyvec_uniform(&s,1,&primes[0],seed,nonce++);

  /* The optimized ternary projection must equal (S1*s + S2*s)/2 exactly. */
  memset(p,0,sizeof(p));
  memset(p1,0,sizeof(p1));
  memset(p2,0,sizeof(p2));
  poly_jlproj_add(p1,&s,mat);
  poly_jlproj_add(p2,&s,mat2);
  polyvec_jlproj_add_ternary(p,&s,1,mat,mat2);
  for(i=0;i<256;i++) {
    const int64_t sum = (int64_t)p1[i]+p2[i];
    if((sum & 1) != 0 || p[i] != sum/2) {
      fprintf(stderr,"ERROR: Ternary projection mismatch at coordinate %zu\n",i);
      return 1;
    }
  }

  /* Opposite sign planes cancel; identical planes recover the dense map. */
  memset(mat,0,sizeof(mat));
  memset(mat2,0xff,sizeof(mat2));
  memset(p,0,sizeof(p));
  polyvec_jlproj_add_ternary(p,&s,1,mat,mat2);
  for(i=0;i<256;i++) if(p[i] != 0) {
    fprintf(stderr,"ERROR: Opposite ternary sign planes did not cancel\n");
    return 2;
  }
  memset(mat2,0,sizeof(mat2));
  memset(p,0,sizeof(p));
  memset(p1,0,sizeof(p1));
  polyvec_jlproj_add_ternary(p,&s,1,mat,mat2);
  poly_jlproj_add(p1,&s,mat);
  if(memcmp(p,p1,sizeof(p)) != 0) {
    fprintf(stderr,"ERROR: Identical ternary sign planes did not recover dense projection\n");
    return 3;
  }

  randombytes(mat,sizeof(mat));
  randombytes(mat2,sizeof(mat2));

  /* Retain the original dense collapse identity as a regression baseline. */
  memset(p1,0,sizeof(p1));
  poly_jlproj_add(p1,&s,mat);
  polxvec_jlproj_collapsmat(&r1,mat,1,buf);
  polx_poly_mul(&tmp,&r1,&s);
  polx_getcoeff(&x,&tmp,0);
  if(!zz_less_than(&x,&modulus.q))
    zz_sub(&x,&x,&modulus.q);
  c = jlproj_collapsproj(p1,buf);
  zz_fromint64(&y,c);
  if(!zz_less_than(&y,&modulus.q))
    zz_sub(&y,&y,&modulus.q);
  if(!zz_equal(&x,&y)) {
    fprintf(stderr,"ERROR: Dense constant coeff doesn't match\n");
    return 4;
  }

  memset(p,0,sizeof(p));
  polyvec_jlproj_add_ternary(p,&s,1,mat,mat2);
  printf("Ternary projected norm: %.2f; expected: %.2f\n",sqrt(jlproj_normsq(p)),sqrt(128)*polyvec_norm(&s,1));

  /* Independently compose the two existing dense collapse kernels modulo q. */
  polxvec_jlproj_collapsmat(&r1,mat,1,buf);
  polxvec_jlproj_collapsmat(&r2,mat2,1,buf);
  polx_add(&tmp,&r1,&r2);
  polx_scale(&tref,&tmp,((UINT64_C(1) << LOGQ)-QOFF+1)/2);
  polxvec_jlproj_collapsmat_ternary(&r,mat,mat2,1,buf);
  polx_sub(&tmp,&r,&tref);
  if(!polx_iszero(&tmp)) {
    fprintf(stderr,"ERROR: Ternary collapse differs from dense reference composition\n");
    return 5;
  }

  polx_poly_mul(&r,&r,&s);
  polx_getcoeff(&x,&r,0);
  if(!zz_less_than(&x,&modulus.q))
    zz_sub(&x,&x,&modulus.q);

  c = jlproj_collapsproj(p,buf);
  zz_fromint64(&y,c);
  if(!zz_less_than(&y,&modulus.q))
    zz_sub(&y,&y,&modulus.q);

  if(!zz_equal(&x,&y)) {
    fprintf(stderr,"ERROR: Ternary constant coeff doesn't match\n");
    return 6;
  }

  memset(p,0,sizeof(p));
  for(i=0;i<21;i++) {
    t[i] = cpucycles();
    poly_jlproj_add(p,&s,mat);
    benchmark_sink ^= (uint32_t)p[i % 256];
  }
  for(i=0;i<20;i++)
    printf("poly_jlproj_dense:    %2lu: %8lld\n", i, t[i+1] - t[i] - overhead);

  memset(p,0,sizeof(p));
  for(i=0;i<21;i++) {
    t[i] = cpucycles();
    polyvec_jlproj_add_ternary(p,&s,1,mat,mat2);
    benchmark_sink ^= (uint32_t)p[i % 256];
  }
  for(i=0;i<20;i++)
    printf("poly_jlproj_ternary:  %2lu: %8lld\n", i, t[i+1] - t[i] - overhead);

  for(i=0;i<21;i++) {
    t[i] = cpucycles();
    polxvec_jlproj_collapsmat(&r,mat,1,buf);
    benchmark_sink ^= (uint16_t)r.vec[0].vec->c[i % N];
  }
  for(i=0;i<20;i++)
    printf("jlproj_collapsmat_dense:    %2lu: %8lld\n", i, t[i+1] - t[i] - overhead);

  for(i=0;i<21;i++) {
    t[i] = cpucycles();
    polxvec_jlproj_collapsmat_ternary(&r,mat,mat2,1,buf);
    benchmark_sink ^= (uint16_t)r.vec[0].vec->c[i % N];
  }
  for(i=0;i<20;i++)
    printf("jlproj_collapsmat_ternary:  %2lu: %8lld\n", i, t[i+1] - t[i] - overhead);

  if(benchmark_sink == UINT64_MAX)
    fprintf(stderr,"benchmark sink: %llu\n",(unsigned long long)benchmark_sink);

  return 0;
}
