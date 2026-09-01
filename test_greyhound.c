#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "malloc.h"
#include "randombytes.h"
#include "fips202.h"
#include "labrador.h"
#include "chihuahua.h"
#include "pack.h"
#include "greyhound.h"

static double wall_time(void) {
  struct timespec t;

  clock_gettime(CLOCK_MONOTONIC,&t);
  return (double)t.tv_sec+1e-9*t.tv_nsec;
}

static void benchmark_seed(uint8_t seed[16], uint8_t domain) {
  static const uint8_t prefix[] = "GREYHOUND-BENCH-SEED-V1";
  const char *value = getenv("GREYHOUND_BENCH_SEED");
  shake128incctx ctx;

  if(value == NULL) {
    randombytes(seed,16);
    return;
  }
  shake128_inc_init(&ctx);
  shake128_inc_absorb(&ctx,prefix,sizeof(prefix)-1);
  shake128_inc_absorb(&ctx,(const uint8_t*)value,strlen(value));
  shake128_inc_absorb(&ctx,&domain,1);
  shake128_inc_finalize(&ctx);
  shake128_inc_squeeze(seed,16,&ctx);
}

static int test_polcom(size_t len) {
  int ret;
  int64_t x,y;
  polz *s = NULL;
  polcomctx ctx = {};
  polcomprf pi = {};
  prncplstmnt st = {};
  witness wt = {};
  __attribute__((aligned(16)))
  uint8_t seed[16];

  printf("Testing Greyhound polynomial commitment scheme\n\n");
  benchmark_seed(seed,0);
  s = _aligned_alloc(64,len*sizeof(polz));
  polzvec_almostuniform(s,len,seed,0);
  polzvec_center(s,len);

  x = 43;
  y = polzvec_eval(s,len,x);

  polcom_commit(&ctx,s,len);
  print_polcomctx_pp(&ctx);

  ret = polcom_eval(&wt,&pi,&ctx,x,y);
  if(ret) {
    printf("ERROR: Greyhound fold grind failed: %d\n",ret);
    goto end;
  }
  free(s);
  s = NULL;
  free_polcomctx(&ctx);
  print_polcomprf_pp(&pi);

  ret = polcom_reduce(&st,&pi);
  if(ret) {
    printf("ERROR: Reduction to Chihuahua statement failed: %d\n",ret);
    goto end;
  }
  free_polcomprf(&pi);
  print_prncplstmnt_pp(&st);
  print_witness_pp(&wt);

  ret = principle_verify(&st,&wt);
  if(ret) {
    printf("ERROR: Verification of Chihuahua statement failed: %d\n",ret);
    goto end;
  }

end:
  free(s);
  free_polcomctx(&ctx);
  free_polcomprf(&pi);
  free_prncplstmnt(&st);
  free_witness(&wt);
  return ret;
}

static int test_pack(size_t len) {
  int ret;
  int64_t x,y;
  double t;
  size_t wire_bytes;
  uint8_t *wire = NULL,*wire2 = NULL;
  polz *s;
  polcomctx ctx = {};
  polcomprf pi = {};
  polcomprf decoded_pi = {};
  composite p = {};
  composite decoded = {};
  __attribute__((aligned(16)))
  uint8_t seed[16];

  printf("Testing Greyhound Pack for degree 2^%.2g\n\n",log2(len << 6));
  benchmark_seed(seed,1);
  s = _aligned_alloc(64,len*sizeof(polz));
  polzvec_almostuniform(s,len,seed,0);
  polzvec_center(s,len);

  x = 43;
  y = polzvec_eval(s,len,x);

  t = wall_time();
  polcom_commit(&ctx,s,len);
  t = wall_time()-t;
  printf("Greyhound Pack commit wall time: %.4fs\n\n",t);
  print_polcomctx_pp(&ctx);

  ret = composite_prove_polcom(&p,&pi,&ctx,x,y);
  if(ret) goto end;
  wire_bytes = greyhound_pack_contextual_serialized_size(&pi,&p);
  wire = malloc(wire_bytes);
  wire2 = malloc(wire_bytes);
  if(!wire || !wire2 || greyhound_pack_serialize_contextual(wire,wire_bytes,&pi,&p) ||
     greyhound_pack_deserialize_contextual(&decoded_pi,&decoded,&pi,&p,wire,wire_bytes) ||
     greyhound_pack_serialize_contextual(wire2,wire_bytes,&decoded_pi,&decoded) ||
     memcmp(wire,wire2,wire_bytes)) {
    printf("ERROR: contextual Greyhound proof round-trip failed\n");
    ret = 1;
    goto end;
  }
  {
    polcomprf bad_pi = {};
    composite bad = {};

    size_t top_bytes = polcomprf_contextual_serialized_size(&pi);
    wire2[top_bytes] = 255;
    if(!greyhound_pack_deserialize_contextual(&bad_pi,&bad,&pi,&p,wire2,wire_bytes)) {
      printf("ERROR: invalid contextual Rice parameter was accepted\n");
      free_polcomprf(&bad_pi);
      free_composite(&bad);
      ret = 1;
      goto end;
    }
    wire2[top_bytes] = wire[top_bytes];
    memset(&wire2[8],0xff,4);
    if(!greyhound_pack_deserialize_contextual(&bad_pi,&bad,&pi,&p,wire2,wire_bytes)) {
      printf("ERROR: invalid Greyhound fold grind nonce was accepted\n");
      free_polcomprf(&bad_pi);
      free_composite(&bad);
      ret = 1;
      goto end;
    }
    memcpy(&wire2[8],&wire[8],4);
    if(!greyhound_pack_deserialize_contextual(&bad_pi,&bad,&pi,&p,wire2,wire_bytes-1)) {
      printf("ERROR: truncated contextual proof was accepted\n");
      free_polcomprf(&bad_pi);
      free_composite(&bad);
      ret = 1;
      goto end;
    }
  }
  printf("Contextual Greyhound proof size: %zu bytes\n\n",wire_bytes);
  ret = composite_verify_polcom(&decoded,&decoded_pi);
  if(ret) {
    printf("ERROR: verify_composite_polcom failed: %d\n",ret);
    goto end;
  }
  {
    uint32_t foldnonce = decoded_pi.foldnonce;
    decoded_pi.foldnonce = foldnonce == FOLD_GRIND_MAX_ATTEMPTS-1
                         ? foldnonce-1 : foldnonce+1;
    if(composite_verify_polcom(&decoded,&decoded_pi) == 0) {
      printf("ERROR: Modified Greyhound fold grind nonce was accepted\n");
      ret = 1;
      goto end;
    }
    decoded_pi.foldnonce = foldnonce;
  }
  ret = composite_verify_polcom(&decoded,&decoded_pi);
  if(ret)
    printf("ERROR: Restored Greyhound fold grind nonce did not verify: %d\n",ret);

end:
  free(wire);
  free(wire2);
  free(s);
  free_polcomctx(&ctx);
  free_polcomprf(&pi);
  free_polcomprf(&decoded_pi);
  free_composite(&p);
  free_composite(&decoded);
  return ret;
}

int main(int argc, char **argv) {
  int ret;
  size_t len;

  len = argc > 1 ? strtoull(argv[1],NULL,0) : (1 << 19);
  if(!len) {
    fprintf(stderr,"Usage: %s [number-of-64-coefficient-polynomials]\n",argv[0]);
    return 1;
  }

  if(!getenv("GREYHOUND_BENCH_PACK_ONLY")) {
    ret = test_polcom(len);
    if(ret) goto end;
  }
  ret = test_pack(len);
  if(ret) goto end;

end:
  free_comkey();
  return ret;
}
