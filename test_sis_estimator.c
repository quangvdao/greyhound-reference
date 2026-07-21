#include <math.h>
#include <stdio.h>
#include "labrador.h"

typedef struct {
  size_t rank;
  size_t width;
  double bound;
  size_t beta;
  size_t dimension;
  double bits;
} test_vector;

int main(void) {
  size_t i,j,k;
  const test_vector vectors[] = {
    { 21,1700,499535789,507,2976,134.355 },
    { 14,592,16228284,486,2394,128.790 },
    { 11,250,1707129,522,2176,138.330 }
  };

  sis_set_security_mode(SIS_SECURITY_INVALID);
  if(sis_secure(21,1700,499535789)) {
    fprintf(stderr,"invalid SIS security mode did not fail closed\n");
    return 1;
  }
  sis_set_security_mode(SIS_SECURITY_L2_QUANTUM128_ADPS16);
  for(i=0;i<sizeof(vectors)/sizeof(vectors[0]);i++) {
    sis_estimate estimate = sis_estimate_l2_core_svp_adps16(
      vectors[i].rank,vectors[i].width,vectors[i].bound);
    if(!estimate.valid || !estimate.finite || estimate.trivially_easy ||
       estimate.beta != vectors[i].beta ||
       estimate.lattice_dimension != vectors[i].dimension ||
       fabs(estimate.quantum_bits-vectors[i].bits) > 1e-9 ||
       !sis_secure(vectors[i].rank,vectors[i].width,vectors[i].bound)) {
      fprintf(stderr,"SIS estimator mismatch in vector %zu\n",i);
      return 2;
    }
  }
  if(sis_secure(1,1,ldexp(1,31)) || sis_secure(0,1,1)) {
    fprintf(stderr,"SIS estimator accepted invalid/trivially-easy input\n");
    return 3;
  }
  if(sis_secure(13,592,16228284)) {
    fprintf(stderr,"SIS estimator accepted a 116.07-bit instance\n");
    return 4;
  }
  {
    const size_t ranks[] = {1,4,9,11,13,21,64};
    const size_t widths[] = {1,48,250,592,1700,4095};
    const double bounds[] = {1,5000,1707129,16228284,499535789,ldexp(1,31)};
    for(i=0;i<sizeof(ranks)/sizeof(ranks[0]);i++)
      for(j=0;j<sizeof(widths)/sizeof(widths[0]);j++)
        for(k=0;k<sizeof(bounds)/sizeof(bounds[0]);k++) {
          sis_estimate estimate = sis_estimate_l2_core_svp_adps16(
            ranks[i],widths[j],bounds[k]);
          int expected = estimate.valid && !estimate.trivially_easy &&
                         (!estimate.finite || estimate.quantum_bits >= 128);
          if(sis_secure(ranks[i],widths[j],bounds[k]) != expected) {
            fprintf(stderr,"Fast SIS predicate mismatch at %zu/%zu/%.0f\n",
                    ranks[i],widths[j],bounds[k]);
            return 5;
          }
        }
  }
  puts("SIS ADPS16 L2 security tests passed");
  return 0;
}
