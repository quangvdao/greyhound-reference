#include <stdint.h>
#include "data.h"
#include "poly.h"

static int16_t center_mod(int64_t a, int16_t p) {
  a %= p;
  if (a > p / 2) a -= p;
  if (a < -p / 2) a += p;
  return (int16_t)a;
}

static int16_t mul_mod(int16_t a, int16_t b, int16_t p) {
  return center_mod((int32_t)a * b, p);
}

static int16_t montgomery_mul(int16_t a, int16_t b, const pdata *prime) {
  int32_t c = (int32_t)a * b;
  int16_t t = (int16_t)c * prime->pinv;
  return (int16_t)((c - (int32_t)t * prime->p) >> 16);
}

void poly_ntt(poly *r, const poly *a, const pdata *prime) {
  if (r != a) polyvec_copy(r, a, 1);
  poly_ntt_ref(r, prime);
  /* The assembly's final butterfly multiplies every lane by prime->f. */
  poly_scale(r, r, prime->f, prime);
}

void poly_invntt(poly *r, const poly *a, const pdata *prime) {
  int16_t work[N];
  int16_t zetainv[N];

  for (size_t i = 0; i < N; i++) {
    work[i] = a->vec->c[i];
    /* zeta has order 2N, hence zeta^-i = -zeta^(N-i) for i != 0. */
    zetainv[i] = i ? -montgomery_mul(1, prime->twist64[N - i], prime) : 1;
  }

  /* Reverse the forward DIF butterflies without dividing by two.  The
   * resulting factor N is intentional and matches the assembly backend. */
  for (size_t len = 1; len <= N / 2; len <<= 1) {
    for (size_t start = 0; start < N; start += 2 * len) {
      size_t k = 0;
      for (size_t j = start; j < start + len; j++) {
        int16_t d = mul_mod(work[j + len], zetainv[k], prime->p);
        int16_t lo = work[j];
        work[j] = center_mod((int32_t)lo + d, prime->p);
        work[j + len] = center_mod((int32_t)lo - d, prime->p);
        k += N / len;
      }
    }
  }

  for (size_t i = 0; i < N; i++)
    r->vec->c[i] = mul_mod(work[i], zetainv[i], prime->p);
}
