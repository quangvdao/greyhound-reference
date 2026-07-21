#ifndef GREYHOUND_H
#define GREYHOUND_H

#include <stdint.h>
#include <stddef.h>
#include "polx.h"
#include "poly.h"
#include "polz.h"
#include "labrador.h"
#include "chihuahua.h"

typedef struct {
  size_t len;
  size_t m;
  size_t n;
  comparams cpp[1];
  const polz *s;
  polx *sx;
  poly *t;
  polz *u1;
  uint8_t h[16];
  uint64_t normsq;  // expected norm
} polcomctx;

typedef struct {
  size_t len;
  size_t m;
  size_t n;
  int64_t x;
  int64_t y;
  comparams cpp[1];
  polz *u1;
  polz *u2;
  uint32_t foldnonce;
  uint64_t normsq;
} polcomprf;

void free_polcomctx(polcomctx *ctx);
void free_polcomprf(polcomprf *pi);

void print_polcomctx_pp(const polcomctx *ctx);
double print_polcomprf_pp(const polcomprf *pi);

int polcom_commit(polcomctx *ctx, const polz *s, size_t len);
int64_t polzvec_eval(const polz *a, size_t len, int64_t x);
int polcom_eval(witness *wt, polcomprf *pi, const polcomctx *ctx, int64_t x, int64_t y);
int polcom_reduce(prncplstmnt *st, const polcomprf *pi);

size_t polcomprf_serialized_size(const polcomprf *pi);
int polcomprf_serialize(uint8_t *out, size_t outlen, const polcomprf *pi);
int polcomprf_deserialize(polcomprf *pi, const uint8_t *in, size_t inlen);
size_t polcomprf_contextual_serialized_size(const polcomprf *pi);
int polcomprf_serialize_contextual(uint8_t *out, size_t outlen, const polcomprf *pi);
int polcomprf_deserialize_contextual(polcomprf *pi, const polcomprf *context,
                                     const uint8_t *in, size_t inlen,
                                     size_t *consumed);

#endif
