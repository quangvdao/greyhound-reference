#ifndef PACK_H
#define PACK_H

#include <stdint.h>
#include <stddef.h>
#include "poly.h"
#include "polx.h"
#include "polz.h"
#include "labrador.h"
#include "chihuahua.h"
#include "dachshund.h"
#include "greyhound.h"

typedef struct {
  size_t l;
  double size;
  double tail_size;
  proof *pi[16];
  witness owt;
} composite;

#define free_composite NAMESPACE(free_composite)
__attribute__((visibility("default")))
void free_composite(composite *proof);

int composite_prove_principle(composite *proof, const prncplstmnt *st, const witness *wt);
#define composite_prove_simple NAMESPACE(composite_prove_simple)
__attribute__((visibility("default")))
int composite_prove_simple(composite *proof, commitment *com, const smplstmnt *st, const witness *wt);
int composite_prove_polcom(composite *proof, polcomprf *ppi, polcomctx *ctx, uint32_t x, uint32_t y);

int composite_verify_principle(const composite *proof, const prncplstmnt *st);
#define composite_verify_simple NAMESPACE(composite_verify_simple)
__attribute__((visibility("default")))
int composite_verify_simple(const composite *proof, const commitment *com, const smplstmnt *st);
int composite_verify_polcom(const composite *proof, const polcomprf *ppi);

/* Canonical complete Greyhound pack: versioned envelope, section lengths,
 * top Greyhound proof, all LBRP folds, and the terminal witness. */
size_t greyhound_pack_serialized_size(const polcomprf *top, const composite *proof);
int greyhound_pack_serialize(uint8_t *out, size_t outlen, const polcomprf *top,
                             const composite *proof);
int greyhound_pack_deserialize(polcomprf *top, composite *proof,
                               const uint8_t *in, size_t inlen);

/* Tight context-dependent proof.  The public commitment and adaptive fold
 * schedule are supplied by the caller and are not repeated on the wire. */
size_t greyhound_pack_contextual_serialized_size(const polcomprf *top,
                                                  const composite *proof);
int greyhound_pack_serialize_contextual(uint8_t *out, size_t outlen,
                                        const polcomprf *top,
                                        const composite *proof);
int greyhound_pack_deserialize_contextual(polcomprf *top, composite *proof,
                                          const polcomprf *top_context,
                                          const composite *proof_context,
                                          const uint8_t *in, size_t inlen);

#endif
