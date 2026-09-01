#ifndef LABRADOR_H
#define LABRADOR_H

#include <stdint.h>
#include <stddef.h>
#include "poly.h"
#include "polx.h"
#include "polz.h"
#include "sparsemat.h"

extern size_t comkey_len;
extern polx *comkey;

typedef struct {
  size_t deg;  // extension degree
  sparsemat a[1];
  polx *phi;
  polx *b;
} constraint;

typedef struct {
  size_t f;       // amortized opening decomposition parts
  size_t fu;      // uniform decomposition parts
  size_t fg;      // quadratic garbage decomposition parts
  size_t b;       // amortized opening decomposition basis
  size_t bu;      // uniform decomposition basis
  size_t bg;      // quadratic garbage decomposition basis
  size_t kappa;   // inner commitment rank
  size_t kappa1;  // outer commitment rank
  size_t u1len;
  size_t u2len;
} comparams;

typedef struct {
  size_t r;          // input witness multiplicity
  size_t *n;         // input witness ranks (r)
  size_t *nu;        // input witness decomposition parts (r)
  int tail;
  comparams cpp[1];  // commitment parameters
  polz *u1;          // outer commitment 1 (kappa1)
  size_t jlnonce;    // JL matrix nonce
  uint32_t foldnonce; // folding-challenge grind nonce
  int32_t p[256];    // JL projection
  polz *bb;          // int to pol lifting pols (LIFTS)
  polz *u2;          // outer commitment 2 (kappa1)
  uint64_t normsq;   // output witness norm
} proof;

typedef struct {
  size_t r;              // total amortized multiplicity
  size_t n;              // rank of amortized opening
  size_t m;              // rank of outer commitment openings (inner comms, garbage)
  int tail;              // true for no outer commitments
  const comparams *cpp;  // commitment parameters
  polx *u1;              // outer commitment 1 (kappa1) or inner commitments (r*kappa) if tail
  polx *u2;              // outer commitment 2 (kappa1) or garbage terms ((r*r+r)/2) if tail
  polx *c;               // challenges (r)
  constraint cnst[1];    // aggregated dot-product constraint
  uint64_t betasq;       // norm bound
  uint8_t h[16];         // hash
} statement;

typedef struct {
  size_t r;
  size_t *n;
  uint64_t *normsq;
  poly **s;
} witness;

typedef enum {
  SIS_SECURITY_LEGACY = 0,
  SIS_SECURITY_L2_QUANTUM128_ADPS16 = 1,
  SIS_SECURITY_INVALID = 2
} sis_security_mode;

#define SIS_MAX_RANK 128
#define FOLD_GRIND_MAX_ATTEMPTS 4096

typedef struct {
  double quantum_bits;
  size_t beta;
  size_t lattice_dimension;
  int valid;
  int finite;
  int trivially_easy;
} sis_estimate;

void sis_set_security_mode(sis_security_mode mode);
sis_security_mode sis_get_security_mode(void);
const char *sis_security_mode_name(void);
sis_estimate sis_estimate_l2_core_svp_adps16(size_t rank, size_t width,
                                             double l2_bound);
int sis_secure(size_t rank, size_t width, double l2_bound);
void print_sis_audit_pp(const char *role, size_t rank, size_t width, double l2_bound);
#define init_comkey NAMESPACE(init_comkey)
__attribute__((visibility("default")))
void init_comkey(size_t n);
#define free_comkey NAMESPACE(free_comkey)
__attribute__((visibility("default")))
void free_comkey(void);

int init_proof(proof *pi, const witness *wt, int quadratic, int tail);
void init_constraint_raw(constraint *cnst, size_t r, size_t n, size_t deg, int quadratic);
void init_constraint(constraint *cnst, const statement *st);
void init_statement(statement *st, const proof *pi, const uint8_t h[16]);
#define init_witness_raw NAMESPACE(init_witness_raw)
__attribute__((visibility("default")))
void init_witness_raw(witness *wt, size_t r, const size_t n[r]);
void witness_merge(witness *wt1, const witness *wt2);
#define set_witness_vector_raw NAMESPACE(set_witness_vector_raw)
__attribute__((visibility("default")))
int set_witness_vector_raw(witness *wt, size_t i, size_t n, size_t deg, const int64_t s[n*deg*N]);
void init_witness(witness *wt, const statement *st);

void free_proof(proof *pi);
void free_constraint(constraint *cnst);
void free_statement(statement *st);
#define free_witness NAMESPACE(free_witness)
__attribute__((visibility("default")))
void free_witness(witness *wt);

double print_proof_pp(const proof *pi);
void print_statement_pp(const statement *pi);
double print_witness_pp(const witness *wt);

/* Canonical proof wire format (little-endian, versioned, self-describing).
 * The size includes proof metadata, the 256 sparse-ternary JL coordinates, and every
 * fixed-width ring-element payload.  proof_deserialize expects a zeroed or
 * previously freed destination. */
size_t proof_serialized_size(const proof *pi);
int proof_serialize(uint8_t *out, size_t outlen, const proof *pi);
int proof_deserialize(proof *pi, const uint8_t *in, size_t inlen);

/* Tight context-dependent encoding.  Shape and parameter data are supplied by
 * a trusted schedule context and are not repeated on the wire. */
size_t proof_contextual_serialized_size(const proof *pi);
int proof_serialize_contextual(uint8_t *out, size_t outlen, const proof *pi);
int proof_deserialize_contextual(proof *pi, const proof *shape,
                                 const uint8_t *in, size_t inlen,
                                 size_t *consumed);

/* Canonical terminal-witness wire format.  Signed int16 coefficients are
 * encoded per witness part with the unique size-minimizing Golomb-Rice
 * parameter; decoding recomputes norms and rejects non-canonical encodings. */
size_t witness_serialized_size(const witness *wt);
int witness_serialize(uint8_t *out, size_t outlen, const witness *wt);
int witness_deserialize(witness *wt, const uint8_t *in, size_t inlen);
size_t witness_contextual_serialized_size(const witness *wt);
int witness_serialize_contextual(uint8_t *out, size_t outlen, const witness *wt);
int witness_deserialize_contextual(witness *wt, const witness *shape,
                                   const uint8_t *in, size_t inlen,
                                   size_t *consumed);

size_t commit_raw(polx *u, poly *t, size_t r, size_t n, const polx s[r][n],
                  size_t off, const comparams *cpp, int tail);
size_t qugarbage_raw(polx *u, poly *g, size_t r, size_t n, const polx s[r][n],
                     size_t off, const comparams *cpp, int tail);
void commit(statement *ost, witness *owt, proof *pi, polx sx[ost->r][ost->n], const witness *iwt);
void reduce_commit(statement *ost, const proof *pi);

int project(statement *ost, proof *pi, uint8_t *jlmat1, uint8_t *jlmat2,
            const witness *iwt);
int reduce_project(statement *ost, uint8_t *jlmat1, uint8_t *jlmat2,
                   const proof *pi, size_t r, uint64_t betasq);

void collaps_jlproj_raw(constraint *cnst, size_t r, size_t n, uint8_t h[16], const int32_t p[256],
                        const uint8_t *jlmat1, const uint8_t *jlmat2);
void collaps_jlproj(constraint *cnst, statement *st, const proof *pi,
                    const uint8_t *jlmat1, const uint8_t *jlmat2);
void lift_aggregate_zqcnst(statement *ost, proof *pi, size_t i, constraint *cnst, const polx sx[ost->r][ost->n]);
void reduce_lift_aggregate_zqcnst(statement *ost, const proof *pi, size_t i, const constraint *cnst);

int amortize(statement *ost, witness *owt, proof *pi, polx sx[ost->r][ost->n]);
int reduce_amortize(statement *ost, const proof *pi);

int prove(statement *ost, witness *owt, proof *pi, const statement *ist, const witness *iwt, int tail);
int reduce(statement *ost, const proof *pi, const statement *ist);
int verify(const statement *st, const witness *wt);

#endif
