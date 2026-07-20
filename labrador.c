#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "randombytes.h"
#include "malloc.h"
#include "aesctr.h"
#include "fips202.h"
#include "data.h"
#include "poly.h"
#include "polx.h"
#include "polz.h"
#include "sparsemat.h"
#include "jlproj.h"
#include "parallel.h"
#include "labrador.h"

size_t comkey_len = 0;
polx *comkey = NULL;

static size_t triangularidx(size_t i,size_t j,size_t r) {
  if(i>j) {  // swap
    j ^= i;
    i ^= j;
    j ^= i;
  }
  i = i*r - (i*i+i)/2 + j;
  return i;
}

static double sis_log2_cutoff(size_t rank) {
  double maxlog;

  maxlog = 2*sqrt(LOGQ*LOGDELTA*N)*sqrt(rank);
  return MIN(LOGQ,maxlog);
}

static sis_security_mode security_mode = (sis_security_mode)-1;

static double root_hermite_factor(size_t beta) {
  double b;

  if(beta <= 40) return 1.01295;
  b = beta;
  return pow(b/(2*M_PI*M_E)*pow(M_PI*b,1/b),1/(2*(b-1)));
}

static size_t beta_from_delta(double target) {
  size_t low = 40;
  size_t high = 1 << 16;

  if(!isfinite(target)) return 0;
  if(root_hermite_factor(40) < target) return 40;
  if(target < root_hermite_factor(high)) return 0;
  while(low < high) {
    size_t mid = low+(high-low)/2;
    if(root_hermite_factor(mid) <= target+1e-13)
      high = mid;
    else
      low = mid+1;
  }
  return low;
}

void sis_set_security_mode(sis_security_mode mode) {
  security_mode = mode;
}

sis_security_mode sis_get_security_mode(void) {
  const char *value;

  if(security_mode != (sis_security_mode)-1) return security_mode;
  value = getenv("LABRADOR_SIS_SECURITY");
  if(value && !strcmp(value,"l2-quantum128-adps16"))
    security_mode = SIS_SECURITY_L2_QUANTUM128_ADPS16;
  else if(value && strcmp(value,"legacy") && strcmp(value,"legacy-heuristic"))
    security_mode = SIS_SECURITY_INVALID;
  else
    security_mode = SIS_SECURITY_LEGACY;
  return security_mode;
}

const char *sis_security_mode_name(void) {
  sis_security_mode mode = sis_get_security_mode();

  if(mode == SIS_SECURITY_L2_QUANTUM128_ADPS16)
    return "l2-quantum128-adps16";
  if(mode == SIS_SECURITY_LEGACY)
    return "legacy-heuristic";
  return "INVALID";
}

sis_estimate sis_estimate_l2_core_svp_adps16(size_t rank, size_t width,
                                             double norm) {
  sis_estimate result = {};
  const double q = ldexp(1,LOGQ)-QOFF;
  const double logq = log2(q);
  double logbound,logdelta,droot,rootvolume,target;
  double lnq,loga2,logb2;
  size_t n,m,d,beta;

  if(!rank || !width || !isfinite(norm) || norm <= 0) return result;
  if(rank > SIZE_MAX/N || width > SIZE_MAX/N) return result;
  result.valid = 1;
  norm = ceil(norm);
  if(norm >= (q-1)/2) {
    result.trivially_easy = 1;
    result.finite = 1;
    return result;
  }
  n = rank*N;
  m = width*N;
  logbound = log2(norm);
  logdelta = logbound*logbound/(4*n*logq);
  droot = sqrt(n*logq/logdelta);
  d = droot >= (double)SIZE_MAX ? SIZE_MAX : (size_t)floor(droot);
  d = MIN(d,m);
  result.lattice_dimension = d;
  if(d <= 1) {
    result.valid = 0;
    return result;
  }

  rootvolume = (double)n/d*logq;
  target = exp2((logbound-rootvolume)/(d-1));
  beta = target >= 1 ? beta_from_delta(target) : 0;
  result.beta = beta ? MIN(beta,d) : d;

  lnq = logq*M_LN2;
  loga2 = log((double)n*lnq);
  logb2 = log((double)d)+2*(double)n/d*lnq;
  if(!beta || beta > d || log(norm) <= 0.5*MIN(loga2,logb2))
    return result;

  result.finite = 1;
  result.quantum_bits = 0.265*beta;
  return result;
}

int sis_secure(size_t rank, size_t width, double norm) {
  sis_estimate estimate;

  if(sis_get_security_mode() == SIS_SECURITY_LEGACY)
    return log2(norm) < sis_log2_cutoff(rank);
  if(sis_get_security_mode() != SIS_SECURITY_L2_QUANTUM128_ADPS16)
    return 0;
  estimate = sis_estimate_l2_core_svp_adps16(rank,width,norm);
  return estimate.valid && !estimate.trivially_easy &&
         (!estimate.finite || estimate.quantum_bits >= 128);
}

void print_sis_audit_pp(const char *role, size_t rank, size_t width, double l2_bound) {
  double log2_bound = log2(l2_bound);
  double cutoff = sis_log2_cutoff(rank);
  sis_estimate estimate = sis_estimate_l2_core_svp_adps16(rank,width,l2_bound);

  printf("    %-18s rank = %zu, width = %zu ring columns, L2 collision bound = %.9g\n",
         role,rank,width,l2_bound);
  printf("      Legacy cutoff: log2(B) %.3f < %.3f (%s; %.3f-bit margin)\n",
         log2_bound,cutoff,log2_bound < cutoff ? "passes" : "FAILS",cutoff-log2_bound);
  if(!estimate.valid)
    printf("      L2 quantum-128 ADPS16: invalid estimator input; FAILS\n");
  else if(estimate.trivially_easy)
    printf("      L2 quantum-128 ADPS16: trivially easy (B >= (q - 1)/2); FAILS\n");
  else if(!estimate.finite)
    printf("      L2 quantum-128 ADPS16: no finite Euclidean reduction found; passes 128-bit floor"
           " (d = %zu)\n",estimate.lattice_dimension);
  else
    printf("      L2 quantum-128 ADPS16: %.3f bits, beta = %zu, d = %zu (%s)\n",
           estimate.quantum_bits,estimate.beta,estimate.lattice_dimension,
           estimate.quantum_bits >= 128 ? "passes" : "FAILS");
  printf("      SIS-AUDIT role=%s rank=%zu width=%zu l2_bound=%.17g\n",
         role,rank,width,l2_bound);
}

typedef struct {
  size_t first;
  uint8_t seed[16];
} comkey_job;

static void init_comkey_chunk(size_t index, void *context) {
  const comkey_job *job = context;
  const size_t chunk = job->first+index;

  polxvec_almostuniform(&comkey[32*chunk],32,job->seed,chunk);
}

void init_comkey(size_t len) {
  size_t chunks;
  polx *buf;

  if(comkey_len >= len)
    return;

  chunks = (len+31)/32;
  buf = _aligned_alloc(64,chunks*32*sizeof(polx));
  if(comkey_len) polxvec_copy(buf,comkey,comkey_len);
  free(comkey);
  comkey = buf;

  comkey_job job = { .first = comkey_len/32, .seed = {} };
  parallel_for(chunks-job.first,32*N,init_comkey_chunk,&job);

  comkey_len = chunks*32;
}

void free_comkey(void) {
  free(comkey);
  comkey = NULL;
  comkey_len = 0;
}

int init_proof(proof *pi, const witness *wt, int quadratic, int tail) {
  size_t i,j,k,t,u;
  size_t nn,rr;
  int ret,decompose;
  double vars,varz,varg;
  void *buf;
  comparams *cpp = pi->cpp;

  pi->r = (quadratic == 2) ? 2*wt->r + 1 : wt->r;
  pi->tail = tail;

  uint64_t normsq[pi->r];
  buf = _malloc(2*pi->r*sizeof(size_t));
  pi->n   = (size_t*)buf;
  pi->nu  = (size_t*)&pi->n[pi->r];

  if(quadratic == 2) {
    for(i=0;i<wt->r;i++) {
      pi->n[i] = wt->n[i];
      pi->n[wt->r+1+i] = wt->n[i];  // sigmam1 / flip
      normsq[i] = wt->normsq[i];
      normsq[wt->r+1+i] = (TAU1+4*TAU2)*wt->normsq[i];
    }
    pi->n[wt->r] = (LOGQ+9)/10*wt->r;  // Zq to Rq liftings
    normsq[wt->r] = ldexp(1,20)/12*pi->n[wt->r]*N;

    for(i=0;i<pi->r;i++)
      pi->nu[i] = 0;
    pi->nu[wt->r-1] = -(size_t)1;
    pi->nu[wt->r] = -(size_t)1;
    pi->nu[2*wt->r] = -(size_t)1;
  }
  else {
    for(i=0;i<pi->r;i++) {
      pi->n[i] = wt->n[i];
      normsq[i] = wt->normsq[i];
      pi->nu[i] = (quadratic) ? -(size_t)1 : 0;
    }
    pi->nu[pi->r - 1] = -(size_t)1;
  }

  for(k=15;k>0;k--) {
    /* decomposition / joining in rank */
    nn = t = 0;
    for(i=0;i<pi->r;i++) {
      t += pi->n[i];
      if(pi->nu[i]) {
        pi->nu[i] = t;
        nn = MAX(nn,t);
        t = 0;
      }
    }
    nn = (nn+k-1)/k;
    for(i=0;i<pi->r;i++)
      if(pi->nu[i])
        pi->nu[i] = (pi->nu[i]+nn-1)/nn;

    rr = 0;
    for(i=0;i<pi->r;i++)
      rr += pi->nu[i];

    /* decomposition in width */
    varz = 0;
    for(i=0;i<pi->r;i++)
      varz += normsq[i];
    varz /= nn*N;
    varz *= TAU1+4*TAU2;
    decompose = !tail && !sis_secure(13,nn,
                                      6*T*SLACK*sqrt(2*(TAU1+4*TAU2)*varz*nn*N));
    decompose = decompose || 64*varz > (1 << 28);
    if(decompose) {
      cpp->f = 2;
      cpp->b = round((log2(12)+log2(varz))/4);  // log2(sqrt(12*var))
    }
    else {
      cpp->f = 1;
      cpp->b = round((log2(12)+log2(varz))/2);
    }

    /* uniform decomposition */
    if(!tail) {
      cpp->fu = (LOGQ+2*cpp->b/3)/cpp->b;
      cpp->bu = (LOGQ+cpp->fu/2)/cpp->fu;
    }
    else {
      cpp->fu = 1;
      cpp->bu = LOGQ;
    }

    /* quadratic garbage decomposition */
    varg = 0;
    cpp->bg = cpp->b;
    if(!quadratic)
      cpp->fg = 0;
    else if(tail)
      cpp->fg = 1;
    else {
      t = u = 0;
      for(i=0;i<pi->r;i++) {
        vars = (double)normsq[i]/(pi->n[i]*N);
        j = pi->n[i];
        while(j >= nn-u) {
          j -= nn-u;
          t += vars*vars*(nn-u);
          varg = MAX(varg,t);
          t = u = 0;
        }
        t += vars*vars*j;
        u += j;
        if(pi->nu[i]) {
          varg = MAX(varg,t);
          t = u = 0;
        }
      }
      varg *= 2*N;
      cpp->fg = ceil((log2(12)+log2(varg))/(2*cpp->b));
      cpp->fg = MAX(1,cpp->fg);
    }

    /* commitment ranks */
    for(cpp->kappa=1;cpp->kappa<=SIS_MAX_RANK;cpp->kappa++) {
      pi->normsq = (ldexp(1,2*cpp->b)/12*(cpp->f-1) + varz/ldexp(1,2*cpp->b*(cpp->f-1)))*nn;
      if(!tail) {
        pi->normsq += (ldexp(1,2*cpp->bu)*(cpp->fu-1)
                       + ldexp(1,2*(LOGQ-(cpp->fu-1)*cpp->bu)))/12*(rr*cpp->kappa+(rr*rr+rr)/2);
      }
      if(!tail && quadratic)
        pi->normsq += (ldexp(1,2*cpp->bg)/12*(cpp->fg-1) + varg/ldexp(1,2*(cpp->fg-1)*cpp->bg))*(rr*rr+rr)/2;
      pi->normsq *= N;
      if(sis_secure(cpp->kappa,nn,
                    6*T*SLACK*ldexp(1,(cpp->f-1)*cpp->b)*sqrt(pi->normsq)))
        break;
    }

    if(!tail) {
      size_t tri = (rr*rr+rr)/2;
      size_t outer1_width = cpp->fu*rr*cpp->kappa+cpp->fg*tri;
      size_t outer2_width = cpp->fu*tri;
      for(cpp->kappa1=1;cpp->kappa1<=SIS_MAX_RANK;cpp->kappa1++)
        if(sis_secure(cpp->kappa1,outer1_width,2*SLACK*sqrt(pi->normsq)) &&
           sis_secure(cpp->kappa1,outer2_width,2*SLACK*sqrt(pi->normsq)))
          break;

      cpp->u1len = cpp->u2len = cpp->kappa1;
      if(cpp->kappa <= SIS_MAX_RANK && cpp->kappa1 <= SIS_MAX_RANK && cpp->fu*rr*cpp->kappa + (cpp->fu+cpp->fg)*(rr*rr+rr)/2 <= 1.1*nn)
        break;
    }
    else {
      cpp->kappa1 = 0;
      cpp->u1len = rr*cpp->kappa;
      if(quadratic) cpp->u1len += (rr*rr+rr)/2;
      cpp->u2len = 2*rr - 1;
      if(cpp->kappa <= SIS_MAX_RANK && (cpp->u1len + cpp->u2len)*LOGQ <= 1.1*nn*(log2(varz)/2+2.05))
        break;
    }
  }

  buf = _aligned_alloc(64,(cpp->u1len+cpp->u2len+LIFTS)*sizeof(polz));
  pi->u1  = (polz*)buf;
  pi->u2  = (polz*)&pi->u1[cpp->u1len];
  pi->bb  = (polz*)&pi->u2[cpp->u2len];

  if(cpp->kappa > SIS_MAX_RANK) {
    fprintf(stderr,"ERROR in init_proof(): Cannot make inner commitments secure!\n");
    ret = 1;
    goto err;
  }
  else if(cpp->kappa1 > SIS_MAX_RANK) {
    fprintf(stderr,"ERROR in init_proof(): Cannot make outer commitments secure!\n");
    ret = 2;
    goto err;
  }
  return 0;

err:
  free_proof(pi);
  return ret;

  /* solve m = fu*kappa*r + (fu+fg)*(r^2+r)/2 = a*r^2 + b*r = n_i/nu_i = nn/r *
   * where r = \sum_i nu_i, nn = \sum_i n_i, and n_i/nu_i = const             */
/*
  double a,b,y,m,r;
  a = (cpp->fu + cpp->fg)/2.0;
  b = cpp->fu*cpp->kappa + a;
  r = 64;
  for(i=0;i<10;i++) {
    y = (a*r+b)*r*r-nn;
    m = (3*a*r+2*b)*r;
    r -= y/m;
  }

  nn = ceil(nn/r);
  for(i=0;i<pi->r;i++) {
    pi->nu[i] = (pi->n[i]+nn/2)/nn;
    pi->nu[i] = MAX(1,pi->nu[i]);
  }
*/
}

void free_proof(proof *pi) {
  free(pi->n);  // one buffer for n,nu
  free(pi->u1);  // one buffer for u1,u2,bb
  pi->n = NULL;
  pi->u1 = NULL;
}

void init_constraint_raw(constraint *cnst, size_t r, size_t n, size_t deg, int quadratic) {
  size_t len;
  void *buf;

  cnst->deg = deg;

  buf = _aligned_alloc(64,r*extlen(n,deg)*sizeof(polx));
  cnst->phi = (polx*)buf;
  buf = _aligned_alloc(64,deg*sizeof(polx));
  cnst->b = (polx*)buf;

  cnst->a->len = 0;
  if(quadratic) {
    len = (r*r+r)/2;
    buf = _malloc(2*len*sizeof(size_t));
    cnst->a->rows = (size_t*)buf;
    cnst->a->cols = &cnst->a->rows[len];

    buf = _aligned_alloc(64,extlen(len,deg)*sizeof(polx));
    cnst->a->coeffs = (polx*)buf;
  }
  else {
    cnst->a->rows = NULL;
    cnst->a->cols = NULL;
    cnst->a->coeffs = NULL;
  }
}

void init_constraint(constraint *cnst, const statement *st) {
  init_constraint_raw(cnst,st->r,st->n,1,st->cpp->fg != 0);
}

void free_constraint(constraint *cnst) {
  free(cnst->phi);  // one buffer for phi (will be realloced)
  free(cnst->b);  // one buffer for b
  free(cnst->a->rows);  // one buffer for rows,cols
  free(cnst->a->coeffs);  // one buffer for coeffs (will be realloced)
  cnst->phi = cnst->b = cnst->a->coeffs = NULL;
  cnst->a->rows = NULL;
}

void init_statement(statement *st, const proof *pi, const uint8_t h[16]) {
  size_t i;
  size_t nn;
  void *buf;

  st->cpp = pi->cpp;
  st->tail = pi->tail;

  st->r = 0;
  for(i=0;i<pi->r;i++)
    st->r += pi->nu[i];

  nn = st->n = 0;
  for(i=0;i<pi->r;i++) {
    nn += pi->n[i];
    if(!pi->nu[i]) continue;
    nn = (nn+pi->nu[i]-1)/pi->nu[i];
    st->n = MAX(st->n,nn);
    nn = 0;
  }

  st->m  = st->r * st->cpp->fu * st->cpp->kappa;  // inner commitments
  st->m += (st->cpp->fu + st->cpp->fg)*(st->r*st->r + st->r)/2;  // garbage

  buf = _aligned_alloc(64,(st->cpp->u1len + st->cpp->u2len + st->r)*sizeof(polx));
  st->u1 = (polx*)buf;
  st->u2 = &st->u1[st->cpp->u1len];
  st->c = &st->u2[st->cpp->u2len];

  nn = st->r*extlen(st->cpp->fu*st->cpp->kappa,st->cpp->kappa1);
  nn += st->cpp->fg*(st->r*st->r + st->r)/2;
  nn = MAX(nn,st->cpp->fu*(st->r*st->r + st->r)/2);
  nn = MAX(nn,st->n);
  init_comkey(nn);
  init_constraint(st->cnst,st);
  memcpy(st->h,h,16);
}

void free_statement(statement *st) {
  free(st->u1);  // one buffer for u1,u2,c,b
  free_constraint(st->cnst);
  st->u1 = NULL;
}

void init_witness_raw(witness *wt, size_t r, const size_t n[r]) {
  size_t i, len;
  void *buf;

  wt->r = r;

  buf = _malloc(r*(sizeof(size_t) + sizeof(uint64_t) + sizeof(poly*)));
  wt->n = (size_t*)buf;
  wt->normsq = (uint64_t*)&wt->n[r];
  wt->s = (poly**)&wt->normsq[r];

  len = 0;
  for(i=0;i<r;i++) {
    wt->n[i] = n[i];
    len += n[i];
  }

  buf = _aligned_alloc(64,len*sizeof(poly));
  for(i=0;i<r;i++) {
    wt->s[i] = (poly*)buf;
    buf = (poly*)buf + n[i];
  }
}

void init_witness(witness *wt, const statement *st) {
  size_t i;
  const size_t r = st->cpp->f + !st->tail;
  size_t n[r];

  for(i=0;i<st->cpp->f;i++)
    n[i] = st->n;
  if(!st->tail)
    n[st->cpp->f] = st->m;

  init_witness_raw(wt,r,n);
}

void witness_merge(witness *wt1, const witness *wt2) {
  size_t i;
  void *buf;
  witness tmp[1];

  tmp->r = wt1->r + wt2->r;

  buf = _malloc(tmp->r*(sizeof(size_t) + sizeof(uint64_t) + sizeof(poly*)));
  tmp->n = (size_t*)buf;
  tmp->normsq = (uint64_t*)&tmp->n[tmp->r];
  tmp->s = (poly**)&tmp->normsq[tmp->r];

  for(i=0;i<wt1->r;i++) {
    tmp->n[i] = wt1->n[i];
    tmp->normsq[i] = wt1->normsq[i];
    tmp->s[i] = wt1->s[i];
  }
  for(i=0;i<wt2->r;i++) {
    tmp->n[wt1->r+i] = wt2->n[i];
    tmp->normsq[wt1->r+i] = wt2->normsq[i];
    tmp->s[wt1->r+i] = wt2->s[i];
  }

  free(wt1->n);
  *wt1 = *tmp;
}

int set_witness_vector_raw(witness *wt, size_t i, size_t n, size_t deg, const int64_t s[n*deg*N]) {
  if(i >= wt->r) {
    fprintf(stderr,"ERROR in set_witness_vector_raw(): Witness vector %zu does not exist\n",i);
    return 1;
  }
  if(n*deg != wt->n[i]) {
    fprintf(stderr,"ERROR in set_witness_vector_raw(): Mismatch of witness vector length\n");
    return 2;
  }

  polyvec_fromint64vec(wt->s[i],n,deg,s);
  wt->normsq[i] = polyvec_sprodz(wt->s[i],wt->s[i],n*deg);
  return 0;
}

void free_witness(witness *wt) {
  if(!wt->n) return;
  free(wt->s[0]);  // one buffer for *s
  free(wt->n);  // one buffer for n,s,normsq
  wt->n = NULL;
}

double print_proof_pp(const proof *pi) {
  size_t i,groups,row_width,outer1_width,outer2_width,tri;
  size_t wire_bytes,metadata_bytes,payload_bytes,accumulated;
  double s, projection_bits, modeled_bits,inner_norm,outer_norm;
  const comparams *cpp = pi->cpp;

  groups = 0;
  for(i=0;i<pi->r;i++) groups += pi->nu[i];
  row_width = accumulated = 0;
  for(i=0;i<pi->r;i++) {
    accumulated += pi->n[i];
    if(pi->nu[i]) {
      row_width = MAX(row_width,(accumulated+pi->nu[i]-1)/pi->nu[i]);
      accumulated = 0;
    }
  }
  tri = (groups*groups+groups)/2;
  outer1_width = groups*cpp->fu*cpp->kappa+cpp->fg*tri;
  outer2_width = cpp->fu*tri;
  wire_bytes = proof_contextual_serialized_size(pi);
  metadata_bytes = 17;
  payload_bytes = (cpp->u1len+cpp->u2len+LIFTS)*N*QBYTES;
  s = sqrt(jlproj_normsq(pi->p));
  projection_bits = (log2(s)-4+2.05)*256;
  modeled_bits = projection_bits+8.0*payload_bytes;
  inner_norm = 6*T*SLACK*ldexp(1,(cpp->f-1)*cpp->b)*sqrt(pi->normsq);
  outer_norm = 2*SLACK*sqrt(pi->normsq);

  printf("Labrador fold proof parameters:\n");
  printf("  SIS security policy: %s\n",sis_security_mode_name());
  printf("  Algebra:\n");
  printf("    Ring: Z_q[X]/(X^%d + 1), q = 2^%d - %d\n",N,LOGQ,QOFF);
  printf("    Wire encoding per ring element: %d coefficients x %d bits = %d bytes\n",
         N,LOGQ,N*QBYTES);
  printf("  Input witness layout:\n");
  printf("    Source parts (r): %zu\n",pi->r);
  printf("    Part ranks n[i] (ring elements in each source vector): ");
  for(i=0;i<MIN(pi->r,10);i++) {
    printf("%zu",pi->n[i]);
    if(i+1<MIN(pi->r,10)) printf(", ");
  }
  if(pi->r>10) printf(", ...");
  printf("\n");
  printf("    Fold blocks nu[i] (0 continues a group; nonzero closes it with this many blocks): ");
  for(i=0;i<MIN(pi->r,10);i++) {
    printf("%zu",pi->nu[i]);
    if(i+1<MIN(pi->r,10)) printf(", ");
  }
  if(pi->r>10) printf(", ...");
  printf("\n");
  printf("    Resulting amortized rows (sum nu): %zu\n",groups);
  printf("  Digit decompositions (base B = 2^logB, depth d):\n");
  printf("    Opening witness:       logB = %zu, B = 2^%zu, d = %zu\n",cpp->b,cpp->b,cpp->f);
  printf("    Uniform commitment:    logB = %zu, B = 2^%zu, d = %zu\n",cpp->bu,cpp->bu,cpp->fu);
  printf("    Quadratic cross-terms: logB = %zu, B = 2^%zu, d = %zu%s\n",
         cpp->bg,cpp->bg,cpp->fg,cpp->fg ? "" : " (not present)");
  printf("  SIS commitment dimensions:\n");
  printf("    Inner module rank kappa:  %zu (tested norm bound %.2e; %s)\n",
         cpp->kappa,inner_norm,sis_secure(cpp->kappa,row_width,inner_norm) ? "passes estimator" : "FAILS estimator");
  printf("    Outer module rank kappa1: %zu",cpp->kappa1);
  if(pi->tail)
    printf(" (tail proof: outer commitment is omitted)\n");
  else
    printf(" (tested norm bound %.2e; %s)\n",outer_norm,
           sis_secure(cpp->kappa1,MAX(outer1_width,outer2_width),outer_norm) ? "passes estimator" : "FAILS estimator");
  printf("    These ranks are matrix row/module dimensions, not bit-security estimates.\n");
  printf("  SIS audit instances (scalar dimensions are n = rank*%d, m = width*%d):\n",N,N);
  print_sis_audit_pp("labrador-inner",cpp->kappa,row_width,inner_norm);
  if(!pi->tail) {
    print_sis_audit_pp("labrador-outer-u1",cpp->kappa1,outer1_width,outer_norm);
    print_sis_audit_pp("labrador-outer-u2",cpp->kappa1,outer2_width,outer_norm);
  }
  printf("  Johnson-Lindenstrauss norm proof:\n");
  printf("    Projection dimension: 256 signed coordinates\n");
  printf("    Matrix retry nonce: %zu (number of deterministic candidates tried)\n",pi->jlnonce);
  printf("    Accepted projected l2 norm: %.2f\n",s);
  printf("  Proof payload:\n");
  printf("    u1 ring elements: %zu (%s)\n",cpp->u1len,
         pi->tail ? "embedded inner commitments" : "first outer commitment");
  printf("    u2 ring elements: %zu (%s)\n",cpp->u2len,
         pi->tail ? "embedded opening/cross terms" : "second outer commitment");
  printf("    Integer-to-ring lift elements: %d\n",LIFTS);
  printf("  Contextual proof encoding (schedule shape out of band):\n");
  printf("    Rice parameter + JL nonce + next norm: %zu bytes\n",metadata_bytes);
  printf("    JL coordinates: %zu bytes (lossless; %.2f bytes under legacy entropy model)\n",
         wire_bytes-metadata_bytes-payload_bytes,projection_bits/8);
  printf("    Ring-element payload: %zu bytes\n",payload_bytes);
  printf("    Exact total: %zu bytes\n",wire_bytes);
  printf("  Legacy entropy model (comparison only; not a wire encoding):\n");
  printf("    Modeled JL coordinates: %.2f bits = %.2f bytes\n",
         projection_bits,projection_bits/8);
  printf("    Modeled payload total: %.2f bits = %.2f bytes\n",
         modeled_bits,modeled_bits/8);
  printf("\n");
  return (double)wire_bytes/1024;
}

void print_statement_pp(const statement *st) {
  size_t i;

  printf("Reduced statement (the next fold's public input) ");
  if(st->tail)
    printf("(tail) ");
  for(i=0;i<5;i++)
    printf("%02hhX",st->h[i]);
  printf(":\n");

  printf("  Amortized rows r: %zu (independent folded witness rows)\n",st->r);
  printf("  Opening width n: %zu ring elements per row\n",st->n);
  printf("  Auxiliary width m: %zu ring elements for commitment openings/cross-terms\n",st->m);
  printf("  Quadratic constraint sparsity: %zu nonzero matrix entries\n",st->cnst->a->len);
  printf("  Accepted witness l2-norm bound beta: %.2f\n",sqrt(st->betasq));
  printf("\n");
}

double print_witness_pp(const witness *wt) {
  size_t i,wire_bytes;
  double s;

  printf("Witness:\n");
  printf("  Witness parts: %zu\n",wt->r);
  printf("  Witness part lengths: ");
  for(i=0;i<MIN(wt->r,10);i++) {
    printf("%zu",wt->n[i]);
    if(i<wt->r-1) printf(", ");
  }
  printf("\n");
  printf("  Witness part norms: ");
  for(i=0;i<MIN(wt->r,10);i++) {
    printf("%.2f",sqrt(wt->normsq[i]));
    if(i<wt->r-1) printf(", ");
  }
  printf("\n");

  wire_bytes = witness_contextual_serialized_size(wt);
  s = (double)wire_bytes/1024;
  printf("  Contextual witness serialization: %zu bytes (one Rice parameter per part; shape out of band)\n",
         wire_bytes);
  printf("\n");
  return s;
}

size_t commit_raw(polx *u, poly *t, size_t r, size_t n, const polx s[r][n],
                  size_t off, const comparams *cpp, int tail)
{
  size_t i;
  const size_t l = cpp->fu*cpp->kappa;
  polx tx[l];

  if(tail) {
    for(i=0;i<r;i++)
      polxvec_mul_extension(&u[off+i*l],comkey,s[i],n,cpp->kappa,1);
    return off+r*l;
  }

  for(i=0;i<r;i++) {
    /* inner commitment */
    polxvec_mul_extension(tx,comkey,s[i],n,cpp->kappa,1);
    polxvec_decompose(&t[i*l],tx,cpp->kappa,cpp->fu,cpp->bu);

    /* outer commitment */
    polxvec_frompolyvec(tx,&t[i*l],l);
    off += polxvec_mul_extension(tx,&comkey[off],tx,l,cpp->kappa1,1);
    polxvec_add(u,u,tx,cpp->kappa1);
  }

  return off;
}

size_t qugarbage_raw(polx *u, poly *g, size_t r, size_t n, const polx s[r][n],
                     size_t off, const comparams *cpp, int tail)
{
  size_t i,j,k,l;
  const size_t len = (r*r+r)/2;
  polx buf[MAX(cpp->fg*len,cpp->kappa1)];
  polx *gx = (tail) ? u+off : buf;

  if(!cpp->fg)
    return off;

  polxvec_setzero(gx,len);
  for(l=0;l<n;l++) {
    k = 0;
    for(i=0;i<r;i++)
      for(j=i;j<r;j++)
        polx_mul_add(&gx[k++],&s[i][l],&s[j][l]);
  }

  if(!tail) {
    polxvec_decompose(g,gx,len,cpp->fg,cpp->bg);
    polxvec_frompolyvec(gx,g,cpp->fg*len);
    off += polxvec_mul_extension(gx,&comkey[off],gx,cpp->fg*len,cpp->kappa1,1);
    polxvec_add(u,u,gx,cpp->kappa1);
  }

  return off;
}

void commit(statement *ost, witness *owt, proof *pi, polx sx[ost->r][ost->n], const witness *iwt) {
  const size_t r = iwt->r;
  const size_t *n = iwt->n;
  const comparams *cpp = ost->cpp;

  /* offsets in aux witness vector v */
  const size_t t = 0;
  const size_t g = t + ost->r*cpp->fu*cpp->kappa;

  size_t i,j,k,l;
  __attribute__((aligned(16)))
  uint8_t hashbuf[16+cpp->u1len*N*QBYTES];
  poly *v = (pi->tail) ? NULL : owt->s[cpp->f];

  polxvec_setzero(ost->u1,cpp->kappa1);

  /* transform and commit to input witness */
  j = k = l = 0;
  for(i=0;i<r;i++) {
    polxvec_frompolyvec(&sx[j][k],iwt->s[i],n[i]);
    k += n[i];
    if(pi->nu[i]) {
      polxvec_setzero(&sx[j][k],pi->nu[i]*ost->n - k);
      l = commit_raw(ost->u1,&v[t+j*cpp->fu*cpp->kappa],pi->nu[i],ost->n,&sx[j],l,cpp,pi->tail);
      j += pi->nu[i];
      k = 0;
    }
  }

  qugarbage_raw(ost->u1,&v[g],ost->r,ost->n,sx,l,cpp,pi->tail);
  polzvec_frompolxvec(pi->u1,ost->u1,cpp->u1len);
  memcpy(hashbuf,ost->h,16);
  polzvec_bitpack(&hashbuf[16],pi->u1,cpp->u1len);
  shake128(ost->h,16,hashbuf,sizeof(hashbuf));
}

void reduce_commit(statement *ost, const proof *pi) {
  __attribute__((aligned(16)))
  uint8_t hashbuf[16+pi->cpp->u1len*N*QBYTES];

  /* first outer commitment resp inner comms + quadratic garbage */
  polzvec_topolxvec(ost->u1,pi->u1,pi->cpp->u1len);
  memcpy(hashbuf,ost->h,16);
  polzvec_bitpack(&hashbuf[16],pi->u1,pi->cpp->u1len);
  shake128(ost->h,16,hashbuf,sizeof(hashbuf));
}

static uint64_t next2power(uint64_t a) {
  a -= 1;
  a |= a >>  1;
  a |= a >>  2;
  a |= a >>  4;
  a |= a >>  8;
  a |= a >> 16;
  a |= a >> 32;
  a += 1;
  return a;
}

int project(statement *ost, proof *pi, uint8_t jlmat[][ost->n][256*N/8], const witness *iwt) {
  const size_t r = iwt->r;
  const size_t *n = iwt->n;

  size_t i,j,k,rep;
  uint64_t normsq = 0, test;
  __attribute__((aligned(16)))
  uint8_t hashbuf[16+1024];
  aes128ctr_ctx aesctx;

  for(i=0;i<r;i++)
    normsq += iwt->normsq[i];
  if(normsq > JLMAXNORMSQ) {
    fprintf(stderr,"ERROR: Total witness norm too big for JL Projection\n");
    return 1;
  }
  test = next2power(4*sqrt(normsq));
  normsq *= 256;
  shake128(hashbuf,32,ost->h,16);
  aes128ctr_init(&aesctx,&hashbuf[16],0);
  pi->jlnonce = 0;
  do {
    aes128ctr_select(&aesctx,++pi->jlnonce);
    memset(pi->p,0,256*sizeof(int32_t));
    j = k = 0;
    for(i=0;i<r;i++) {
      aes128ctr_squeezeblocks(jlmat[j][k],n[i]*256*N/8/AES128CTR_BLOCKBYTES,&aesctx);
      polyvec_jlproj_add(pi->p,iwt->s[i],n[i],jlmat[j][k]);
      k += n[i];
      if(pi->nu[i]) {
        memset(jlmat[j][k],0,(pi->nu[i]*ost->n-k)*256*N/8);
        j += pi->nu[i];
        k = 0;
      }
    }
    rep = 0;
    for(i=0;i<256;i++)  // TODO: vectorize
      if((uint64_t)labs(pi->p[i]) >= test)
        rep = 1;
  } while(rep || jlproj_normsq(pi->p) > normsq);

  memcpy(&hashbuf[16],pi->p,1024);
  shake128(ost->h,16,hashbuf,16+1024);
  return 0;
}

int reduce_project(statement *ost, uint8_t jlmat[][ost->n][256*N/8], const proof *pi, size_t r, uint64_t betasq) {
  size_t i,j,k;
  __attribute__((aligned(16)))
  uint8_t hashbuf[16+1024];
  aes128ctr_ctx aesctx;

  if(jlproj_normsq(pi->p) > 256*MIN(JLMAXNORMSQ,betasq)) {
    fprintf(stderr,"ERROR in reduce_project(): Witness projection longer than bound\n");
    return 1;
  }

  shake128(hashbuf,32,ost->h,16);
  aes128ctr_init(&aesctx,&hashbuf[16],pi->jlnonce);
  memcpy(&hashbuf[16],pi->p,1024);
  shake128(ost->h,16,hashbuf,sizeof(hashbuf));

  /* JL Matrix */
  j = k = 0;
  for(i=0;i<r;i++) {
    aes128ctr_squeezeblocks(jlmat[j][k],pi->n[i]*256*N/8/AES128CTR_BLOCKBYTES,&aesctx);
    k += pi->n[i];
    if(pi->nu[i]) {
      memset(jlmat[j][k],0,(pi->nu[i]*ost->n-k)*256*N/8);
      j += pi->nu[i];
      k = 0;
    }
  }

  return 0;
}

void collaps_jlproj_raw(constraint *cnst, size_t r, size_t n, uint8_t h[16], const int32_t p[256],
                       const uint8_t jlmat[r][n][256*N/8])
{
  __attribute__((aligned(32)))
  uint8_t hashbuf[32+QBYTES*256+24];  // additional 24 bytes because of vector loads
  int64_t x;

  shake128(hashbuf,sizeof(hashbuf),h,16);
  memcpy(h,hashbuf,16);
  polxvec_jlproj_collapsmat(cnst->phi,**jlmat,r*n,&hashbuf[32]);
  x = jlproj_collapsproj(p,&hashbuf[32]);
  polx_monomial(cnst->b,x,0);
}

void collaps_jlproj(constraint *cnst, statement *st, const proof *pi, const uint8_t jlmat[st->r][st->n][256*N/8]) {
  collaps_jlproj_raw(cnst,st->r,st->n,st->h,pi->p,jlmat);
}

void lift_aggregate_zqcnst(statement *ost, proof *pi, size_t i, constraint *cnst, const polx sx[ost->r][ost->n]) {
  __attribute__((aligned(16)))
  uint8_t hashbuf[16+N*QBYTES];
  polx alpha[1];

  polxvec_sprod(cnst->b,cnst->phi,*sx,ost->r*ost->n);
  polz_frompolx(&pi->bb[i],cnst->b);
  polz_setcoeff_fromint64(&pi->bb[i],0,0);
  memcpy(hashbuf,ost->h,16);
  polz_bitpack(&hashbuf[16],&pi->bb[i]);
  shake128(hashbuf,32,hashbuf,sizeof(hashbuf));
  memcpy(ost->h,hashbuf,16);
  polxvec_quarternary(alpha,1,&hashbuf[16],0);
  if(!i) {
    polxvec_polx_mul(ost->cnst->phi,alpha,cnst->phi,ost->r*ost->n);
    polx_mul(ost->cnst->b,alpha,cnst->b);
    sparsemat_polx_mul(ost->cnst->a,alpha,cnst->a);
  }
  else {
    polxvec_polx_mul_add(ost->cnst->phi,alpha,cnst->phi,ost->r*ost->n);
    polx_mul_add(ost->cnst->b,alpha,cnst->b);
    sparsemat_polx_mul_add(ost->cnst->a,alpha,cnst->a);
  }
}

void reduce_lift_aggregate_zqcnst(statement *ost, const proof *pi, size_t i, const constraint *cnst) {
  __attribute__((aligned(16)))
  uint8_t hashbuf[16+N*QBYTES];
  polz b[1];
  polx alpha[1];
  zz c0[1];

  polz_frompolx(b,cnst->b);
  polz_getcoeff(c0,b,0);
  polzvec_copy(b,&pi->bb[i],1);
  polz_setcoeff(b,c0,0);
  polz_topolx(cnst->b,b);

  memcpy(hashbuf,ost->h,16);
  polz_bitpack(&hashbuf[16],&pi->bb[i]);
  shake128(hashbuf,32,hashbuf,sizeof(hashbuf));
  memcpy(ost->h,hashbuf,16);
  polxvec_quarternary(alpha,1,&hashbuf[16],0);
  if(!i) {
    polxvec_polx_mul(ost->cnst->phi,alpha,cnst->phi,ost->r*ost->n);
    polx_mul(ost->cnst->b,alpha,cnst->b);
    sparsemat_polx_mul(ost->cnst->a,alpha,cnst->a);
  }
  else {
    polxvec_polx_mul_add(ost->cnst->phi,alpha,cnst->phi,ost->r*ost->n);
    polx_mul_add(ost->cnst->b,alpha,cnst->b);
    sparsemat_polx_mul_add(ost->cnst->a,alpha,cnst->a);
  }
}

static void aggregate(statement *ost, const proof *pi, const statement *ist) {
  const size_t n = ist->n;
  const size_t m = ist->m;
  const size_t r = ist->r;
  const comparams *cpp = ist->cpp;

  const size_t l = cpp->fu*cpp->kappa;
  /* offsets in aux witness vector v */
  const size_t t = 0;
  const size_t g = t + r*l;
  const size_t h = g + cpp->fg*(r*r+r)/2;

  size_t i,j,k;
  polx *tmp = (polx*)_aligned_alloc(64,MAX(n,m)*sizeof(polx));
  polx (*phi)[ost->n] = (polx(*)[ost->n])ost->cnst->phi;
  j = k = 0;
  for(i=0;i<cpp->f;i++) {
    j += pi->nu[i];
    k  = (pi->nu[i]) ? 0 : k+n;
  }
  polx *phiv = &phi[j][k];

  __attribute__((aligned(16)))
  uint8_t hashbuf[32];
  polx chalbuf[2*cpp->kappa1+cpp->kappa+2];
  const polx *alpha = chalbuf;
  const polx *beta = &alpha[cpp->kappa1];
  const polx *gamma = &beta[cpp->kappa1];
  const polx *delta = &gamma[cpp->kappa];

  shake128(hashbuf,32,ost->h,16);
  memcpy(ost->h,hashbuf,16);
  polxvec_quarternary(chalbuf,sizeof(chalbuf)/sizeof(polx),&hashbuf[16],0);

  /* Bv1 = u1 */
  j = 0;
  for(i=0;i<r;i++)
    j += polxvec_collaps_add_extension(&phiv[t+i*l],alpha,&comkey[j],l,cpp->kappa1,1);
  polxvec_collaps_add_extension(&phiv[g],alpha,&comkey[j],h-g,cpp->kappa1,1);
  polxvec_sprod_add(ost->cnst->b,alpha,ist->u1,cpp->kappa1);

  /* Bv2 = u2 */
  polxvec_collaps_add_extension(&phiv[h],beta,comkey,m-h,cpp->kappa1,1);
  polxvec_sprod_add(ost->cnst->b,beta,ist->u2,cpp->kappa1);

  /*    gamma:         - Az      + \sum_i  c_i t_i                 = 0 */
  /* delta[0]: - <z,z>           + \sum_ij c_i c_j g_ij            = 0 */
  /*        1:         - <phi,z> + \sum_ij c_i c_j h_ij            = 0 */
  /* delta[1]:                   + \sum_ij a_ij g_ij + \sum_i h_ii = b */

  /* -Az, -<phi,z> */
  polxvec_copy(tmp,ist->cnst->phi,n);
  polxvec_collaps_add_extension(tmp,gamma,comkey,n,cpp->kappa,1);
  j = k = 0;
  for(i=0;i<cpp->f;i++) {
    if(i) polxvec_scale(tmp,tmp,n,(int64_t)1 << cpp->b);
    polxvec_sub(&phi[j][k],&phi[j][k],tmp,n);
    j += pi->nu[i];
    k = (pi->nu[i]) ? 0 : k+n;
  }

  /* \sum_i c_i t_i */
  for(i=0;i<r;i++) {
    polxvec_polx_mul(tmp,&ist->c[i],gamma,cpp->kappa);
    for(j=0;j<cpp->fu;j++) {
      if(j) polxvec_scale(tmp,tmp,cpp->kappa,(int64_t)1 << cpp->bu);
      polxvec_add(&phiv[t+i*l+j*cpp->kappa],&phiv[t+i*l+j*cpp->kappa],tmp,cpp->kappa);
    }
  }

  /* \sum_ij c_i c_j h_ij */
  k = 0;
  for(i=0;i<r;i++)
    for(j=i;j<r;j++)
      polx_mul(&tmp[k++],&ist->c[i],&ist->c[j]);

  if(ist->cnst->a->len) {
    /* \sum_ij c_i c_j g_ij */
    polx_scale(chalbuf,&delta[0],2);
    j = 0;
    for(i=0;i<r;i++) {
      polx_mul(&tmp[k+j],&delta[0],&tmp[j]);
      polxvec_polx_mul(&tmp[k+j+1],chalbuf,&tmp[j+1],r-i-1);
      j += r-i;
    }

    /* \sum_ij a_ij g_ij */
    for(i=0;i<ist->cnst->a->len;i++) {
      j = triangularidx(ist->cnst->a->rows[i],ist->cnst->a->cols[i],r);
      polx_mul_add(&tmp[k+j],&delta[1],&ist->cnst->a->coeffs[i]);
    }
  }

  /* \sum_i h_ii */
  j = 0;
  for(i=0;i<r;i++) {
    polx_add(&tmp[j],&delta[1],&tmp[j]);
    j += r-i;
  }

  for(i=0;i<cpp->fg;i++) {
    if(i) polxvec_scale(&tmp[k],&tmp[k],k,(int64_t)1 << cpp->bg);
    polxvec_add(&phiv[g+i*k],&phiv[g+i*k],&tmp[k],k);
  }
  for(i=0;i<cpp->fu;i++) {
    if(i) polxvec_scale(&tmp[0],&tmp[0],k,(int64_t)1 << cpp->bu);
    polxvec_add(&phiv[h+i*k],&phiv[h+i*k],&tmp[0],k);
  }

  /* b */
  polx_mul_add(ost->cnst->b,&delta[1],ist->cnst->b);

  /* a */
  /* Recall: ost->r = cpp->f*pi->nu[0] + pi->nu[cpp->f] */
  if(ist->cnst->a->len) {
    for(i=0;i<cpp->f;i++) {
      for(j=0;j<pi->nu[0];j++) {
        for(k=i;k<cpp->f;k++) {
          int64_t s = (k==i) ? -1 : -2;
          s <<= (i+k)*cpp->b;
          polx_scale(&ost->cnst->a->coeffs[ost->cnst->a->len],&delta[0],s);
          ost->cnst->a->rows[ost->cnst->a->len] = i*pi->nu[0]+j;
          ost->cnst->a->cols[ost->cnst->a->len] = k*pi->nu[0]+j;
          ost->cnst->a->len += 1;
        }
      }
    }
    // reducing size; alignment should be preserved
    ost->cnst->a->coeffs = realloc(ost->cnst->a->coeffs,ost->cnst->a->len*sizeof(polx));
  }

  free(tmp);
}

static void amortize_tail(statement *ost, witness *owt, proof *pi, polx sx[ost->r][ost->n]) {
  const size_t r = ost->r;
  const size_t n = ost->n;
  const comparams *cpp = ost->cpp;

  size_t i;
  __attribute__((aligned(16)))
  uint8_t hashbuf[16+2*N*QBYTES];
  polx (*phi)[n] = (polx(*)[n])ost->cnst->phi;
  polx *hx = ost->u2;
  polz *hz = pi->u2;

  polxvec_sprod(&hx[0],phi[0],sx[0],n);
  polz_frompolx(&hz[0],&hx[0]);

  memcpy(hashbuf,ost->h,16);
  polzvec_bitpack(&hashbuf[16],&hz[0],1);
  shake128(hashbuf,32,hashbuf,16+N*QBYTES);
  polxvec_challenge(&ost->c[0],1,&hashbuf[16],0);

  polxvec_polx_mul(sx[0],&ost->c[0],sx[0],n);
  polxvec_polx_mul(phi[0],&ost->c[0],phi[0],n);

  for(i=1;i<r;i++) {
    polxvec_sprod(&hx[2*i-1],phi[i],sx[0],n);
    polxvec_sprod_add(&hx[2*i-1],phi[0],sx[i],n);
    polxvec_sprod(&hx[2*i],phi[i],sx[i],n);
    polzvec_frompolxvec(&hz[2*i-1],&hx[2*i-1],2);
    polzvec_bitpack(&hashbuf[16],&hz[2*i-1],2);
    shake128(hashbuf,32,hashbuf,16+2*N*QBYTES);
    polxvec_challenge(&ost->c[i],1,&hashbuf[16],0);
    polxvec_polx_mul_add(sx[0],&ost->c[i],sx[i],n);
    polxvec_polx_mul_add(phi[0],&ost->c[i],phi[i],n);
  }

  memcpy(ost->h,hashbuf,16);
  polxvec_decompose(owt->s[0],sx[0],n,cpp->f,cpp->b);

  ost->betasq = 0;
  for(i=0;i<cpp->f;i++) {
    owt->normsq[i] = polyvec_sprodz(owt->s[i],owt->s[i],n);
    ost->betasq += owt->normsq[i];
  }
  pi->normsq = ost->betasq;

  ost->cnst->phi = realloc(ost->cnst->phi,n*sizeof(polx));
}

void amortize(statement *ost, witness *owt, proof *pi, polx sx[ost->r][ost->n]) {
  const size_t r = ost->r;
  const size_t n = ost->n;
  const size_t m = ost->m;
  const comparams *cpp = ost->cpp;

  const size_t t = 0;
  const size_t g = t + r*cpp->fu*cpp->kappa;
  const size_t h = g + cpp->fg*(r*r+r)/2;

  if(pi->tail) {
    amortize_tail(ost,owt,pi,sx);
    return;
  }

  size_t i,j,k,l;
  poly *vh = (poly*)&owt->s[cpp->f][h];
  polx (*phi)[n] = (polx(*)[n])ost->cnst->phi;
  polx *hx = (polx*)_aligned_alloc(64,(m-h)*sizeof(polx));

  /* linear garbage */
  polxvec_setzero(hx,(r*r+r)/2);
  for(l=0;l<n;l++) {
    k = 0;
    for(i=0;i<r;i++) {
      polx_mul_add(&hx[k++],&phi[i][l],&sx[i][l]);
      for(j=i+1;j<r;j++) {
        polx_mul_add(&hx[k],&phi[i][l],&sx[j][l]);
        polx_mul_add(&hx[k++],&phi[j][l],&sx[i][l]);
      }
    }
  }
  polxvec_decompose(vh,hx,(r*r+r)/2,cpp->fu,cpp->bu);
  polxvec_frompolyvec(hx,vh,m-h);

  /* second outer commitment */
  polxvec_mul_extension(ost->u2,comkey,hx,m-h,cpp->kappa1,1);
  polzvec_frompolxvec(pi->u2,ost->u2,cpp->u2len);

  /* amortization */
  __attribute__((aligned(16)))
  uint8_t hashbuf[16+cpp->kappa1*N*QBYTES];
  memcpy(hashbuf,ost->h,16);
  polzvec_bitpack(&hashbuf[16],pi->u2,cpp->kappa1);
  shake128(hashbuf,32,hashbuf,sizeof(hashbuf));
  memcpy(ost->h,hashbuf,16);
  polxvec_challenge(ost->c,r,&hashbuf[16],0);

  polxvec_polx_mul(sx[0],&ost->c[0],sx[0],n);
  polxvec_polx_mul(phi[0],&ost->c[0],phi[0],n);
  for(i=1;i<r;i++) {
    polxvec_polx_mul_add(sx[0],&ost->c[i],sx[i],n);
    polxvec_polx_mul_add(phi[0],&ost->c[i],phi[i],n);
  }
  polxvec_decompose(owt->s[0],sx[0],n,cpp->f,cpp->b);

  ost->betasq = 0;
  for(i=0;i<cpp->f;i++) {
    owt->normsq[i] = polyvec_sprodz(owt->s[i],owt->s[i],n);
    ost->betasq += owt->normsq[i];
  }
  owt->normsq[cpp->f] = polyvec_sprodz(owt->s[cpp->f],owt->s[cpp->f],m);
  ost->betasq += owt->normsq[cpp->f];
  pi->normsq = ost->betasq;

  ost->cnst->phi = realloc(ost->cnst->phi,n*sizeof(polx));
  free(hx);
}

int reduce_amortize(statement *ost, const proof *pi) {
  const size_t r = ost->r;
  const size_t n = ost->n;
  const comparams *cpp = ost->cpp;

  size_t i;
  polx (*phi)[n] = (polx(*)[n])ost->cnst->phi;

  ost->betasq = pi->normsq;
  if(!sis_secure(cpp->kappa,ost->n,
                 6*T*SLACK*ldexp(1,(cpp->f-1)*cpp->b)*sqrt(ost->betasq))) {
    fprintf(stderr,"ERROR in reduce_amortize(): Inner commitments not secure\n");
    return 1;
  }
  if(!pi->tail &&
     (!sis_secure(cpp->kappa1,cpp->fu*ost->r*cpp->kappa+
                              cpp->fg*(ost->r*ost->r+ost->r)/2,
                  2*SLACK*sqrt(ost->betasq)) ||
      !sis_secure(cpp->kappa1,cpp->fu*(ost->r*ost->r+ost->r)/2,
                  2*SLACK*sqrt(ost->betasq)))) {
    fprintf(stderr,"ERROR in reduce_amortize(): Outer commitments not secure\n");
    return 2;
  }

  /* second outer commitment resp garbage terms */
  polzvec_topolxvec(ost->u2,pi->u2,cpp->u2len);

  __attribute__((aligned(16)))
  uint8_t hashbuf[16+cpp->u2len*N*QBYTES];
  memcpy(hashbuf,ost->h,16);

  if(pi->tail) {
    polzvec_bitpack(&hashbuf[16],pi->u2,1);
    shake128(hashbuf,32,hashbuf,16+N*QBYTES);
    polxvec_challenge(&ost->c[0],1,&hashbuf[16],0);
    for(i=1;i<r;i++) {
      polzvec_bitpack(&hashbuf[16],&pi->u2[2*i-1],2);
      shake128(hashbuf,32,hashbuf,16+2*N*QBYTES);
      polxvec_challenge(&ost->c[i],1,&hashbuf[16],0);
    }
  }
  else {
    polzvec_bitpack(&hashbuf[16],pi->u2,cpp->u2len);
    shake128(hashbuf,32,hashbuf,sizeof(hashbuf));
    polxvec_challenge(ost->c,r,&hashbuf[16],0);
  }
  memcpy(ost->h,hashbuf,16);

  polxvec_polx_mul(phi[0],&ost->c[0],phi[0],n);
  for(i=1;i<r;i++)
    polxvec_polx_mul_add(phi[0],&ost->c[i],phi[i],n);

  ost->cnst->phi = realloc(ost->cnst->phi,n*sizeof(polx));
  return 0;
}

int prove(statement *ost, witness *owt, proof *pi, const statement *ist, const witness *iwt, int tail) {
  int ret;
  size_t i;
  constraint cnst[1] = {};
  void *buf = NULL;

  ret = init_proof(pi,iwt,ist->cpp->fg != 0,tail);
  if(ret) // commitments not secure (1/2)
    return ret;
  init_statement(ost,pi,ist->h);
  init_witness(owt,ost);
  printf("Predicted witness norm: %.2f\n\n",sqrt(pi->normsq));

  {
    buf = _aligned_alloc(64,ost->r*ost->n*(sizeof(polx)+256*N/8));
    polx (*sx)[ost->n] = (polx(*)[ost->n])buf;
    uint8_t (*jlmat)[ost->n][256*N/8] = (uint8_t(*)[ost->n][256*N/8])sx[ost->r];
    commit(ost,owt,pi,sx,iwt);
    ret = project(ost,pi,jlmat,iwt);
    if(ret) {
      ret += 10;
      goto err;
    }

    init_constraint_raw(cnst,ost->r,ost->n,1,0);
    for(i=0;i<LIFTS;i++) {
      collaps_jlproj(cnst,ost,pi,jlmat);
      lift_aggregate_zqcnst(ost,pi,i,cnst,sx);
    }
    free_constraint(cnst);

    aggregate(ost,pi,ist);
    amortize(ost,owt,pi,sx);
    free(buf);
    buf = NULL;
  }

  polx_refresh(ost->cnst->b);
  polxvec_refresh(ost->cnst->phi,ost->n);
  return 0;

err:
  free_proof(pi);
  free_statement(ost);
  free_witness(owt);
  free(buf);
  free_constraint(cnst);
  return ret;
}

int reduce(statement *ost, const proof *pi, const statement *ist) {
  size_t i;
  int ret;
  uint8_t (*jlmat)[ost->n][256*N/8];
  constraint cnst[1];

  init_statement(ost,pi,ist->h);
  init_constraint(cnst,ost);
  jlmat = _aligned_alloc(64,ost->r*ost->n*256*N/8);

  reduce_commit(ost,pi);
  ret = reduce_project(ost,jlmat,pi,pi->r,ist->betasq);
  if(ret) goto err;  // projection too long

  for(i=0;i<LIFTS;i++) {
    collaps_jlproj(cnst,ost,pi,jlmat);
    reduce_lift_aggregate_zqcnst(ost,pi,i,cnst);
  }
  free_constraint(cnst);
  free(jlmat);
  jlmat = NULL;

  aggregate(ost,pi,ist);
  ret = reduce_amortize(ost,pi);
  if(ret) {  // commitments not secure (1/2)
    ret += 10;
    goto err;
  }

  polx_refresh(ost->cnst->b);
  polxvec_refresh(ost->cnst->phi,ost->n);
  return 0;

err:
  free_statement(ost);
  free_constraint(cnst);
  free(jlmat);
  return ret;
}

int verify(const statement *st, const witness *wt) {
  const comparams *cpp = st->cpp;
  const size_t r = st->r;
  const size_t n = st->n;
  const size_t m = st->m;

  const size_t l = cpp->fu*cpp->kappa;
  const size_t t = 0;
  const size_t g = t + r*l;
  const size_t h = g + cpp->fg*(r*r+r)/2;

  size_t i,j,k;
  int ret = 0;
  uint64_t normsq = 0;
  polx tmp0[cpp->kappa];
  polx tmp1[cpp->kappa1];

  void *buf = _aligned_alloc(64,(n+m)*sizeof(polx));
  polx *z = (polx*)buf;
  polx *v = &z[n];

  for(i=0;i<wt->r;i++)
    normsq += polyvec_sprodz(wt->s[i],wt->s[i],wt->n[i]);
  if(normsq > st->betasq) {
    fprintf(stderr,"ERROR in verify(): Total witness norm bigger than bound\n");
    ret = 1;
    goto end;
  }

  if(!st->tail) {
    /* Bv1 = u1 */
    polxvec_copy(tmp0,st->u1,cpp->kappa1);
    j = 0;
    for(i=0;i<r;i++) {
      polxvec_frompolyvec(&v[t+i*l],&wt->s[cpp->f][t+i*l],l);
      j += polxvec_mul_extension(tmp1,&comkey[j],&v[t+i*l],l,cpp->kappa1,1);
      polxvec_sub(tmp0,tmp0,tmp1,cpp->kappa1);
    }
    polxvec_frompolyvec(&v[g],&wt->s[cpp->f][g],h-g);
    polxvec_mul_extension(tmp1,&comkey[j],&v[g],h-g,cpp->kappa1,1);
    polxvec_sub(tmp0,tmp0,tmp1,cpp->kappa1);
    if(!polxvec_iszero(tmp0,cpp->kappa1)) {
      fprintf(stderr,"ERROR in verify(): First outer commitment opening wrong\n");
      ret = 2;
      goto end;
    }

    /* Bv2 = u2 */
    polxvec_frompolyvec(&v[h],&wt->s[cpp->f][h],m-h);
    polxvec_mul_extension(tmp0,comkey,&v[h],m-h,cpp->kappa1,1);
    polxvec_sub(tmp0,tmp0,st->u2,cpp->kappa1);
    if(!polxvec_iszero(tmp0,cpp->kappa1)) {
      fprintf(stderr,"ERROR in verify(): Second outer commitment opening wrong\n");
      ret = 3;
      goto end;
    }

    /* reconstruct inner coms and garbage */
    for(i=0;i<r;i++)
      for(j=1;j<cpp->fu;j++)
        polxvec_scale_add(&v[t+i*l],&v[t+i*l+j*cpp->kappa],cpp->kappa,(int64_t)1 << j*cpp->bu);
    for(i=1;i<cpp->fg;i++)
      polxvec_scale_add(&v[g],&v[g+i*(r*r+r)/2],(r*r+r)/2,(int64_t)1 << i*cpp->bg);
    for(i=1;i<cpp->fu;i++)
      polxvec_scale_add(&v[h],&v[h+i*(r*r+r)/2],(r*r+r)/2,(int64_t)1 << i*cpp->bu);
  }
  else {
    polxvec_copy(&v[t],st->u1,cpp->u1len);
    polxvec_copy(&v[h],st->u2,cpp->u2len);
  }

  /* reconstruct z */
  polxvec_reconstruct(z,wt->s[0],n,cpp->f,cpp->b);

  /* Az = \sum_i c_i t_i */
  polxvec_polx_mul(&v[t],&st->c[0],&v[t],cpp->kappa);
  for(i=1;i<r;i++)
    polxvec_polx_mul_add(&v[t],&st->c[i],&v[t+i*l],cpp->kappa);
  polxvec_mul_extension(tmp0,comkey,z,n,cpp->kappa,1);
  polxvec_sub(tmp0,tmp0,&v[t],cpp->kappa);
  if(!polxvec_iszero(tmp0,cpp->kappa)) {
    fprintf(stderr,"ERROR in verify(): Amortized (inner commitment) opening wrong\n");
    ret = 4;
    goto end;
  }

  /* \sum_ij a_ij g_ij + \sum_i h_ii = b */
  polx_neg(tmp0,st->cnst->b);
  j = 0;
  for(i=0;i<r;i++) {
    polx_add(tmp0,tmp0,&v[h+j]);
    j += (st->tail) ? 2 : r-i;
  }
  for(i=0;i<st->cnst->a->len;i++) {
    j = triangularidx(st->cnst->a->rows[i],st->cnst->a->cols[i],r);
    polx_mul_add(tmp0,&st->cnst->a->coeffs[i],&v[g+j]);
  }
  if(!polxvec_iszero(tmp0,1)) {
    fprintf(stderr,"ERROR in verify(): Aggregated dot-product constraint doesn't hold\n");
    ret = 5;
    goto end;
  }

  /* <z,z> = \sum_ij c_i c_j g_ij */
  if(st->cnst->a->len) {
    k = 0;
    for(i=0;i<r;i++) {
      polx_mul(tmp0,&st->c[i],&v[g+k]);
      polxvec_sprod(&tmp0[1],&st->c[i+1],&v[g+k+1],r-1-i);
      polx_scale_add(tmp0,&tmp0[1],2);
      if(i) polx_mul_add(&v[g],&st->c[i],tmp0);
      else polx_mul(&v[g],&st->c[i],tmp0);
      k += r-i;
    }
    polxvec_sprod(tmp0,z,z,n);
    polx_sub(tmp0,tmp0,&v[g]);
    if(!polx_iszero(tmp0)) {
      fprintf(stderr,"ERROR in verify(): Quadratic garbage polynomials wrong\n");
      ret = 6;
      goto end;
    }
  }

  /* <phi,z> = \sum_ij c_i c_j h_ij */
  if(!st->tail) {
    k = 0;
    for(i=0;i<r;i++) {
      polxvec_sprod(tmp0,&st->c[i],&v[h+k],r-i);
      if(i) polx_mul_add(&v[h],&st->c[i],tmp0);
      else polx_mul(&v[h],&st->c[i],tmp0);
      k += r-i;
    }
  }
  else {
    polx_mul(&v[h],&st->c[0],&v[h]);
    polx_mul(&v[h],&st->c[0],&v[h]);
    for(i=1;i<r;i++) {
      polx_mul_add(&v[h+2*i-1],&st->c[i],&v[h+2*i]);
      polx_mul_add(&v[h],&st->c[i],&v[h+2*i-1]);
    }
  }
  polxvec_sprod(tmp0,st->cnst->phi,z,n);
  polx_sub(tmp0,tmp0,&v[h]);
  if(!polxvec_iszero(tmp0,1)) {
    fprintf(stderr,"ERROR in verify(): Linear garbage polynomials wrong\n");
    ret = 7;
    goto end;
  }

end:
  free(buf);
  return ret;
}
