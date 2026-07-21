#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "malloc.h"
#include "labrador.h"
#include "chihuahua.h"
#include "dachshund.h"
#include "greyhound.h"
#include "pack.h"

static double wall_time(void) {
  struct timespec t;

  clock_gettime(CLOCK_MONOTONIC,&t);
  return (double)t.tv_sec+1e-9*t.tv_nsec;
}

static void print_pack_size_pp(const char *label, const composite *p) {
  double total_bytes = p->size*1024;
  double tail_bytes = p->tail_size*1024;

  printf("%s proof components: total = %.0f bytes, fold = %.0f bytes, tail = %.0f bytes\n",
         label,total_bytes,total_bytes-tail_bytes,tail_bytes);
}

void free_composite(composite *p) {
  size_t i;

  free_witness(&p->owt);
  for(i=0;i<p->l;i++) {
    free_proof(p->pi[i]);
    free(p->pi[i]);
    p->pi[i] = NULL;
  }
  p->l = 0;
}

static int composite_prove(composite *p, statement *tst, witness *twt, double *twtsize) {
  int ret;
  size_t i = 0;
  double pisize;

  while(p->l < 16) {
    p->pi[p->l] = _malloc(sizeof(proof));
    ret = prove(&tst[i^1],&twt[i^1],p->pi[p->l],&tst[i],&twt[i],0);
    if(ret) return ret;
    printf("Labrador fold %zu candidate:\n\n",p->l+1);
    pisize = print_proof_pp(p->pi[p->l]);
    print_statement_pp(&tst[i^1]);
    twtsize[i^1] = print_witness_pp(&twt[i^1]);
    if(pisize + twtsize[i^1] >= twtsize[i]) {
      printf("Labrador fold %zu rejected: %.2f + %.2f >= %.2f bytes\n\n",
             p->l+1,pisize*1024,twtsize[i^1]*1024,twtsize[i]*1024);
      free_proof(p->pi[p->l]);
      free_statement(&tst[i^1]);
      free_witness(&twt[i^1]);
      break;
    }

    printf("Labrador fold %zu accepted: %.2f + %.2f < %.2f bytes\n\n",
           p->l+1,pisize*1024,twtsize[i^1]*1024,twtsize[i]*1024);

    free_statement(&tst[i]);
    free_witness(&twt[i]);
    p->size += pisize;
    p->l += 1;
    i ^= 1;
  }

  if(p->l < 16) {
    ret = prove(&tst[i^1],&twt[i^1],p->pi[p->l],&tst[i],&twt[i],1);
    if(ret) return ret;
    printf("Labrador fold %zu tail candidate:\n\n",p->l+1);
    pisize = print_proof_pp(p->pi[p->l]);
    print_statement_pp(&tst[i^1]);
    twtsize[i^1] = print_witness_pp(&twt[i^1]);
    if(pisize + twtsize[i^1] >= twtsize[i]) {
      printf("Labrador fold %zu tail rejected: %.2f + %.2f >= %.2f bytes\n\n",
             p->l+1,pisize*1024,twtsize[i^1]*1024,twtsize[i]*1024);
      free_proof(p->pi[p->l]);
      free(p->pi[p->l]);
    }
    else {
      printf("Labrador fold %zu tail accepted: %.2f + %.2f < %.2f bytes\n\n",
             p->l+1,pisize*1024,twtsize[i^1]*1024,twtsize[i]*1024);
      p->size += pisize;
      p->l += 1;
      i ^= 1;
    }

    free_statement(&tst[i^1]);
    free_witness(&twt[i^1]);
  }

  free_statement(&tst[i]);
  p->owt = twt[i];
  p->tail_size = twtsize[i];
  p->size += twtsize[i];
  return 0;
}

int composite_prove_principle(composite *p, const prncplstmnt *st, const witness *wt) {
  int ret;
  statement tst[2] = {};
  witness twt[2] = {};
  double twtsize[2];
  double t;

  p->l = 0;
  p->tail_size = 0;
  memset(&p->owt,0,sizeof(witness));
  p->pi[p->l] = _malloc(sizeof(proof));

  t = wall_time();
  ret = principle_prove(tst,twt,p->pi[p->l],st,wt,0);
  if(ret)
    goto err;
  printf("Labrador fold 1 candidate:\n\n");
  p->size = print_proof_pp(p->pi[p->l]);
  print_statement_pp(tst);
  twtsize[0] = print_witness_pp(twt);
  printf("Labrador fold 1 accepted as initial reduction\n\n");
  p->l += 1;
  ret = composite_prove(p,tst,twt,twtsize);
  if(ret) {
    ret += 10;
    goto err;
  }
  t = wall_time()-t;

  printf("Commitment key length: %zu\n",comkey_len);
  printf("Chihuahua Pack members: %zu\n",p->l);
  print_pack_size_pp("Chihuahua",p);
  printf("Chihuahua Pack proving wall time: %.4fs\n",t);
  return 0;

err:
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  free_witness(&twt[0]);
  free_witness(&twt[1]);
  free(p->pi[p->l]);
  p->pi[p->l] = NULL;
  free_composite(p);
  return ret;
}

int composite_prove_simple(composite *p, commitment *com, const smplstmnt *st, const witness *wt) {
  int ret;
  statement tst[2] = {};
  witness twt[2] = {};
  double twtsize[2];
  double t;

  p->l = 0;
  p->tail_size = 0;
  memset(&p->owt,0,sizeof(witness));
  p->pi[p->l] = _malloc(sizeof(proof));

  t = wall_time();
  ret = simple_prove(tst,twt,p->pi[p->l],com,st,wt,0);
  if(ret)
    goto err;
  printf("Labrador fold 1 candidate:\n\n");
  p->size = print_proof_pp(p->pi[p->l]);
  print_statement_pp(tst);
  twtsize[0] = print_witness_pp(twt);
  printf("Labrador fold 1 accepted as initial reduction\n\n");
  p->l += 1;

  ret = composite_prove(p,tst,twt,twtsize);
  if(ret) {
    ret += 10;
    goto err;
  }
  t = wall_time()-t;

  printf("Commitment key length: %zu\n",comkey_len);
  printf("Dachshund Pack members: %zu\n",p->l);
  print_pack_size_pp("Dachshund",p);
  printf("Dachshund Pack proving wall time: %.4fs\n",t);
  return 0;

err:
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  free_witness(&twt[0]);
  free_witness(&twt[1]);
  free(p->pi[p->l]);
  p->pi[p->l] = NULL;
  free_composite(p);
  return ret;
}

int composite_prove_polcom(composite *p, polcomprf *ppi, polcomctx *ctx, uint32_t x, uint32_t y) {
  int ret;
  size_t wire_bytes;
  prncplstmnt tst0[1] = {};
  statement tst[2] = {};
  witness twt[2] = {};
  double twtsize[2];
  double t;

  p->l = 0;
  p->tail_size = 0;
  memset(&p->owt,0,sizeof(witness));
  p->pi[p->l] = _malloc(sizeof(proof));

  t = wall_time();
  ret = polcom_eval(&twt[1],ppi,ctx,x,y);
  if(ret)
    goto err;
  ret = polcom_reduce(tst0,ppi);
  if(ret)
    goto err;
  p->size = print_polcomprf_pp(ppi);
  print_prncplstmnt_pp(tst0);
  twtsize[1] = print_witness_pp(&twt[1]);

  ret = principle_prove(tst,twt,p->pi[p->l],tst0,&twt[1],0);
  if(ret) {
    ret += 10;
    goto err;
  }
  printf("Labrador fold 1 candidate:\n\n");
  p->size += print_proof_pp(p->pi[p->l]);
  print_statement_pp(tst);
  twtsize[0] = print_witness_pp(twt);
  printf("Labrador fold 1 accepted as initial reduction\n\n");
  free_prncplstmnt(tst0);
  free_witness(&twt[1]);
  p->l += 1;

  ret = composite_prove(p,tst,twt,twtsize);
  if(ret) {
    ret += 20;
    goto err;
  }
  t = wall_time()-t;

  wire_bytes = greyhound_pack_contextual_serialized_size(ppi,p);
  if(!wire_bytes) { ret = 30; goto err; }
  p->tail_size = (double)witness_contextual_serialized_size(&p->owt)/1024;
  p->size = (double)wire_bytes/1024;

  printf("Commitment key length: %zu\n",comkey_len);
  printf("Greyhound Pack members: %zu\n",p->l);
  print_pack_size_pp("Greyhound",p);
  printf("Greyhound Pack proving wall time: %.4fs\n",t);
  return 0;

err:
  free_prncplstmnt(tst0);
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  free_witness(&twt[0]);
  free_witness(&twt[1]);
  free(p->pi[p->l]);
  p->pi[p->l] = NULL;
  free_composite(p);
  return ret;
}

static int composite_verify(const composite *p, statement *tst) {
  size_t i,j;
  int ret;

  i = 0;
  for(j=1;j<p->l;j++) {
    ret = reduce(&tst[i^1],p->pi[j],&tst[i]);
    free_statement(&tst[i]);
    i ^= 1;
    if(ret)  // projection too long or commitments not secure (1/2/3)
      return ret + 4*j;
  }

  ret = verify(&tst[i],&p->owt);
  if(ret)
    return ret + 100;

  return 0;
}

int composite_verify_principle(const composite *p, const prncplstmnt *st) {
  int ret;
  statement tst[2] = {};
  double t;

  t = wall_time();
  ret = principle_reduce(tst,p->pi[0],st);
  if(ret)  // projection too long or commitments not secure (1/2/3)
    goto err;

  ret = composite_verify(p,tst);
  if(ret) {
    ret += 10;
    goto err;
  }
  t = wall_time()-t;
  printf("Chihuahua Pack verification wall time: %.4fs\n",t);
  return 0;

err:
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  return ret;
}

int composite_verify_simple(const composite *p, const commitment *com, const smplstmnt *st) {
  int ret;
  statement tst[2] = {};
  double t;

  t = wall_time();
  ret = simple_reduce(tst,p->pi[0],com,st);
  if(ret)
    goto err;

  ret = composite_verify(p,tst);
  if(ret) {
    ret += 10;
    goto err;
  }
  t = wall_time()-t;
  printf("Dachshund Pack verification wall time: %.4fs\n",t);
  return 0;

err:
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  return ret;
}

int composite_verify_polcom(const composite *p, const polcomprf *ppi) {
  int ret;
  prncplstmnt tst0[1] = {};
  statement tst[2] = {};
  double t;

  t = wall_time();
  ret = polcom_reduce(tst0,ppi);
  if(ret)
    goto err;

  ret = principle_reduce(tst,p->pi[0],tst0);
  free_prncplstmnt(tst0);
  if(ret) {
    ret += 10;
    goto err;
  }

  ret = composite_verify(p,tst);
  if(ret) {
    ret += 20;
    goto err;
  }
  t = wall_time()-t;
  printf("Greyhound Pack verification wall time: %.4fs\n",t);
  return 0;

err:
  free_prncplstmnt(tst0);
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  return ret;
}
