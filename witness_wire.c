#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "labrador.h"
#include "rice.h"

#define WITNESS_WIRE_HEADER_BYTES 16
#define WITNESS_WIRE_PART_BYTES 24
#define WITNESS_WIRE_VERSION 1
#define WITNESS_MAX_RICE_K 15

static void put16(uint8_t *out, uint16_t x) {
  out[0] = (uint8_t)x; out[1] = (uint8_t)(x >> 8);
}

static void put64(uint8_t *out, uint64_t x) {
  size_t i; for(i=0;i<8;i++) out[i] = (uint8_t)(x >> (8*i));
}

static uint16_t get16(const uint8_t *in) {
  return (uint16_t)in[0] | (uint16_t)in[1] << 8;
}

static uint64_t get64(const uint8_t *in) {
  size_t i; uint64_t x = 0;
  for(i=0;i<8;i++) x |= (uint64_t)in[i] << (8*i);
  return x;
}

static int part_values(int32_t *values, const poly *part, size_t n) {
  size_t i,j;
  if(n > SIZE_MAX/N) return 1;
  for(i=0;i<n;i++)
    for(j=0;j<N;j++) values[i*N+j] = part[i].vec->c[j];
  return 0;
}

size_t witness_serialized_size(const witness *wt) {
  size_t i,total,count,bytes;
  uint64_t bits;
  int32_t *values;

  if(!wt || !wt->r || !wt->n || !wt->s) return 0;
  if(wt->r > (SIZE_MAX-WITNESS_WIRE_HEADER_BYTES)/WITNESS_WIRE_PART_BYTES) return 0;
  total = WITNESS_WIRE_HEADER_BYTES+WITNESS_WIRE_PART_BYTES*wt->r;
  for(i=0;i<wt->r;i++) {
    if(!wt->n[i] || wt->n[i] > SIZE_MAX/N) return 0;
    count = wt->n[i]*N;
    values = malloc(count*sizeof(*values));
    if(!values || part_values(values,wt->s[i],wt->n[i])) { free(values); return 0; }
    rice_i32_parameter(values,count,WITNESS_MAX_RICE_K,&bits);
    free(values);
    bytes = (size_t)((bits+7)/8);
    if(bytes > SIZE_MAX-total) return 0;
    total += bytes;
  }
  return total;
}

int witness_serialize(uint8_t *out, size_t outlen, const witness *wt) {
  size_t i,j,off,count,bytes;
  uint64_t bits;
  unsigned k;
  int32_t *values;

  if(!out || outlen != witness_serialized_size(wt)) return 1;
  memcpy(out,"LBTW",4);
  put16(&out[4],WITNESS_WIRE_VERSION);
  put16(&out[6],N);
  put64(&out[8],wt->r);
  off = WITNESS_WIRE_HEADER_BYTES+WITNESS_WIRE_PART_BYTES*wt->r;
  for(i=0;i<wt->r;i++) {
    count = wt->n[i]*N;
    values = malloc(count*sizeof(*values));
    if(!values || part_values(values,wt->s[i],wt->n[i])) { free(values); return 2; }
    k = rice_i32_parameter(values,count,WITNESS_MAX_RICE_K,&bits);
    bytes = (size_t)((bits+7)/8);
    j = WITNESS_WIRE_HEADER_BYTES+i*WITNESS_WIRE_PART_BYTES;
    put64(&out[j],wt->n[i]);
    put64(&out[j+8],bytes);
    out[j+16] = (uint8_t)k;
    memset(&out[j+17],0,7);
    if(rice_i32_encode(&out[off],bytes,values,count,k)) { free(values); return 3; }
    free(values);
    off += bytes;
  }
  return off == outlen ? 0 : 4;
}

int witness_deserialize(witness *wt, const uint8_t *in, size_t inlen) {
  size_t i,j,off,r,count,bytes,total,sum_n = 0;
  uint64_t r64,n64,b64,bits;
  unsigned k,best_k;
  size_t *shape = NULL;
  int32_t *values = NULL;

  if(!wt || !in || inlen < WITNESS_WIRE_HEADER_BYTES) return 1;
  if(memcmp(in,"LBTW",4) || get16(&in[4]) != WITNESS_WIRE_VERSION ||
     get16(&in[6]) != N) return 2;
  r64 = get64(&in[8]);
  if(!r64 || r64 > SIZE_MAX || r64 > (inlen-WITNESS_WIRE_HEADER_BYTES)/WITNESS_WIRE_PART_BYTES) return 3;
  r = (size_t)r64;
  shape = malloc(r*sizeof(*shape));
  if(!shape) return 4;
  off = WITNESS_WIRE_HEADER_BYTES+WITNESS_WIRE_PART_BYTES*r;
  total = off;
  for(i=0;i<r;i++) {
    j = WITNESS_WIRE_HEADER_BYTES+i*WITNESS_WIRE_PART_BYTES;
    n64 = get64(&in[j]); b64 = get64(&in[j+8]); k = in[j+16];
    if(!n64 || n64 > SIZE_MAX/N || b64 > SIZE_MAX || k > WITNESS_MAX_RICE_K) goto err;
    for(bytes=17;bytes<24;bytes++) if(in[j+bytes]) goto err;
    shape[i] = (size_t)n64;
    count = shape[i]*N;
    /* Every Rice-coded coefficient needs at least its unary terminator.  This
     * prevents a tiny malformed input from claiming an enormous allocation. */
    if((size_t)b64 < count/8+(count%8 != 0) || shape[i] > SIZE_MAX-sum_n) goto err;
    sum_n += shape[i];
    if((size_t)b64 > inlen-total) goto err;
    total += (size_t)b64;
  }
  if(total != inlen || sum_n > SIZE_MAX/sizeof(poly)) goto err;
  memset(wt,0,sizeof(*wt));
  init_witness_raw(wt,r,shape);
  off = WITNESS_WIRE_HEADER_BYTES+WITNESS_WIRE_PART_BYTES*r;
  for(i=0;i<r;i++) {
    j = WITNESS_WIRE_HEADER_BYTES+i*WITNESS_WIRE_PART_BYTES;
    bytes = (size_t)get64(&in[j+8]); k = in[j+16]; count = shape[i]*N;
    values = malloc(count*sizeof(*values));
    if(!values || rice_i32_decode(values,count,&in[off],bytes,k,UINT16_MAX)) goto decode_err;
    best_k = rice_i32_parameter(values,count,WITNESS_MAX_RICE_K,&bits);
    if(best_k != k || (bits+7)/8 != bytes) goto decode_err;
    for(j=0;j<shape[i];j++)
      for(total=0;total<N;total++) wt->s[i][j].vec->c[total] = (int16_t)values[j*N+total];
    free(values); values = NULL;
    wt->normsq[i] = polyvec_sprodz(wt->s[i],wt->s[i],wt->n[i]);
    off += bytes;
  }
  free(shape);
  return 0;

decode_err:
  free(values);
  free_witness(wt);
err:
  free(shape);
  return 5;
}

size_t witness_contextual_serialized_size(const witness *wt) {
  size_t i,total,count,bytes;
  uint64_t bits;
  int32_t *values;

  if(!wt || !wt->r || !wt->n || !wt->s) return 0;
  total=wt->r;
  for(i=0;i<wt->r;i++) {
    if(!wt->n[i] || wt->n[i]>SIZE_MAX/N) return 0;
    count=wt->n[i]*N; values=malloc(count*sizeof(*values));
    if(!values || part_values(values,wt->s[i],wt->n[i])) { free(values); return 0; }
    rice_i32_parameter(values,count,WITNESS_MAX_RICE_K,&bits); free(values);
    bytes=(size_t)((bits+7)/8);
    if(bytes>SIZE_MAX-total) return 0;
    total+=bytes;
  }
  return total;
}

int witness_serialize_contextual(uint8_t *out, size_t outlen,
                                 const witness *wt) {
  size_t i,off=0,count,bytes;
  uint64_t bits;
  unsigned k;
  int32_t *values;

  if(!out || outlen!=witness_contextual_serialized_size(wt)) return 1;
  for(i=0;i<wt->r;i++) {
    count=wt->n[i]*N; values=malloc(count*sizeof(*values));
    if(!values || part_values(values,wt->s[i],wt->n[i])) { free(values); return 2; }
    k=rice_i32_parameter(values,count,WITNESS_MAX_RICE_K,&bits);
    bytes=(size_t)((bits+7)/8); out[off++]=(uint8_t)k;
    if(rice_i32_encode(&out[off],bytes,values,count,k)) { free(values); return 3; }
    free(values); off+=bytes;
  }
  return off==outlen ? 0 : 4;
}

int witness_deserialize_contextual(witness *wt, const witness *shape,
                                   const uint8_t *in, size_t inlen,
                                   size_t *consumed) {
  size_t i,j,k_index,off=0,count,bytes,sum_n=0;
  uint64_t bits;
  unsigned k,best_k;
  size_t *n=NULL;
  int32_t *values=NULL;

  if(!wt || !shape || !in || !consumed || !shape->r || !shape->n) return 1;
  n=malloc(shape->r*sizeof(*n));
  if(!n) return 2;
  for(i=0;i<shape->r;i++) {
    if(!shape->n[i] || shape->n[i]>SIZE_MAX/N || shape->n[i]>SIZE_MAX-sum_n) goto err;
    n[i]=shape->n[i]; sum_n+=n[i];
  }
  if(sum_n>SIZE_MAX/sizeof(poly)) goto err;
  memset(wt,0,sizeof(*wt)); init_witness_raw(wt,shape->r,n);
  for(i=0;i<shape->r;i++) {
    if(off>=inlen) goto decode_err;
    k=in[off++]; if(k>WITNESS_MAX_RICE_K) goto decode_err;
    count=n[i]*N; values=malloc(count*sizeof(*values));
    if(!values || rice_i32_decode_prefix(values,count,&in[off],inlen-off,k,
                                          UINT16_MAX,&bytes)) goto decode_err;
    best_k=rice_i32_parameter(values,count,WITNESS_MAX_RICE_K,&bits);
    if(best_k!=k || (bits+7)/8!=bytes) goto decode_err;
    for(j=0;j<n[i];j++)
      for(k_index=0;k_index<N;k_index++)
        wt->s[i][j].vec->c[k_index]=(int16_t)values[j*N+k_index];
    free(values); values=NULL; off+=bytes;
    wt->normsq[i]=polyvec_sprodz(wt->s[i],wt->s[i],wt->n[i]);
  }
  free(n); *consumed=off; return 0;
decode_err:
  free(values); free_witness(wt);
err:
  free(n); return 3;
}
