#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "malloc.h"
#include "greyhound.h"

#define POLCOM_WIRE_HEADER_BYTES 112
#define POLCOM_WIRE_VERSION 2

static void put16(uint8_t *out, uint16_t x) {
  out[0] = (uint8_t)x; out[1] = (uint8_t)(x >> 8);
}
static void put32(uint8_t *out, uint32_t x) {
  size_t i; for(i=0;i<4;i++) out[i] = (uint8_t)(x >> (8*i));
}
static void put64(uint8_t *out, uint64_t x) {
  size_t i; for(i=0;i<8;i++) out[i] = (uint8_t)(x >> (8*i));
}
static uint16_t get16(const uint8_t *in) {
  return (uint16_t)in[0] | (uint16_t)in[1] << 8;
}
static uint32_t get32(const uint8_t *in) {
  size_t i; uint32_t x = 0;
  for(i=0;i<4;i++) x |= (uint32_t)in[i] << (8*i);
  return x;
}
static uint64_t get64(const uint8_t *in) {
  size_t i; uint64_t x = 0;
  for(i=0;i<8;i++) x |= (uint64_t)in[i] << (8*i);
  return x;
}

size_t polcomprf_serialized_size(const polcomprf *pi) {
  size_t payload;
  if(!pi || !pi->u1 || !pi->u2 || !pi->cpp->kappa1) return 0;
  if(pi->cpp->kappa1 > (SIZE_MAX-POLCOM_WIRE_HEADER_BYTES)/(2*N*QBYTES)) return 0;
  payload = 2*pi->cpp->kappa1*N*QBYTES;
  return POLCOM_WIRE_HEADER_BYTES+payload;
}

int polcomprf_serialize(uint8_t *out, size_t outlen, const polcomprf *pi) {
  size_t i,off;
  const size_t values[6] = {pi->cpp->f,pi->cpp->fu,pi->cpp->b,
                            pi->cpp->bu,pi->cpp->kappa,pi->cpp->kappa1};
  if(!out || outlen != polcomprf_serialized_size(pi) ||
     pi->foldnonce >= FOLD_GRIND_MAX_ATTEMPTS) return 1;
  memcpy(out,"GHPR",4);
  put16(&out[4],POLCOM_WIRE_VERSION); put16(&out[6],N); put16(&out[8],LOGQ);
  memset(&out[10],0,2); put32(&out[12],pi->foldnonce);
  put64(&out[16],pi->len); put64(&out[24],pi->m); put64(&out[32],pi->n);
  put64(&out[40],(uint64_t)pi->x); put64(&out[48],(uint64_t)pi->y);
  off = 56;
  for(i=0;i<6;i++,off+=8) put64(&out[off],values[i]);
  put64(&out[104],pi->normsq);
  off = POLCOM_WIRE_HEADER_BYTES;
  polzvec_bitpack(&out[off],pi->u1,pi->cpp->kappa1);
  off += pi->cpp->kappa1*N*QBYTES;
  polzvec_bitpack(&out[off],pi->u2,pi->cpp->kappa1);
  return 0;
}

int polcomprf_deserialize(polcomprf *pi, const uint8_t *in, size_t inlen) {
  size_t i,off,kappa1,payload;
  uint64_t value[6],len,m,n;
  if(!pi || !in || inlen < POLCOM_WIRE_HEADER_BYTES) return 1;
  if(memcmp(in,"GHPR",4) || get16(&in[4]) != POLCOM_WIRE_VERSION ||
     get16(&in[6]) != N || get16(&in[8]) != LOGQ) return 2;
  if(in[10] || in[11] || get32(&in[12]) >= FOLD_GRIND_MAX_ATTEMPTS) return 2;
  len=get64(&in[16]); m=get64(&in[24]); n=get64(&in[32]);
  if(!len || !m || !n || len > SIZE_MAX || m > SIZE_MAX || n > SIZE_MAX) return 3;
  off=56;
  for(i=0;i<6;i++,off+=8) { value[i]=get64(&in[off]); if(value[i] > SIZE_MAX) return 3; }
  kappa1=(size_t)value[5];
  if(!value[0] || !value[1] || !value[2] || !value[3] || !value[4] || !kappa1 ||
     value[2] > LOGQ || value[3] > LOGQ || value[4] > SIS_MAX_RANK || kappa1 > SIS_MAX_RANK)
    return 4;
  if(kappa1 > (SIZE_MAX-POLCOM_WIRE_HEADER_BYTES)/(2*N*QBYTES)) return 4;
  payload=2*kappa1*N*QBYTES;
  if(inlen != POLCOM_WIRE_HEADER_BYTES+payload) return 5;
  if(value[4] > SIZE_MAX/value[1] || (size_t)(value[4]*value[1]) > SIZE_MAX/(size_t)n) return 4;
  memset(pi,0,sizeof(*pi));
  pi->len=(size_t)len; pi->m=(size_t)m; pi->n=(size_t)n;
  pi->x=(int64_t)get64(&in[40]); pi->y=(int64_t)get64(&in[48]);
  pi->cpp->f=(size_t)value[0]; pi->cpp->fu=(size_t)value[1];
  pi->cpp->b=(size_t)value[2]; pi->cpp->bu=(size_t)value[3];
  pi->cpp->kappa=(size_t)value[4]; pi->cpp->kappa1=kappa1;
  pi->cpp->u1len=pi->cpp->kappa*pi->cpp->fu*pi->n;
  if(pi->cpp->fu > SIZE_MAX/pi->n) return 4;
  pi->cpp->u2len=pi->cpp->fu*pi->n;
  pi->foldnonce=get32(&in[12]);
  pi->normsq=get64(&in[104]);
  pi->u1=_aligned_alloc(64,kappa1*sizeof(polz));
  pi->u2=_aligned_alloc(64,kappa1*sizeof(polz));
  if(!pi->u1 || !pi->u2) { free(pi->u1); free(pi->u2); memset(pi,0,sizeof(*pi)); return 6; }
  off=POLCOM_WIRE_HEADER_BYTES;
  for(i=0;i<kappa1;i++,off+=N*QBYTES) polz_bitunpack(&pi->u1[i],&in[off]);
  for(i=0;i<kappa1;i++,off+=N*QBYTES) polz_bitunpack(&pi->u2[i],&in[off]);
  return 0;
}

size_t polcomprf_contextual_serialized_size(const polcomprf *pi) {
  if(!pi || !pi->u2 || !pi->cpp->kappa1) return 0;
  if(pi->cpp->kappa1 > (SIZE_MAX-12)/(N*QBYTES)) return 0;
  return 12+pi->cpp->kappa1*N*QBYTES;
}

int polcomprf_serialize_contextual(uint8_t *out, size_t outlen,
                                   const polcomprf *pi) {
  if(!out || outlen != polcomprf_contextual_serialized_size(pi) ||
     pi->foldnonce >= FOLD_GRIND_MAX_ATTEMPTS) return 1;
  put64(out,pi->normsq);
  put32(&out[8],pi->foldnonce);
  polzvec_bitpack(&out[12],pi->u2,pi->cpp->kappa1);
  return 0;
}

int polcomprf_deserialize_contextual(polcomprf *pi, const polcomprf *context,
                                     const uint8_t *in, size_t inlen,
                                     size_t *consumed) {
  size_t size,kappa1,i,off;

  if(!pi || !context || !in || !consumed || !context->u1) return 1;
  kappa1=context->cpp->kappa1;
  if(!kappa1 || kappa1 > (SIZE_MAX-12)/(N*QBYTES) || inlen<12 ||
     get32(&in[8]) >= FOLD_GRIND_MAX_ATTEMPTS) return 2;
  size=12+kappa1*N*QBYTES;
  if(inlen<size) return 3;
  memset(pi,0,sizeof(*pi));
  pi->len=context->len; pi->m=context->m; pi->n=context->n;
  pi->x=context->x; pi->y=context->y; *pi->cpp=*context->cpp;
  pi->foldnonce=get32(&in[8]);
  pi->normsq=get64(in);
  pi->u1=_aligned_alloc(64,kappa1*sizeof(polz));
  pi->u2=_aligned_alloc(64,kappa1*sizeof(polz));
  if(!pi->u1 || !pi->u2) {
    free(pi->u1); free(pi->u2); memset(pi,0,sizeof(*pi)); return 4;
  }
  polzvec_copy(pi->u1,context->u1,kappa1);
  off=12;
  for(i=0;i<kappa1;i++,off+=N*QBYTES) polz_bitunpack(&pi->u2[i],&in[off]);
  *consumed=size;
  return 0;
}
