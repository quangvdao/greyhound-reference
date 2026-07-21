#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "malloc.h"
#include "labrador.h"

#define PROOF_WIRE_HEADER_BYTES 124
#define PROOF_WIRE_VERSION 3
#define JL_COORDINATES 256
#define MAX_RICE_K 31

typedef struct {
  uint8_t *bytes;
  size_t bit;
  size_t capacity;
} bitwriter;

typedef struct {
  const uint8_t *bytes;
  size_t bit;
  size_t capacity;
} bitreader;

static void put16(uint8_t *out, uint16_t x) {
  out[0] = (uint8_t)x;
  out[1] = (uint8_t)(x >> 8);
}

static void put32(uint8_t *out, uint32_t x) {
  size_t i;
  for(i=0;i<4;i++) out[i] = (uint8_t)(x >> (8*i));
}

static void put64(uint8_t *out, uint64_t x) {
  size_t i;
  for(i=0;i<8;i++) out[i] = (uint8_t)(x >> (8*i));
}

static uint16_t get16(const uint8_t *in) {
  return (uint16_t)in[0] | (uint16_t)in[1] << 8;
}

static uint32_t get32(const uint8_t *in) {
  size_t i;
  uint32_t x = 0;
  for(i=0;i<4;i++) x |= (uint32_t)in[i] << (8*i);
  return x;
}

static uint64_t get64(const uint8_t *in) {
  size_t i;
  uint64_t x = 0;
  for(i=0;i<8;i++) x |= (uint64_t)in[i] << (8*i);
  return x;
}

static int add_size(size_t *total, size_t count, size_t width) {
  if(count > (SIZE_MAX-*total)/width) return 1;
  *total += count*width;
  return 0;
}

static uint32_t zigzag32(int32_t x) {
  if(x >= 0) return (uint32_t)x << 1;
  return (uint32_t)(2*(-(int64_t)x)-1);
}

static int32_t unzigzag32(uint32_t x) {
  if((x & 1) == 0) return (int32_t)(x >> 1);
  return (int32_t)(-(int64_t)(x >> 1)-1);
}

static uint64_t rice_bit_count(const int32_t p[JL_COORDINATES], unsigned k) {
  size_t i;
  uint64_t bits = 0;

  for(i=0;i<JL_COORDINATES;i++)
    bits += (uint64_t)(zigzag32(p[i]) >> k)+1+k;
  return bits;
}

static unsigned rice_parameter(const int32_t p[JL_COORDINATES], uint64_t *bits) {
  unsigned k,best_k = 0;
  uint64_t n,best = rice_bit_count(p,0);

  for(k=1;k<=MAX_RICE_K;k++) {
    n = rice_bit_count(p,k);
    if(n < best) {
      best = n;
      best_k = k;
    }
  }
  if(bits != NULL) *bits = best;
  return best_k;
}

static int put_bit(bitwriter *w, unsigned value) {
  if(w->bit >= w->capacity) return 1;
  if(value != 0) w->bytes[w->bit >> 3] |= (uint8_t)(1u << (w->bit & 7));
  w->bit++;
  return 0;
}

static int get_bit(bitreader *r, unsigned *value) {
  if(r->bit >= r->capacity) return 1;
  *value = (r->bytes[r->bit >> 3] >> (r->bit & 7)) & 1u;
  r->bit++;
  return 0;
}

static int rice_encode(uint8_t *out, size_t outlen,
                       const int32_t p[JL_COORDINATES], unsigned k) {
  size_t i,j;
  uint32_t x,q;
  bitwriter w = {out,0,outlen*8};

  memset(out,0,outlen);
  for(i=0;i<JL_COORDINATES;i++) {
    x = zigzag32(p[i]);
    q = x >> k;
    for(j=0;j<q;j++) if(put_bit(&w,1)) return 1;
    if(put_bit(&w,0)) return 1;
    for(j=0;j<k;j++) if(put_bit(&w,(x >> j) & 1u)) return 1;
  }
  return w.capacity-w.bit >= 8;
}

static int rice_decode_prefix(int32_t p[JL_COORDINATES], const uint8_t *in,
                              size_t inlen, unsigned k, size_t *consumed) {
  size_t i,j;
  uint64_t q,x;
  unsigned bit;
  bitreader r = {in,0,inlen*8};

  for(i=0;i<JL_COORDINATES;i++) {
    q = 0;
    do {
      if(get_bit(&r,&bit)) return 1;
      if(bit != 0 && ++q > (UINT32_MAX >> k)) return 1;
    } while(bit != 0);
    x = q << k;
    for(j=0;j<k;j++) {
      if(get_bit(&r,&bit)) return 1;
      x |= (uint64_t)bit << j;
    }
    if(x > UINT32_MAX) return 1;
    p[i] = unzigzag32((uint32_t)x);
  }
  *consumed = (r.bit+7)/8;
  while(r.bit < *consumed*8) {
    if(get_bit(&r,&bit) || bit != 0) return 1;
  }
  return 0;
}

static int rice_decode(int32_t p[JL_COORDINATES], const uint8_t *in,
                       size_t inlen, unsigned k) {
  size_t consumed;

  if(rice_decode_prefix(p,in,inlen,k,&consumed)) return 1;
  return consumed != inlen;
}

size_t proof_serialized_size(const proof *pi) {
  size_t total = PROOF_WIRE_HEADER_BYTES;
  size_t payload;
  uint64_t jlbits,jlbytes;

  if(pi == NULL) return 0;
  if(add_size(&total,pi->r,16)) return 0;
  if(pi->cpp->u1len > SIZE_MAX-pi->cpp->u2len-LIFTS) return 0;
  payload = pi->cpp->u1len+pi->cpp->u2len+LIFTS;
  rice_parameter(pi->p,&jlbits);
  jlbytes = (jlbits+7)/8;
  if(jlbytes > SIZE_MAX || add_size(&total,(size_t)jlbytes,1)) return 0;
  if(add_size(&total,payload,N*QBYTES)) return 0;
  return total;
}

int proof_serialize(uint8_t *out, size_t outlen, const proof *pi) {
  size_t i,off,size,jlbytes;
  size_t params[10];
  uint64_t jlbits;
  unsigned rice_k;

  if(pi == NULL || pi->n == NULL || pi->u1 == NULL ||
     pi->foldnonce >= FOLD_GRIND_MAX_ATTEMPTS) return 1;
  const size_t values[10] = {
    pi->cpp->f,pi->cpp->fu,pi->cpp->fg,
    pi->cpp->b,pi->cpp->bu,pi->cpp->bg,
    pi->cpp->kappa,pi->cpp->kappa1,pi->cpp->u1len,pi->cpp->u2len
  };
  memcpy(params,values,sizeof(params));

  size = proof_serialized_size(pi);
  if(size == 0 || out == NULL || outlen != size) return 1;
  rice_k = rice_parameter(pi->p,&jlbits);
  jlbytes = (size_t)((jlbits+7)/8);

  memcpy(out,"LBRP",4);
  put16(&out[4],PROOF_WIRE_VERSION);
  put16(&out[6],N);
  put16(&out[8],LOGQ);
  put16(&out[10],LIFTS);
  out[12] = (uint8_t)(pi->tail != 0);
  out[13] = (uint8_t)rice_k;
  memset(&out[14],0,2);
  put64(&out[16],pi->r);
  off = 24;
  for(i=0;i<10;i++,off+=8) put64(&out[off],params[i]);
  put64(&out[104],pi->jlnonce);
  put64(&out[112],pi->normsq);
  put32(&out[120],pi->foldnonce);

  off = PROOF_WIRE_HEADER_BYTES;
  for(i=0;i<pi->r;i++,off+=16) {
    put64(&out[off],pi->n[i]);
    put64(&out[off+8],pi->nu[i]);
  }
  if(rice_encode(&out[off],jlbytes,pi->p,rice_k)) return 1;
  off += jlbytes;
  polzvec_bitpack(&out[off],pi->u1,pi->cpp->u1len);
  off += pi->cpp->u1len*N*QBYTES;
  polzvec_bitpack(&out[off],pi->u2,pi->cpp->u2len);
  off += pi->cpp->u2len*N*QBYTES;
  polzvec_bitpack(&out[off],pi->bb,LIFTS);
  return 0;
}

int proof_deserialize(proof *pi, const uint8_t *in, size_t inlen) {
  size_t i,off,size,payload,fixed,jlbytes;
  uint64_t r64,param[10];
  unsigned rice_k;
  void *buf;

  if(pi == NULL || in == NULL || inlen < PROOF_WIRE_HEADER_BYTES) return 1;
  if(memcmp(in,"LBRP",4) != 0 || get16(&in[4]) != PROOF_WIRE_VERSION) return 2;
  if(get16(&in[6]) != N || get16(&in[8]) != LOGQ || get16(&in[10]) != LIFTS) return 3;
  if(in[12] > 1 || in[13] > MAX_RICE_K || in[14] != 0 || in[15] != 0) return 4;
  rice_k = in[13];

  r64 = get64(&in[16]);
  if(r64 == 0 || r64 > SIZE_MAX) return 5;
  off = 24;
  for(i=0;i<10;i++,off+=8) {
    param[i] = get64(&in[off]);
    if(param[i] > SIZE_MAX) return 5;
  }
  if(param[0] == 0 || param[1] == 0 || param[3] > LOGQ || param[4] > LOGQ ||
     param[5] > LOGQ || param[6] == 0 || param[6] > SIS_MAX_RANK ||
     param[7] > SIS_MAX_RANK)
    return 6;

  memset(pi,0,sizeof(*pi));
  pi->r = (size_t)r64;
  pi->tail = in[12];
  pi->cpp->f = (size_t)param[0];
  pi->cpp->fu = (size_t)param[1];
  pi->cpp->fg = (size_t)param[2];
  pi->cpp->b = (size_t)param[3];
  pi->cpp->bu = (size_t)param[4];
  pi->cpp->bg = (size_t)param[5];
  pi->cpp->kappa = (size_t)param[6];
  pi->cpp->kappa1 = (size_t)param[7];
  pi->cpp->u1len = (size_t)param[8];
  pi->cpp->u2len = (size_t)param[9];
  pi->jlnonce = (size_t)get64(&in[104]);
  pi->normsq = get64(&in[112]);
  pi->foldnonce = get32(&in[120]);
  if(pi->foldnonce >= FOLD_GRIND_MAX_ATTEMPTS)
    return 6;

  fixed = PROOF_WIRE_HEADER_BYTES;
  if(add_size(&fixed,pi->r,16)) return 7;
  if(pi->cpp->u1len > SIZE_MAX-pi->cpp->u2len-LIFTS) return 7;
  payload = pi->cpp->u1len+pi->cpp->u2len+LIFTS;
  if(add_size(&fixed,payload,N*QBYTES) || fixed >= inlen) return 7;
  jlbytes = inlen-fixed;
  buf = _malloc(2*pi->r*sizeof(size_t));
  pi->n = buf;
  pi->nu = &pi->n[pi->r];
  buf = _aligned_alloc(64,payload*sizeof(polz));
  pi->u1 = buf;
  pi->u2 = &pi->u1[pi->cpp->u1len];
  pi->bb = &pi->u2[pi->cpp->u2len];

  off = PROOF_WIRE_HEADER_BYTES;
  for(i=0;i<pi->r;i++,off+=16) {
    pi->n[i] = (size_t)get64(&in[off]);
    pi->nu[i] = (size_t)get64(&in[off+8]);
  }
  if(rice_decode(pi->p,&in[off],jlbytes,rice_k) || rice_parameter(pi->p,NULL) != rice_k) {
    free_proof(pi);
    return 8;
  }
  size = proof_serialized_size(pi);
  if(size != inlen) {
    free_proof(pi);
    return 8;
  }
  off += jlbytes;
  for(i=0;i<payload;i++,off+=N*QBYTES) polz_bitunpack(&pi->u1[i],&in[off]);
  return 0;
}

size_t proof_contextual_serialized_size(const proof *pi) {
  size_t total=21,payload;
  uint64_t bits;

  if(!pi || !pi->n || !pi->u1) return 0;
  if(pi->cpp->u1len > SIZE_MAX-pi->cpp->u2len-LIFTS) return 0;
  payload=pi->cpp->u1len+pi->cpp->u2len+LIFTS;
  rice_parameter(pi->p,&bits);
  if((bits+7)/8 > SIZE_MAX-total) return 0;
  total+=(size_t)((bits+7)/8);
  if(add_size(&total,payload,N*QBYTES)) return 0;
  return total;
}

int proof_serialize_contextual(uint8_t *out, size_t outlen, const proof *pi) {
  size_t off,jlbytes;
  uint64_t bits;
  unsigned k;

  if(!out || outlen != proof_contextual_serialized_size(pi) ||
     pi->foldnonce >= FOLD_GRIND_MAX_ATTEMPTS) return 1;
  k=rice_parameter(pi->p,&bits); jlbytes=(size_t)((bits+7)/8);
  out[0]=(uint8_t)k; put64(&out[1],pi->jlnonce); put64(&out[9],pi->normsq);
  off=17;
  put32(&out[off],pi->foldnonce);
  off+=4;
  if(rice_encode(&out[off],jlbytes,pi->p,k)) return 2;
  off+=jlbytes;
  polzvec_bitpack(&out[off],pi->u1,pi->cpp->u1len);
  off+=pi->cpp->u1len*N*QBYTES;
  polzvec_bitpack(&out[off],pi->u2,pi->cpp->u2len);
  off+=pi->cpp->u2len*N*QBYTES;
  polzvec_bitpack(&out[off],pi->bb,LIFTS);
  return 0;
}

int proof_deserialize_contextual(proof *pi, const proof *shape,
                                 const uint8_t *in, size_t inlen,
                                 size_t *consumed) {
  size_t i,jlbytes,off,payload,total;
  uint64_t nonce;
  unsigned k;
  void *buf;

  if(!pi || !shape || !in || !consumed || !shape->r || !shape->n ||
     !shape->nu || inlen<21) return 1;
  k=in[0]; nonce=get64(&in[1]);
  if(k>MAX_RICE_K || nonce>SIZE_MAX ||
     shape->r>SIZE_MAX/(2*sizeof(size_t)) ||
     shape->cpp->u1len > SIZE_MAX-shape->cpp->u2len-LIFTS) return 2;
  payload=shape->cpp->u1len+shape->cpp->u2len+LIFTS;
  off=17;
  if(get32(&in[off]) >= FOLD_GRIND_MAX_ATTEMPTS) return 2;
  off+=4;
  if(payload > (SIZE_MAX-off)/(N*QBYTES)) return 2;

  memset(pi,0,sizeof(*pi)); pi->r=shape->r; pi->tail=shape->tail;
  *pi->cpp=*shape->cpp; pi->jlnonce=(size_t)nonce; pi->normsq=get64(&in[9]);
  pi->foldnonce=get32(&in[17]);
  buf=_malloc(2*pi->r*sizeof(size_t)); pi->n=buf; pi->nu=&pi->n[pi->r];
  memcpy(pi->n,shape->n,pi->r*sizeof(size_t));
  memcpy(pi->nu,shape->nu,pi->r*sizeof(size_t));
  if(rice_decode_prefix(pi->p,&in[off],inlen-off,k,&jlbytes) ||
     rice_parameter(pi->p,NULL)!=k) goto err;
  if(payload*N*QBYTES>SIZE_MAX-off-jlbytes) goto err;
  total=off+jlbytes+payload*N*QBYTES;
  if(total>inlen) goto err;
  buf=_aligned_alloc(64,payload*sizeof(polz));
  pi->u1=buf; pi->u2=&pi->u1[pi->cpp->u1len]; pi->bb=&pi->u2[pi->cpp->u2len];
  off+=jlbytes;
  for(i=0;i<payload;i++,off+=N*QBYTES) polz_bitunpack(&pi->u1[i],&in[off]);
  *consumed=total;
  return 0;
err:
  free_proof(pi); return 3;
}
