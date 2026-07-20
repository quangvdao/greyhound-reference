#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "pack.h"

#define PACK_WIRE_HEADER_BYTES 24
#define PACK_WIRE_VERSION 1

static void put16(uint8_t *out, uint16_t x) { out[0]=(uint8_t)x; out[1]=(uint8_t)(x>>8); }
static void put64(uint8_t *out, uint64_t x) { size_t i; for(i=0;i<8;i++) out[i]=(uint8_t)(x>>(8*i)); }
static uint16_t get16(const uint8_t *in) { return (uint16_t)in[0] | (uint16_t)in[1]<<8; }
static uint64_t get64(const uint8_t *in) { size_t i; uint64_t x=0; for(i=0;i<8;i++) x|=(uint64_t)in[i]<<(8*i); return x; }

size_t greyhound_pack_serialized_size(const polcomprf *top, const composite *p) {
  size_t i,total,n;
  if(!top || !p || !p->l || p->l > 16) return 0;
  total=PACK_WIRE_HEADER_BYTES+8*(p->l+2);
  n=polcomprf_serialized_size(top);
  if(!n || n>SIZE_MAX-total) return 0;
  total+=n;
  for(i=0;i<p->l;i++) {
    n=proof_serialized_size(p->pi[i]);
    if(!n || n>SIZE_MAX-total) return 0;
    total+=n;
  }
  n=witness_serialized_size(&p->owt);
  if(!n || n>SIZE_MAX-total) return 0;
  return total+n;
}

int greyhound_pack_serialize(uint8_t *out, size_t outlen, const polcomprf *top,
                             const composite *p) {
  size_t i,off,table,n;
  if(!out || outlen != greyhound_pack_serialized_size(top,p)) return 1;
  memcpy(out,"GHCP",4); put16(&out[4],PACK_WIRE_VERSION); put16(&out[6],N);
  put16(&out[8],LOGQ); memset(&out[10],0,6); put64(&out[16],p->l);
  table=PACK_WIRE_HEADER_BYTES; off=table+8*(p->l+2);
  n=polcomprf_serialized_size(top); put64(&out[table],n); table+=8;
  if(polcomprf_serialize(&out[off],n,top)) return 2; off+=n;
  for(i=0;i<p->l;i++) {
    n=proof_serialized_size(p->pi[i]); put64(&out[table],n); table+=8;
    if(proof_serialize(&out[off],n,p->pi[i])) return 3; off+=n;
  }
  n=witness_serialized_size(&p->owt); put64(&out[table],n);
  if(witness_serialize(&out[off],n,&p->owt)) return 4;
  return off+n == outlen ? 0 : 5;
}

int greyhound_pack_deserialize(polcomprf *top, composite *p,
                               const uint8_t *in, size_t inlen) {
  size_t i,l,off,table,n,total;
  uint64_t l64,n64;
  if(!top || !p || !in || inlen < PACK_WIRE_HEADER_BYTES) return 1;
  if(memcmp(in,"GHCP",4) || get16(&in[4])!=PACK_WIRE_VERSION ||
     get16(&in[6])!=N || get16(&in[8])!=LOGQ) return 2;
  for(i=10;i<16;i++) if(in[i]) return 2;
  l64=get64(&in[16]);
  if(!l64 || l64>16 || l64>SIZE_MAX) return 3;
  l=(size_t)l64; table=PACK_WIRE_HEADER_BYTES;
  if(l+2 > (inlen-table)/8) return 3;
  off=table+8*(l+2); total=off;
  for(i=0;i<l+2;i++) {
    n64=get64(&in[table+8*i]);
    if(!n64 || n64>SIZE_MAX || (size_t)n64>inlen-total) return 4;
    total+=(size_t)n64;
  }
  if(total!=inlen) return 4;
  memset(top,0,sizeof(*top)); memset(p,0,sizeof(*p));
  n=(size_t)get64(&in[table]);
  if(polcomprf_deserialize(top,&in[off],n)) goto err;
  off+=n; table+=8;
  for(i=0;i<l;i++) {
    n=(size_t)get64(&in[table]); table+=8;
    p->pi[i]=malloc(sizeof(*p->pi[i]));
    if(!p->pi[i]) goto err;
    memset(p->pi[i],0,sizeof(*p->pi[i]));
    if(proof_deserialize(p->pi[i],&in[off],n)) goto err;
    p->l++; off+=n;
  }
  n=(size_t)get64(&in[table]);
  if(witness_deserialize(&p->owt,&in[off],n)) goto err;
  p->tail_size=(double)n/1024;
  p->size=(double)inlen/1024;
  return 0;
err:
  free_polcomprf(top); free_composite(p); return 5;
}

size_t greyhound_pack_contextual_serialized_size(const polcomprf *top,
                                                const composite *p) {
  size_t i,total,n;

  if(!top || !p || !p->l || p->l>16) return 0;
  total=polcomprf_contextual_serialized_size(top);
  if(!total) return 0;
  for(i=0;i<p->l;i++) {
    n=proof_contextual_serialized_size(p->pi[i]);
    if(!n || n>SIZE_MAX-total) return 0;
    total+=n;
  }
  n=witness_contextual_serialized_size(&p->owt);
  if(!n || n>SIZE_MAX-total) return 0;
  return total+n;
}

int greyhound_pack_serialize_contextual(uint8_t *out, size_t outlen,
                                      const polcomprf *top,
                                      const composite *p) {
  size_t i,off=0,n;

  if(!out || outlen!=greyhound_pack_contextual_serialized_size(top,p)) return 1;
  n=polcomprf_contextual_serialized_size(top);
  if(polcomprf_serialize_contextual(&out[off],n,top)) return 2;
  off+=n;
  for(i=0;i<p->l;i++) {
    n=proof_contextual_serialized_size(p->pi[i]);
    if(proof_serialize_contextual(&out[off],n,p->pi[i])) return 3;
    off+=n;
  }
  n=witness_contextual_serialized_size(&p->owt);
  if(witness_serialize_contextual(&out[off],n,&p->owt)) return 4;
  return off+n==outlen ? 0 : 5;
}

int greyhound_pack_deserialize_contextual(polcomprf *top, composite *p,
                                        const polcomprf *top_context,
                                        const composite *p_context,
                                        const uint8_t *in, size_t inlen) {
  size_t i,off=0,n,tail_start;

  if(!top || !p || !top_context || !p_context || !in ||
     !p_context->l || p_context->l>16) return 1;
  memset(top,0,sizeof(*top)); memset(p,0,sizeof(*p));
  if(polcomprf_deserialize_contextual(top,top_context,in,inlen,&n)) goto err;
  off+=n;
  for(i=0;i<p_context->l;i++) {
    p->pi[i]=malloc(sizeof(*p->pi[i]));
    if(!p->pi[i]) goto err;
    memset(p->pi[i],0,sizeof(*p->pi[i])); p->l=i+1;
    if(proof_deserialize_contextual(p->pi[i],p_context->pi[i],&in[off],
                                  inlen-off,&n)) goto err;
    off+=n;
  }
  tail_start=off;
  if(witness_deserialize_contextual(&p->owt,&p_context->owt,&in[off],
                                  inlen-off,&n)) goto err;
  off+=n;
  if(off!=inlen) goto err;
  p->tail_size=(double)(off-tail_start)/1024;
  p->size=(double)off/1024;
  return 0;
err:
  free_polcomprf(top); free_composite(p); return 2;
}
