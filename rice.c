#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include "rice.h"

typedef struct { uint8_t *bytes; size_t bit,capacity; } bitwriter;
typedef struct { const uint8_t *bytes; size_t bit,capacity; } bitreader;

static uint32_t zigzag32(int32_t x) {
  return x >= 0 ? (uint32_t)x << 1 : (uint32_t)(2*(-(int64_t)x)-1);
}

static int32_t unzigzag32(uint32_t x) {
  return (x & 1) ? (int32_t)(-(int64_t)(x >> 1)-1) : (int32_t)(x >> 1);
}

static int put_bit(bitwriter *w, unsigned value) {
  if(w->bit >= w->capacity) return 1;
  if(value) w->bytes[w->bit >> 3] |= (uint8_t)(1u << (w->bit & 7));
  w->bit++;
  return 0;
}

static int get_bit(bitreader *r, unsigned *value) {
  if(r->bit >= r->capacity) return 1;
  *value = (r->bytes[r->bit >> 3] >> (r->bit & 7)) & 1u;
  r->bit++;
  return 0;
}

static uint64_t bit_count(const int32_t *values, size_t len, unsigned k) {
  size_t i;
  uint64_t bits = 0;
  for(i=0;i<len;i++) bits += (uint64_t)(zigzag32(values[i]) >> k)+1+k;
  return bits;
}

unsigned rice_i32_parameter(const int32_t *values, size_t len, unsigned max_k,
                            uint64_t *bits) {
  unsigned k,best_k = 0;
  uint64_t n,best = bit_count(values,len,0);
  for(k=1;k<=max_k;k++) {
    n = bit_count(values,len,k);
    if(n < best) { best = n; best_k = k; }
  }
  if(bits) *bits = best;
  return best_k;
}

int rice_i32_encode(uint8_t *out, size_t outlen, const int32_t *values,
                    size_t len, unsigned k) {
  size_t i,j;
  uint32_t x,q;
  bitwriter w = {out,0,outlen*8};
  memset(out,0,outlen);
  for(i=0;i<len;i++) {
    x = zigzag32(values[i]);
    q = x >> k;
    for(j=0;j<q;j++) if(put_bit(&w,1)) return 1;
    if(put_bit(&w,0)) return 1;
    for(j=0;j<k;j++) if(put_bit(&w,(x >> j) & 1u)) return 1;
  }
  return w.capacity-w.bit >= 8;
}

int rice_i32_decode_prefix(int32_t *values, size_t len, const uint8_t *in,
                           size_t inlen, unsigned k, uint32_t max_zigzag,
                           size_t *consumed) {
  size_t i,j;
  uint32_t q,x;
  unsigned bit;
  bitreader r = {in,0,inlen*8};
  for(i=0;i<len;i++) {
    q = 0;
    do {
      if(get_bit(&r,&bit)) return 1;
      if(bit && ++q > (max_zigzag >> k)) return 1;
    } while(bit);
    x = q << k;
    for(j=0;j<k;j++) {
      if(get_bit(&r,&bit)) return 1;
      x |= (uint32_t)bit << j;
    }
    if(x > max_zigzag) return 1;
    values[i] = unzigzag32(x);
  }
  *consumed = (r.bit+7)/8;
  while(r.bit < *consumed*8) if(get_bit(&r,&bit) || bit) return 1;
  return 0;
}

int rice_i32_decode(int32_t *values, size_t len, const uint8_t *in,
                    size_t inlen, unsigned k, uint32_t max_zigzag) {
  size_t consumed;

  if(rice_i32_decode_prefix(values,len,in,inlen,k,max_zigzag,&consumed)) return 1;
  return consumed != inlen;
}
