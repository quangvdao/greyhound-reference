#ifndef RICE_H
#define RICE_H

#include <stddef.h>
#include <stdint.h>

unsigned rice_i32_parameter(const int32_t *values, size_t len, unsigned max_k,
                            uint64_t *bits);
int rice_i32_encode(uint8_t *out, size_t outlen, const int32_t *values,
                    size_t len, unsigned k);
int rice_i32_decode(int32_t *values, size_t len, const uint8_t *in,
                    size_t inlen, unsigned k, uint32_t max_zigzag);
int rice_i32_decode_prefix(int32_t *values, size_t len, const uint8_t *in,
                           size_t inlen, unsigned k, uint32_t max_zigzag,
                           size_t *consumed);

#endif
