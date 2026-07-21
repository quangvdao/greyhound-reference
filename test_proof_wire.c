#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "labrador.h"

int main(void) {
  size_t i,size,consumed,payload = 3+2+LIFTS;
  uint8_t seed[16] = {1};
  proof a = {},b = {};
  uint8_t *wire1,*wire2,*wire3;

  a.r = 3;
  a.tail = 1;
  a.n = malloc(2*a.r*sizeof(size_t));
  a.nu = &a.n[a.r];
  a.u1 = aligned_alloc(64,payload*sizeof(polz));
  if(a.n == NULL || a.u1 == NULL) return 1;
  a.u2 = &a.u1[3];
  a.bb = &a.u2[2];
  a.cpp->f = 2;
  a.cpp->fu = 6;
  a.cpp->fg = 1;
  a.cpp->b = 6;
  a.cpp->bu = 5;
  a.cpp->bg = 6;
  a.cpp->kappa = 14;
  a.cpp->kappa1 = 3;
  a.cpp->u1len = 3;
  a.cpp->u2len = 2;
  a.jlnonce = 17;
  a.foldnonce = 23;
  a.normsq = UINT64_C(0x123456789abcdef0);
  for(i=0;i<a.r;i++) {
    a.n[i] = 100+i;
    a.nu[i] = i == a.r-1 ? 7 : 0;
  }
  for(i=0;i<256;i++) a.p[i] = (int32_t)i-128;
  a.p[0] = INT32_MIN;
  a.p[1] = INT32_MAX;
  polzvec_almostuniform(a.u1,payload,seed,9);

  size = proof_serialized_size(&a);
  if(size >= 124+16*a.r+4*256+payload*N*QBYTES) return 2;
  wire1 = malloc(size);
  wire2 = malloc(size);
  if(wire1 == NULL || wire2 == NULL) return 3;
  if(proof_serialize(wire1,size,&a) != 0) return 4;
  if(proof_deserialize(&b,wire1,size) != 0) return 5;
  if(proof_serialized_size(&b) != size) return 6;
  if(proof_serialize(wire2,size,&b) != 0) return 7;
  if(memcmp(wire1,wire2,size) != 0) return 8;

  wire1[4]++;
  if(proof_deserialize(&(proof){0},wire1,size) == 0) return 9;
  wire1[4]--;

  wire1[13] = (uint8_t)((wire1[13]+1)%32);
  if(proof_deserialize(&(proof){0},wire1,size) == 0) return 10;
  wire1[13] = wire2[13];

  memset(&wire1[120],0xff,4);
  if(proof_deserialize(&(proof){0},wire1,size) == 0) return 14;
  memcpy(&wire1[120],&wire2[120],4);

  wire3 = malloc(size+1);
  if(wire3 == NULL) return 11;
  memcpy(wire3,wire1,size);
  wire3[size] = 0;
  if(proof_deserialize(&(proof){0},wire3,size+1) == 0) return 12;
  if(proof_deserialize(&(proof){0},wire1,size-1) == 0) return 13;

  free_proof(&b);
  size = proof_contextual_serialized_size(&a);
  if(proof_serialize_contextual(wire1,size,&a) != 0) return 15;
  if(proof_deserialize_contextual(&b,&a,wire1,size,&consumed) != 0 ||
     consumed != size || b.foldnonce != a.foldnonce) return 16;
  if(proof_serialize_contextual(wire2,size,&b) != 0 ||
     memcmp(wire1,wire2,size) != 0) return 17;

  free(wire1);
  free(wire2);
  free(wire3);
  free_proof(&a);
  free_proof(&b);
  printf("Proof contextual wire round-trip: %zu bytes\n",size);
  return 0;
}
