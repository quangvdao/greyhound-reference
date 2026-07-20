CC ?= /usr/bin/cc
CFLAGS += -std=c2x -Wall -Wextra -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wpointer-arith -Wno-unused-function -flto=auto \
  -fwrapv -march=native -mtune=native -O3
RM = /bin/rm
OPENSSL_CFLAGS ?= $(shell pkg-config --cflags openssl 2>/dev/null)
OPENSSL_LDFLAGS ?= $(shell pkg-config --libs-only-L openssl 2>/dev/null)
GMP_CFLAGS ?= $(shell pkg-config --cflags gmp 2>/dev/null)
GMP_LDFLAGS ?= $(shell pkg-config --libs-only-L gmp 2>/dev/null)

ARCH := $(shell uname -m)
BACKEND ?= auto

ifneq ($(ARCH),x86_64)
CFLAGS += -Wno-pass-failed
endif

ifeq ($(BACKEND),portable)
NTT_SOURCES = ntt_portable.c
CFLAGS += -DGREYHOUND_PORTABLE_BACKEND
else ifeq ($(BACKEND),avx512)
ifneq ($(ARCH),x86_64)
$(error BACKEND=avx512 is only supported on x86-64)
endif
NTT_SOURCES = ntt.S invntt.S
else ifeq ($(BACKEND),auto)
ifeq ($(ARCH),x86_64)
NTT_SOURCES = ntt.S invntt.S
else
NTT_SOURCES = ntt_portable.c
CFLAGS += -DGREYHOUND_PORTABLE_BACKEND
endif
else
$(error BACKEND must be auto, portable, or avx512)
endif

SOURCES = pack.c pack_wire.c greyhound.c greyhound_wire.c dachshund.c chihuahua.c labrador.c proof_wire.c witness_wire.c rice.c \
  data.c jlproj.c polx.c poly.c polz.c sparsemat.c \
  aesctr.c fips202.c randombytes.c cpucycles.c parallel.c $(NTT_SOURCES)
HEADERS = pack.h greyhound.h dachshund.h chihuahua.h labrador.h \
  data.h jlproj.h polx.h poly.h polz.h sparsemat.h fq.inc shuffle.inc \
  aesctr.h fips202.h randombytes.h malloc.h cpucycles.h parallel.h rice.h

.PHONY: all clean

all: \
  test_aesctr \
  test_ntt \
  test_poly \
  test_polz \
  test_jlproj \
  test_sis_estimator \
  test_proof_wire \
  test_chihuahua \
  test_dachshund \
  test_greyhound

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test_aesctr: test_aesctr.c aesctr.c aesctr.h randombytes.c randombytes.h cpucycles.c cpucycles.h
	$(CC) $(CFLAGS) $(OPENSSL_CFLAGS) test_aesctr.c aesctr.c randombytes.c cpucycles.c -o test_aesctr $(OPENSSL_LDFLAGS) -lcrypto

test_ntt: test_ntt.c data.c data.h poly.c poly.h $(NTT_SOURCES) fq.inc shuffle.inc aesctr.c aesctr.h fips202.c fips202.h randombytes.c randombytes.h cpucycles.c cpucycles.h
	$(CC) $(CFLAGS) test_ntt.c data.c poly.c $(NTT_SOURCES) aesctr.c fips202.c randombytes.c cpucycles.c -o test_ntt -lm

test_poly: test_poly.c data.c data.h poly.c poly.h $(NTT_SOURCES) fq.inc shuffle.inc aesctr.c aesctr.h fips202.c fips202.h randombytes.c randombytes.h
	$(CC) $(CFLAGS) test_poly.c data.c poly.c $(NTT_SOURCES) aesctr.c fips202.c randombytes.c -o test_poly -lm

test_polz: test_polz.c data.c data.h polx.c polx.h poly.c poly.h polz.c polz.h $(NTT_SOURCES) fq.inc shuffle.inc aesctr.c aesctr.h fips202.c fips202.h randombytes.c randombytes.h cpucycles.c cpucycles.h
	$(CC) $(CFLAGS) $(GMP_CFLAGS) test_polz.c data.c polx.c poly.c polz.c parallel.c $(NTT_SOURCES) aesctr.c fips202.c randombytes.c cpucycles.c -o test_polz -lm $(GMP_LDFLAGS) -lgmp

test_jlproj: test_jlproj.c data.c data.h jlproj.c jlproj.h polx.c polx.h poly.c poly.h polz.c polz.h $(NTT_SOURCES) fq.inc shuffle.inc aesctr.c aesctr.h fips202.c fips202.h randombytes.c randombytes.h cpucycles.c cpucycles.h
	$(CC) $(CFLAGS) test_jlproj.c jlproj.c data.c polx.c poly.c polz.c parallel.c $(NTT_SOURCES) aesctr.c fips202.c randombytes.c cpucycles.c -o test_jlproj -lm

test_sis_estimator: test_sis_estimator.c $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) test_sis_estimator.c $(SOURCES) -o $@ -lm

test_proof_wire: test_proof_wire.c $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) test_proof_wire.c $(SOURCES) -o $@ -lm

test_chihuahua: test_chihuahua.c $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) test_chihuahua.c $(SOURCES) -o $@ -lm

test_dachshund: test_dachshund.c $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) test_dachshund.c $(SOURCES) -o $@ -lm

test_greyhound: test_greyhound.c $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) test_greyhound.c $(SOURCES) -o $@ -lm

libdogs.so: $(SOURCES) $(HEADERS)
	$(CC) -shared -fPIC -fvisibility=hidden $(CFLAGS) -o $@ $(SOURCES)

clean:
	-$(RM) -f *.o *.so
	-$(RM) -f test_aesctr
	-$(RM) -f test_ntt
	-$(RM) -f test_poly
	-$(RM) -f test_polz
	-$(RM) -f test_jlproj
	-$(RM) -f test_sis_estimator
	-$(RM) -f test_proof_wire
	-$(RM) -f test_chihuahua
	-$(RM) -f test_dachshund
	-$(RM) -f test_greyhound
