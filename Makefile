CC ?= /usr/bin/cc
CFLAGS += -Ofast -march=native -fopenmp
NISTFLAGS += 
SOURCES = sign.c packing.c polyvec.c poly.c ntt.c reduce.c polymatrix.c
HEADERS = config.h params.h api.h sign.h packing.h polymatrix.h polyvec.h poly.h ntt.h \
  reduce.h symmetric.h randombytes.h
KECCAK_SOURCES = $(SOURCES) fips202.c symmetric-shake.c
KECCAK_HEADERS = $(HEADERS) fips202.h

.PHONY: all speed shared clean

all: attack

kats: \
  PQgenKAT_sign2 \

tests: \
  test/test_eaglesign2 \

test_vectors: \
  test/test_vectors2 \

speed: \
  test/test_speed2 \

shared: \
  libpq_eaglesign2_ref.so \
  libpq_fips202_ref.so \
  libpq_aes256ctr_ref.so \

attack: \
  test/key_recovery2 test/key_recoveryF2 \

libpq_fips202_ref.so: fips202.c fips202.h
	$(CC) -shared -fPIC $(CFLAGS) -o $@ $<

libpq_aes256ctr_ref.so: aes256ctr.c aes256ctr.h
	$(CC) -shared -fPIC $(CFLAGS) -o $@ $<

libpq_eaglesign2_ref.so: $(SOURCES) $(HEADERS) symmetric-shake.c
	$(CC) -shared -fPIC $(CFLAGS) -DEAGLESIGN_MODE=2 \
	  -o $@ $(SOURCES) symmetric-shake.c

libpq_eaglesign5_ref.so: $(SOURCES) $(HEADERS) symmetric-shake.c
	$(CC) -shared -fPIC $(CFLAGS) -DEAGLESIGN_MODE=5 \
	  -o $@ $(SOURCES) symmetric-shake.c

libpq_eaglesign52_ref.so: $(SOURCES) $(HEADERS) symmetric-shake.c
	$(CC) -shared -fPIC $(CFLAGS) -DEAGLESIGN_MODE=52 \
	  -o $@ $(SOURCES) symmetric-shake.c

test/test_eaglesign2: test/test_eaglesign.c rng.c  $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DEAGLESIGN_MODE=2 \
	  -o $@ $< rng.c $(KECCAK_SOURCES) -lcrypto

test/test_eaglesign5: test/test_eaglesign.c rng.c  $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DEAGLESIGN_MODE=5 \
	  -o $@ $< rng.c $(KECCAK_SOURCES) -lcrypto

test/test_eaglesign52: test/test_eaglesign.c rng.c  $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DEAGLESIGN_MODE=52 \
	  -o $@ $< rng.c $(KECCAK_SOURCES) -lcrypto

test/test_vectors2: test/test_vectors.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DEAGLESIGN_MODE=2 \
	  -o $@ $< $(KECCAK_SOURCES)

test/test_vectors5: test/test_vectors.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DEAGLESIGN_MODE=5 \
	  -o $@ $< $(KECCAK_SOURCES)

test/test_vectors52: test/test_vectors.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DEAGLESIGN_MODE=52 \
	  -o $@ $< $(KECCAK_SOURCES)

test/test_speed2: test/test_speed.c test/speed_print.c test/speed_print.h \
  test/cpucycles.c test/cpucycles.h randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DEAGLESIGN_MODE=2 \
	  -o $@ $< test/speed_print.c test/cpucycles.c randombytes.c \
	  $(KECCAK_SOURCES)

test/test_speed5: test/test_speed.c test/speed_print.c test/speed_print.h \
  test/cpucycles.c test/cpucycles.h randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DEAGLESIGN_MODE=5 \
	  -o $@ $< test/speed_print.c test/cpucycles.c randombytes.c \
	  $(KECCAK_SOURCES)

test/test_speed52: test/test_speed.c test/speed_print.c test/speed_print.h \
  test/cpucycles.c test/cpucycles.h randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DEAGLESIGN_MODE=52 \
	  -o $@ $< test/speed_print.c test/cpucycles.c randombytes.c \
	  $(KECCAK_SOURCES)

PQgenKAT_sign2: PQgenKAT_sign.c rng.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(NISTFLAGS) -DEAGLESIGN_MODE=2 \
	  -o $@ $< rng.c $(KECCAK_SOURCES) $(LDFLAGS) -lcrypto

PQgenKAT_sign5: PQgenKAT_sign.c rng.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(NISTFLAGS) -DEAGLESIGN_MODE=5 \
	  -o $@ $< rng.c $(KECCAK_SOURCES) $(LDFLAGS) -lcrypto

PQgenKAT_sign52: PQgenKAT_sign.c rng.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(NISTFLAGS) -DEAGLESIGN_MODE=52 \
	  -o $@ $< rng.c $(KECCAK_SOURCES) $(LDFLAGS) -lcrypto

test/key_recovery2: test/key_recovery.c test/speed_print.c test/speed_print.h \
  test/cpucycles.c test/cpucycles.h randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DEAGLESIGN_MODE=2 \
	  -o $@ $< test/speed_print.c test/cpucycles.c randombytes.c \
	  $(KECCAK_SOURCES) -lm

test/key_recovery5: test/key_recovery.c test/speed_print.c test/speed_print.h \
  test/cpucycles.c test/cpucycles.h randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DEAGLESIGN_MODE=5 \
	  -o $@ $< test/speed_print.c test/cpucycles.c randombytes.c \
	  $(KECCAK_SOURCES) -lm

test/key_recovery52: test/key_recovery.c test/speed_print.c test/speed_print.h \
  test/cpucycles.c test/cpucycles.h randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DEAGLESIGN_MODE=52 \
	  -o $@ $< test/speed_print.c test/cpucycles.c randombytes.c \
	  $(KECCAK_SOURCES) -lm

test/key_recoveryF2: test/key_recoveryF.c test/speed_print.c test/speed_print.h \
  test/cpucycles.c test/cpucycles.h randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DEAGLESIGN_MODE=2 \
	  -o $@ $< test/speed_print.c test/cpucycles.c randombytes.c \
	  $(KECCAK_SOURCES) -lm

test/key_recoveryF5: test/key_recoveryF.c test/speed_print.c test/speed_print.h \
  test/cpucycles.c test/cpucycles.h randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DEAGLESIGN_MODE=5 \
	  -o $@ $< test/speed_print.c test/cpucycles.c randombytes.c \
	  $(KECCAK_SOURCES) -lm

test/key_recoveryF52: test/key_recoveryF.c test/speed_print.c test/speed_print.h \
  test/cpucycles.c test/cpucycles.h randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DEAGLESIGN_MODE=52 \
	  -o $@ $< test/speed_print.c test/cpucycles.c randombytes.c \
	  $(KECCAK_SOURCES) -lm


clean:
	rm -f *~ test/*~ *.gcno *.gcda *.lcov
	rm -f libpq_eaglesign2_ref.so
	rm -f libpq_eaglesign5_ref.so
	rm -f libpq_eaglesign52_ref.so
	rm -f libpq_fips202_ref.so
	rm -f libpq_aes256ctr_ref.so
	rm -f test/test_eaglesign2
	rm -f test/test_eaglesign5
	rm -f test/test_eaglesign52
	rm -f test/test_vectors2
	rm -f test/test_vectors5
	rm -f test/test_vectors52
	rm -f test/test_speed2
	rm -f test/test_speed5
	rm -f test/test_speed52
	rm -f PQgenKAT_sign2
	rm -f PQgenKAT_sign5
	rm -f PQgenKAT_sign52
	rm -f test/key_recovery2
	rm -f test/key_recovery5
	rm -f test/key_recovery52
	rm -f test/key_recoveryF2
	rm -f test/key_recoveryF5
	rm -f test/key_recoveryF52
