#ifndef API_H
#define API_H

#include <stddef.h>
#include <stdint.h>

#define pq_eaglesign2_PUBLICKEYBYTES 2720
#define pq_eaglesign2_SECRETKEYBYTES 3520
#define pq_eaglesign2_BYTES 2336

#define pq_eaglesign2_ref_PUBLICKEYBYTES pq_eaglesign2_PUBLICKEYBYTES
#define pq_eaglesign2_ref_SECRETKEYBYTES pq_eaglesign2_SECRETKEYBYTES
#define pq_eaglesign2_ref_BYTES pq_eaglesign2_BYTES

int pq_eaglesign2_ref_keypair(uint8_t *pk, uint8_t *sk);

int pq_eaglesign2_ref_signature(uint8_t *sig, size_t *siglen,
                                const uint8_t *m, size_t mlen,
                                const uint8_t *sk);

int pq_eaglesign2_ref(uint8_t *sm, size_t *smlen,
                      const uint8_t *m, size_t mlen,
                      const uint8_t *sk);

int pq_eaglesign2_ref_verify(const uint8_t *sig, size_t siglen,
                             const uint8_t *m, size_t mlen,
                             const uint8_t *pk);

int pq_eaglesign2_ref_open(uint8_t *m, size_t *mlen,
                           const uint8_t *sm, size_t smlen,
                           const uint8_t *pk);

#define pq_eaglesign2aes_ref_PUBLICKEYBYTES pq_eaglesign2_ref_PUBLICKEYBYTES
#define pq_eaglesign2aes_ref_SECRETKEYBYTES pq_eaglesign2_ref_SECRETKEYBYTES
#define pq_eaglesign2aes_ref_BYTES pq_eaglesign2_ref_BYTES

int pq_eaglesign2aes_ref_keypair(uint8_t *pk, uint8_t *sk);

int pq_eaglesign2aes_ref_signature(uint8_t *sig, size_t *siglen,
                                   const uint8_t *m, size_t mlen,
                                   const uint8_t *sk);

int pq_eaglesign2aes_ref(uint8_t *sm, size_t *smlen,
                         const uint8_t *m, size_t mlen,
                         const uint8_t *sk);

int pq_eaglesign2aes_ref_verify(const uint8_t *sig, size_t siglen,
                                const uint8_t *m, size_t mlen,
                                const uint8_t *pk);

int pq_eaglesign2aes_ref_open(uint8_t *m, size_t *mlen,
                              const uint8_t *sm, size_t smlen,
                              const uint8_t *pk);


#define pq_eaglesign3aes_ref_PUBLICKEYBYTES pq_eaglesign3_ref_PUBLICKEYBYTES
#define pq_eaglesign3aes_ref_SECRETKEYBYTES pq_eaglesign3_ref_SECRETKEYBYTES
#define pq_eaglesign3aes_ref_BYTES pq_eaglesign3_ref_BYTES

int pq_eaglesign3aes_ref_keypair(uint8_t *pk, uint8_t *sk);

int pq_eaglesign3aes_ref_signature(uint8_t *sig, size_t *siglen,
                                   const uint8_t *m, size_t mlen,
                                   const uint8_t *sk);

int pq_eaglesign3aes_ref(uint8_t *sm, size_t *smlen,
                         const uint8_t *m, size_t mlen,
                         const uint8_t *sk);

int pq_eaglesign3aes_ref_verify(const uint8_t *sig, size_t siglen,
                                const uint8_t *m, size_t mlen,
                                const uint8_t *pk);

int pq_eaglesign3aes_ref_open(uint8_t *m, size_t *mlen,
                              const uint8_t *sm, size_t smlen,
                              const uint8_t *pk);

#define pq_eaglesign5_PUBLICKEYBYTES 6432
#define pq_eaglesign5_SECRETKEYBYTES 8000
#define pq_eaglesign5_BYTES 5408

#define pq_eaglesign5_ref_PUBLICKEYBYTES pq_eaglesign5_PUBLICKEYBYTES
#define pq_eaglesign5_ref_SECRETKEYBYTES pq_eaglesign5_SECRETKEYBYTES
#define pq_eaglesign5_ref_BYTES pq_eaglesign5_BYTES

int pq_eaglesign5_ref_keypair(uint8_t *pk, uint8_t *sk);

int pq_eaglesign5_ref_signature(uint8_t *sig, size_t *siglen,
                                const uint8_t *m, size_t mlen,
                                const uint8_t *sk);

int pq_eaglesign5_ref(uint8_t *sm, size_t *smlen,
                      const uint8_t *m, size_t mlen,
                      const uint8_t *sk);

int pq_eaglesign5_ref_verify(const uint8_t *sig, size_t siglen,
                             const uint8_t *m, size_t mlen,
                             const uint8_t *pk);

int pq_eaglesign5_ref_open(uint8_t *m, size_t *mlen,
                           const uint8_t *sm, size_t smlen,
                           const uint8_t *pk);

#define pq_eaglesign5aes_ref_PUBLICKEYBYTES pq_eaglesign5_ref_PUBLICKEYBYTES
#define pq_eaglesign5aes_ref_SECRETKEYBYTES pq_eaglesign5_ref_SECRETKEYBYTES
#define pq_eaglesign5aes_ref_BYTES pq_eaglesign5_ref_BYTES

int pq_eaglesign5aes_ref_keypair(uint8_t *pk, uint8_t *sk);

int pq_eaglesign5aes_ref_signature(uint8_t *sig, size_t *siglen,
                                   const uint8_t *m, size_t mlen,
                                   const uint8_t *sk);

int pq_eaglesign5aes_ref(uint8_t *sm, size_t *smlen,
                         const uint8_t *m, size_t mlen,
                         const uint8_t *sk);

int pq_eaglesign5aes_ref_verify(const uint8_t *sig, size_t siglen,
                                const uint8_t *m, size_t mlen,
                                const uint8_t *pk);

int pq_eaglesign5aes_ref_open(uint8_t *m, size_t *mlen,
                              const uint8_t *sm, size_t smlen,
                              const uint8_t *pk);

#define pq_eaglesign52_PUBLICKEYBYTES 17696
#define pq_eaglesign52_SECRETKEYBYTES 20544
#define pq_eaglesign52_BYTES 5408

#define pq_eaglesign52_ref_PUBLICKEYBYTES pq_eaglesign52_PUBLICKEYBYTES
#define pq_eaglesign52_ref_SECRETKEYBYTES pq_eaglesign52_SECRETKEYBYTES
#define pq_eaglesign52_ref_BYTES pq_eaglesign52_BYTES

int pq_eaglesign52_ref_keypair(uint8_t *pk, uint8_t *sk);

int pq_eaglesign52_ref_signature(uint8_t *sig, size_t *siglen,
                                 const uint8_t *m, size_t mlen,
                                 const uint8_t *sk);

int pq_eaglesign52_ref(uint8_t *sm, size_t *smlen,
                       const uint8_t *m, size_t mlen,
                       const uint8_t *sk);

int pq_eaglesign52_ref_verify(const uint8_t *sig, size_t siglen,
                              const uint8_t *m, size_t mlen,
                              const uint8_t *pk);

int pq_eaglesign52_ref_open(uint8_t *m, size_t *mlen,
                            const uint8_t *sm, size_t smlen,
                            const uint8_t *pk);

#define pq_eaglesign52aes_ref_PUBLICKEYBYTES pq_eaglesign52_ref_PUBLICKEYBYTES
#define pq_eaglesign52aes_ref_SECRETKEYBYTES pq_eaglesign52_ref_SECRETKEYBYTES
#define pq_eaglesign52aes_ref_BYTES pq_eaglesign52_ref_BYTES

int pq_eaglesign52aes_ref_keypair(uint8_t *pk, uint8_t *sk);

int pq_eaglesign52aes_ref_signature(uint8_t *sig, size_t *siglen,
                                    const uint8_t *m, size_t mlen,
                                    const uint8_t *sk);

int pq_eaglesign52aes_ref(uint8_t *sm, size_t *smlen,
                          const uint8_t *m, size_t mlen,
                          const uint8_t *sk);

int pq_eaglesign52aes_ref_verify(const uint8_t *sig, size_t siglen,
                                 const uint8_t *m, size_t mlen,
                                 const uint8_t *pk);

int pq_eaglesign52aes_ref_open(uint8_t *m, size_t *mlen,
                               const uint8_t *sm, size_t smlen,
                               const uint8_t *pk);

#endif
