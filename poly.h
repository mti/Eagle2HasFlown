/*
 * Implementors: EagleSign Team
 * This implementation is highly inspired from Dilithium and
 * Falcon Signatures' implementations
 */

#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "params.h"

typedef struct
{
  S_Q_SIZE coeffs[N];
} poly;

#define poly_add EAGLESIGN_NAMESPACE(poly_add)
void poly_add(poly *c, const poly *a, const poly *b);
#define poly_sub EAGLESIGN_NAMESPACE(poly_sub)
void poly_sub(poly *c, const poly *a, const poly *b);

#define poly_ntt EAGLESIGN_NAMESPACE(poly_ntt)
void poly_ntt(poly *a);
#define poly_invntt_tomont EAGLESIGN_NAMESPACE(poly_invntt_tomont)
void poly_invntt_tomont(poly *a);

#define poly_pointwise_montgomery EAGLESIGN_NAMESPACE(poly_pointwise_montgomery)
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b);

#define poly_uniform EAGLESIGN_NAMESPACE(poly_uniform)
void poly_uniform(poly *a,
                  const uint8_t seed[SEEDBYTES],
                  uint16_t nonce);
#define poly_uniform_eta EAGLESIGN_NAMESPACE(poly_uniform_eta)
void poly_uniform_eta(poly *a,
                      const uint8_t seed[CRHBYTES],
                      uint16_t nonce,
                      int param);

#define poly_challenge EAGLESIGN_NAMESPACE(poly_challenge)
void poly_challenge(poly *c, const uint8_t seed[SEEDBYTES],
                    uint16_t nonce, int param);

#define poly_decompose EAGLESIGN_NAMESPACE(_poly_decompose)
void poly_decompose(poly *a1, poly *a0, const poly *a);

#define poly_pack_ETA EAGLESIGN_NAMESPACE(poly_pack_ETA)
void poly_pack_ETA(uint8_t *r, const poly *a, unsigned int logeta);

#define poly_unpack_ETA EAGLESIGN_NAMESPACE(poly_unpack_ETA)
void poly_unpack_ETA(poly *r, const uint8_t *a, unsigned int logeta);

#define poly_chknorm EAGLESIGN_NAMESPACE(poly_chknorm)
int poly_chknorm(const poly *a, int32_t B);

#endif
