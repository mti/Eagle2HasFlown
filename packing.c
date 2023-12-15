/*
 * Implementors: EagleSign Team
 * This implementation is highly inspired from Dilithium and
 * Falcon Signatures' implementations
 */

#include "params.h"
#include "packing.h"
#include "polyvec.h"
#include "poly.h"

/*************************************************
 * Name:        pack_pk
 *
 * Description: Bit-pack public key pk = (rho, E).
 *
 * Arguments:   - uint8_t pk[]: output byte array
 *              - const uint8_t rho[]: byte array containing rho
 *              - const polyvecl E[K]: array containing the matrix E
 **************************************************/
void pack_pk(uint8_t pk[CRYPTO_EAGLESIGN_PUBLICKEYBYTES],
             const uint8_t rho[SEEDBYTES],
             const polyvecl E[K])
{
  unsigned int i, j;

  for (i = 0; i < SEEDBYTES; ++i)
    pk[i] = rho[i];
  pk += SEEDBYTES;

  for (i = 0; i < K; ++i)
    for (j = 0; j < L; ++j)
      poly_pack_ETA(pk + (i * L + j) * NBYTES * LOGQ, &E[i].vec[j], LOGQ);
}

/*************************************************
 * Name:        unpack_pk
 *
 * Description: Unpack public key pk = (rho, E).
 *
 * Arguments:   - uint8_t rho[]: output byte array for rho
 *              - polyvecl E[K]: array containing the matrix E
 *              - const uint8_t pk[]: byte array containing bit-packed pk
 **************************************************/
void unpack_pk(
    uint8_t rho[SEEDBYTES],
    polyvecl E[K],
    const uint8_t pk[CRYPTO_EAGLESIGN_PUBLICKEYBYTES])
{
  unsigned int i, j;

  for (i = 0; i < SEEDBYTES; ++i)
    rho[i] = pk[i];
  pk += SEEDBYTES;

  for (i = 0; i < K; ++i)
    for (j = 0; j < L; ++j)
      poly_unpack_ETA(&E[i].vec[j], pk + (i * L + j) * NBYTES * LOGQ, LOGQ);
}

/*************************************************
 * Name:        pack_sk
 *
 * Description: Bit-pack secret key sk = (rho, tr, g, D, F, Es).
 *
 * Arguments:   - uint8_t sk[]: output byte array
 *              - const uint8_t rho[]: byte array containing rho
 *              - const uint8_t tr[]: byte array containing tr
 *              - const poly g: the polynomial g
 *              - const polyvecl D[]: array containing the matrix D
 *              - const polyvecl F[]: array containing the matrix F
 *              - const polyvecl Es[]: array containing the matrix Es
 **************************************************/
void pack_sk(uint8_t sk[CRYPTO_EAGLESIGN_SECRETKEYBYTES],
             const uint8_t rho[SEEDBYTES],
             const uint8_t tr[SEEDBYTES],
             const poly g,
             const polyvecl D[K],
             const polyvecl F[L],
             const polyvecl Es[K])
{
  unsigned int i, j;

  for (i = 0; i < SEEDBYTES; ++i)
    sk[i] = rho[i];
  sk += SEEDBYTES;

  for (i = 0; i < SEEDBYTES; ++i)
    sk[i] = tr[i];
  sk += SEEDBYTES;

  poly_pack_ETA(sk, &g, 2);
  sk += NBYTES * 2;

  for (i = 0; i < K; ++i)
    for (j = 0; j < L; ++j)
      poly_pack_ETA(sk + (i * L + j) * NBYTES * 2, &D[i].vec[j], 2);

  sk += K * L * NBYTES * 2;

  for (i = 0; i < L; ++i)
    for (j = 0; j < L; ++j)
      poly_pack_ETA(sk + (i * L + j) * NBYTES * LOGETAF, &F[i].vec[j], LOGETAF);

  sk += L * L * NBYTES * LOGETAF;

  for (i = 0; i < K; ++i)
    for (j = 0; j < L; ++j)
      poly_pack_ETA(sk + (i * L + j) * NBYTES * LOGQ, &Es[i].vec[j], LOGQ);
}

/*************************************************
 * Name:        unpack_sk
 *
 * Description: Unpack secret key sk = (rho, tr, g, D, F, Es).
 *
 * Arguments:   - uint8_t rho[]: output byte array for rho
 *              - uint8_t tr[]: byte array containing tr
 *              - const poly g: the polynomial g
 *              - const polyvecl D[]: array containing the matrix D
 *              - const polyvecl F[]: array containing the matrix F
 *              - const polyvecl Es[]: array containing the matrix Es
 *              - const uint8_t sk[]: output byte array
 **************************************************/
void unpack_sk(
    uint8_t rho[SEEDBYTES],
    uint8_t tr[SEEDBYTES],
    poly *g,
    polyvecl D[K],
    polyvecl F[L],
    polyvecl Es[K],
    const uint8_t sk[CRYPTO_EAGLESIGN_SECRETKEYBYTES])
{
  unsigned int i, j;

  for (i = 0; i < SEEDBYTES; ++i)
    rho[i] = sk[i];
  sk += SEEDBYTES;

  for (i = 0; i < SEEDBYTES; ++i)
    tr[i] = sk[i];
  sk += SEEDBYTES;

  poly_unpack_ETA(g, sk, 2);
  sk += NBYTES * 2;

  for (i = 0; i < K; ++i)
    for (j = 0; j < L; ++j)
      poly_unpack_ETA(&D[i].vec[j], sk + (i * L + j) * NBYTES * 2, 2);

  sk += K * L * NBYTES * 2;

  for (i = 0; i < L; ++i)
    for (j = 0; j < L; ++j)
      poly_unpack_ETA(&F[i].vec[j], sk + (i * L + j) * NBYTES * LOGETAF, LOGETAF);

  sk += L * L * NBYTES * LOGETAF;

  for (i = 0; i < K; ++i)
    for (j = 0; j < L; ++j)
      poly_unpack_ETA(&Es[i].vec[j], sk + (i * L + j) * NBYTES * LOGQ, LOGQ);
}

/*************************************************
 * Name:        pack_sig
 *
 * Description: Bit-pack signature sig = (r, z).
 *
 * Arguments:   - uint8_t sig[]: output byte array
 *              - const uint8_t r[]: input byte array for r
 *              - const polyvecl *z: pointer to vector z
 **************************************************/
void pack_sig(uint8_t sig[CRYPTO_EAGLESIGN_BYTES],
              const uint8_t r[SEEDBYTES],
              const polyvecl *z)
{
  unsigned int i;

  for (i = 0; i < SEEDBYTES; ++i)
    sig[i] = r[i];
  sig += SEEDBYTES;

  for (i = 0; i < L; ++i)
    poly_pack_ETA(sig + i * NBYTES * Z_SIZE, &z->vec[i], Z_SIZE);
}

/*************************************************
 * Name:        unpack_sig
 *
 * Description: Unpack signature sig = (r, z).
 *
 * Arguments:   - uint8_t r[]: output byte array for r
 *              - polyvecl *z: pointer to vector z
 *              - const uint8_t sig[]: byte array containing
 *                bit-packed signature
 *
 * Returns 1 in case of malformed signature; otherwise 0.
 **************************************************/
int unpack_sig(
    uint8_t r[SEEDBYTES],
    polyvecl *z,
    const uint8_t sig[CRYPTO_EAGLESIGN_BYTES])
{
  unsigned int i;

  for (i = 0; i < SEEDBYTES; ++i)
    r[i] = sig[i];
  sig += SEEDBYTES;

  for (i = 0; i < L; ++i)
    poly_unpack_ETA(&z->vec[i], sig + i * NBYTES * Z_SIZE, Z_SIZE);
}
