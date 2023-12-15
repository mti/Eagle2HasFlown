/*
 * Implementors: EagleSign Team
 * This implementation is highly inspired from Dilithium and
 * Falcon Signatures' implementations
 */

#include <stdint.h>
#include <stdio.h>
#include "params.h"
#include "sign.h"
#include "packing.h"
#include "polymatrix.h"
#include "polyvec.h"
#include "poly.h"
#include "randombytes.h"
#include "symmetric.h"
#include "fips202.h"

/*************************************************
 * Name:        crypto_sign_keypair
 *
 * Description: Generates public and private key.
 *
 * Arguments:   - uint8_t *pk: pointer to output public key (allocated
 *                             array of CRYPTO_EAGLESIGN_PUBLICKEYBYTES bytes)
 *              - uint8_t *sk: pointer to output private key (allocated
 *                             array of CRYPTO_EAGLESIGN_SECRETKEYBYTES bytes)
 *
 * Returns 0 (success)
 **************************************************/

int crypto_sign_keypair(uint8_t *pk, uint8_t *sk)
{
  uint8_t seedbuf[3 * SEEDBYTES + CRHBYTES];
  uint8_t tr[SEEDBYTES], *rhoprime1, *rhoprime3;
  const uint8_t *rho, *rhoprime2;
  poly g, g_inv;
  polyvecl A[K], F[L], D[K], F_INV[L], E[K], Es[K];
  polyveck tmp2[L];
  keccak_state state;

  /* Get randomness for rho and rhoprime */
  randombytes(seedbuf, SEEDBYTES);
  shake256(seedbuf, 3 * SEEDBYTES + CRHBYTES, seedbuf, SEEDBYTES);
  rho = seedbuf;
  rhoprime1 = seedbuf + SEEDBYTES;
  rhoprime2 = seedbuf + 2 * SEEDBYTES;
  rhoprime3 = seedbuf + 3 * SEEDBYTES;

  /* Expand matrix A in NTT form*/
  polymatrix_expand(A, rho);

  /* Small and Sparse polynomials based matrix D*/
  polymatrix_k_l_expand_d(D, rhoprime2);

  /* Small and Sparse polynomial g*/
retry_g:
  shake256(rhoprime1, SEEDBYTES, rhoprime1, SEEDBYTES);
  if (poly_expand_g_invertible(&g_inv, &g, rhoprime1))
    goto retry_g;

  /* Small and Uniform polynomials based invertible matrix F and its inverse F_INV */
retry_F:
  shake256(rhoprime3, CRHBYTES, rhoprime3, CRHBYTES);
  if (polymatrix_l_expand_f_invertible(F_INV, F, rhoprime3))
    goto retry_F;

  /* Compute Es = (AF^{-1} + D)*/
  polymatrix_pointwise_product(tmp2, A, F_INV);
  polymatrix_reformat(Es, tmp2);
  polymatrix_pointwise_add(Es, Es, D);

  /* Compute E = Es*g^{-1} */
  poly_pointwise_matrix_product(E, g_inv, Es);

  /*Bring back elements from ntt*/
  polymatrix_invntt_tomont_k_l(E);
  polymatrix_invntt_tomont_k_l(Es);
  poly_invntt_tomont(&g);
  polymatrix_invntt_tomont_l_l(F);
  polymatrix_invntt_tomont_k_l(D);

  /* Extracting public key pk = (rho, E) */
  pack_pk(pk, rho, E);

  /* Compute tr = H(pk)*/
  // tr = CRH1(pk)
  shake256(tr, SEEDBYTES, pk, CRYPTO_EAGLESIGN_PUBLICKEYBYTES);

  /* Extracting secret key sk = (rho, tr, g, D, F, Es) */
  pack_sk(sk, rho, tr, g, D, F, Es);
  return 0;
}

/*************************************************
 * Name:        crypto_sign_signature
 *
 * Description: Computes signature.
 *
 * Arguments:   - uint8_t *sig:   pointer to output signature (of length CRYPTO_EAGLESIGN_BYTES)
 *              - size_t *siglen: pointer to output length of signature
 *              - uint8_t *m:     pointer to message to be signed
 *              - size_t mlen:    length of message
 *              - uint8_t *sk:    pointer to bit-packed secret key
 *
 * Returns 0 (success)
 **************************************************/
int crypto_sign_signature(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen,
                          const uint8_t *sk)
{
  uint8_t tmp[K * NBYTES * LOGQ],
      rho[SEEDBYTES], tr[SEEDBYTES],
      mu[CRHBYTES], rhoprime[CRHBYTES],
      r[SEEDBYTES], cp[SEEDBYTES];
  polyvecl A[K], F[L], E[K], Es[K], D[K], y, c, u, z, h;
  polyveck p, p_highbits, p_lowbits, p_prime, p_prime_highbits, p_prime_lowbits, hp;
  poly g;
  uint16_t nonce = 0, nonce_c = 0;
  keccak_state state;

  unpack_sk(rho, tr, &g, D, F, Es, sk);

  poly_ntt(&g);
  polymatrix_ntt_k_l(Es);
  polymatrix_ntt_l_l(F);
  polymatrix_ntt_k_l(D);

  /* Expand matrix A in NTT form*/
  polymatrix_expand(A, rho);

  /* Generating ephemeral secret keys y */
  randombytes(rhoprime, CRHBYTES);
rej:
  polyveck_uniform_eta_y(&y, rhoprime, &nonce);

  /* Computing P = Es*y */
  polymatrix_pointwise_montgomery(&p, Es, &y);
  polyveck_invntt_tomont(&p);

  /* Decomposition p_highbits, p_lowbits = HighBits(p)*/
  polyveck_decompose(&p_highbits, &p_lowbits, &p);

  /* Compute r = G(p_highbits)*/
  polyveck_pack_P(tmp, &p_highbits);
  shake256_init(&state);
  shake256_absorb(&state, tmp, K * NBYTES * LOGQ);
  shake256_finalize(&state);
  shake256_squeeze(r, SEEDBYTES, &state);

  /* Compute mu = CRH(tr, m) */
  shake256_init(&state);
  shake256_absorb(&state, tr, SEEDBYTES);
  shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state);

  /* Call the random oracle and Compute c = H(mu,r)*/
  shake256_init(&state);
  shake256_absorb(&state, r, SEEDBYTES);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_finalize(&state);
  shake256_squeeze(cp, SEEDBYTES, &state);
  nonce_c = 0;
  polyvecl_challenge(&c, cp, &nonce_c, 2);

  /* Compute u = y + Fc */
  polymatrix_pointwise_montgomery_l_l(&h, F, &c);
  polyvecl_add(&u, &y, &h);

  /* Compute z = gu, reject if it reveals secret (Zero Knowledge Property)*/
  polyvecl_pointwise_poly_montgomery(&z, &g, &u);
  polyvecl_invntt_tomont(&z);
  if (polyvecl_chknorm(&z, TG * (GAMMA1 - BETA)))
    goto rej;

  /* Check that adding DFc does not change high bits of p and low bits
   * do not reveal secret information */
  polymatrix_pointwise_montgomery(&hp, D, &h);
  polyveck_invntt_tomont(&hp);
  polyveck_add(&p_prime, &p, &hp);
  if (polyveck_chknorm(&p_prime, (Q >> 1) - L * TD * BETA))
    goto rej;
  polyveck_decompose(&p_prime_highbits, &p_prime_lowbits, &p_prime);
  if (polyveck_chknorm(&p_prime_lowbits, GAMMA2 - L * TD * BETA))
    goto rej;

  /* Packing Signature (r, z)*/
  pack_sig(sig, r, &z);

  *siglen = CRYPTO_EAGLESIGN_BYTES;
  return 0;
}

/*************************************************
 * Name:        crypto_sign
 *
 * Description: Compute signed message.
 *
 * Arguments:   - uint8_t *sm: pointer to output signed message (allocated
 *                             array with CRYPTO_EAGLESIGN_BYTES + mlen bytes),
 *                             can be equal to m
 *              - size_t *smlen: pointer to output length of signed
 *                               message
 *              - const uint8_t *m: pointer to message to be signed
 *              - size_t mlen: length of message
 *              - const uint8_t *sk: pointer to bit-packed secret key
 *
 * Returns 0 (success)
 **************************************************/
int crypto_sign(uint8_t *sm,
                size_t *smlen,
                const uint8_t *m,
                size_t mlen,
                const uint8_t *sk)
{

  size_t i;

  for (i = 0; i < mlen; ++i)
    sm[CRYPTO_EAGLESIGN_BYTES + mlen - 1 - i] = m[mlen - 1 - i];
  crypto_sign_signature(sm, smlen, sm + CRYPTO_EAGLESIGN_BYTES, mlen, sk);

  *smlen += mlen;
  return 0;
}

/*************************************************
 * Name:        crypto_sign_verify
 *
 * Description: Verifies signature.
 *
 * Arguments:   - uint8_t *m: pointer to input signature
 *              - size_t siglen: length of signature
 *              - const uint8_t *m: pointer to message
 *              - size_t mlen: length of message
 *              - const uint8_t *pk: pointer to bit-packed public key
 *
 * Returns 0 if signature could be verified correctly and -1 otherwise
 **************************************************/
int crypto_sign_verify(const uint8_t *sig,
                       size_t siglen,
                       const uint8_t *m,
                       size_t mlen,
                       const uint8_t *pk)
{
  unsigned int i, j;
  uint8_t rho[SEEDBYTES], tmp[K * NBYTES * LOGQ],
      mu[CRHBYTES], cp[SEEDBYTES], r[SEEDBYTES],
      r_prime[SEEDBYTES];
  polyvecl A[K], E[K], c, c_prime, z;
  polyveck v, v_highbits, v_lowbits, v1;
  keccak_state state;
  uint16_t nonce_c = 0;

  unpack_pk(rho, E, pk);
  unpack_sig(r, &z, sig);

  if (polyvecl_chknorm(&z, TG * (GAMMA1 - BETA)))
    return -1;

  /* Applying NTT Transformation*/
  polymatrix_ntt_k_l(E);
  polyvecl_ntt(&z);

  /* Compute mu = CRH(H(pk), msg) */
  shake256(mu, SEEDBYTES, pk, CRYPTO_EAGLESIGN_PUBLICKEYBYTES);
  shake256_init(&state);
  shake256_absorb(&state, mu, SEEDBYTES);
  shake256_absorb(&state, m, mlen);
  shake256_finalize(&state);
  shake256_squeeze(mu, CRHBYTES, &state);

  /* Call the random oracle and Compute C = H(mu,r) dans B_tau^l*/
  shake256_init(&state);
  shake256_absorb(&state, r, SEEDBYTES);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_finalize(&state);
  shake256_squeeze(cp, SEEDBYTES, &state);
  nonce_c = 0;
  polyvecl_challenge(&c, cp, &nonce_c, 2);

  /* Expand matrix A in NTT form*/
  polymatrix_expand(A, rho);

  /* Compute V = Ez − Ac */
  polymatrix_pointwise_montgomery(&v, E, &z);
  polymatrix_pointwise_montgomery(&v1, A, &c);
  polyveck_sub(&v, &v, &v1);
  polyveck_invntt_tomont(&v);

  /* Decomposition v_highbits, v_lowbits = HighBits(v)*/
  polyveck_decompose(&v_highbits, &v_lowbits, &v);

  /* Compute r_prime = G(V) */
  polyveck_pack_P(tmp, &v_highbits);

  shake256_init(&state);
  shake256_absorb(&state, tmp, K * NBYTES * LOGQ);
  shake256_finalize(&state);
  shake256_squeeze(r_prime, SEEDBYTES, &state);

  /* Compute c_prime = H(mu,r_prime) */
  shake256_init(&state);
  shake256_absorb(&state, r_prime, SEEDBYTES);
  shake256_absorb(&state, mu, CRHBYTES);
  shake256_finalize(&state);
  shake256_squeeze(cp, SEEDBYTES, &state);
  nonce_c = 0;
  polyvecl_challenge(&c_prime, cp, &nonce_c, 2);

  /* Comparing c and c_prime */
  for (i = 0; i < L; ++i)
    for (j = 0; j < N; ++j)
      if (c.vec[i].coeffs[j] != c_prime.vec[i].coeffs[j])
        return -1;

  return 0;
}

/*************************************************
 * Name:        crypto_sign_open
 *
 * Description: Verify signed message.
 *
 * Arguments:   - uint8_t *m: pointer to output message (allocated
 *                            array with smlen bytes), can be equal to sm
 *              - size_t *mlen: pointer to output length of message
 *              - const uint8_t *sm: pointer to signed message
 *              - size_t smlen: length of signed message
 *              - const uint8_t *pk: pointer to bit-packed public key
 *
 * Returns 0 if signed message could be verified correctly and -1 otherwise
 **************************************************/
int crypto_sign_open(uint8_t *m,
                     size_t *mlen,
                     const uint8_t *sm,
                     size_t smlen,
                     const uint8_t *pk)
{

  size_t i;

  if (smlen < CRYPTO_EAGLESIGN_BYTES)
    goto badsig;

  *mlen = smlen - CRYPTO_EAGLESIGN_BYTES;

  if (crypto_sign_verify(sm, CRYPTO_EAGLESIGN_BYTES, sm + CRYPTO_EAGLESIGN_BYTES, *mlen, pk))
    goto badsig;
  else
  {
    /* All good, copy msg, return 0 */
    for (i = 0; i < *mlen; ++i)
      m[i] = sm[CRYPTO_EAGLESIGN_BYTES + i];
    return 0;
  }

badsig:
  /* Signature verification failed */
  *mlen = -1;
  for (i = 0; i < smlen; ++i)
    m[i] = 0;

  return -1;
}
