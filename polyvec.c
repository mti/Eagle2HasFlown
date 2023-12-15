/*
 * Implementors: EagleSign Team
 * This implementation is highly inspired from Dilithium and
 * Falcon Signatures' implementations
 */

#include <stdint.h>
#include <stdio.h>
#include "params.h"
#include "polyvec.h"
#include "poly.h"
#include "fips202.h"

/**************************************************************/
/************ Vectors of polynomials of length L **************/
/**************************************************************/

void polyvecl_uniform_eta_f(polyvecl *v, const uint8_t seed[CRHBYTES], uint16_t *nonce)
{
  unsigned int i;

  for (i = 0; i < L; ++i)
    poly_uniform_eta(&v->vec[i], seed, (*nonce)++, 0);
}

void polyvecl_challenge(polyvecl *v, const uint8_t seed[SEEDBYTES], uint16_t *nonce, int param)
{
  unsigned int i;

  for (i = 0; i < L; ++i)
    poly_challenge(&v->vec[i], seed, (*nonce)++, param);
}

void polyvecl_add(polyvecl *w, const polyvecl *u, const polyvecl *v)
{
  unsigned int i;

  for (i = 0; i < L; ++i)
    poly_add(&w->vec[i], &u->vec[i], &v->vec[i]);
}

void polyvecl_ntt(polyvecl *v)
{
  unsigned int i;

  for (i = 0; i < L; ++i)
    poly_ntt(&v->vec[i]);
}

void polyvecl_invntt_tomont(polyvecl *v)
{
  unsigned int i;

  for (i = 0; i < L; ++i)
    poly_invntt_tomont(&v->vec[i]);
}

void polyvecl_pointwise_poly_montgomery(polyvecl *r, const poly *a, const polyvecl *v)
{
  unsigned int i;

  for (i = 0; i < L; ++i)
    poly_pointwise_montgomery(&r->vec[i], a, &v->vec[i]);
}

void polyvecl_pointwise_acc_montgomery(poly *w,
                                       const polyvecl *u,
                                       const polyvecl *v)
{
  unsigned int i;
  poly t;

  poly_pointwise_montgomery(w, &u->vec[0], &v->vec[0]);
  for (i = 1; i < L; ++i)
  {
    poly_pointwise_montgomery(&t, &u->vec[i], &v->vec[i]);
    poly_add(w, w, &t);
  }
}

/**************************************************************/
/************ Vectors of polynomials of length K **************/
/**************************************************************/

void polyveck_uniform_eta_y(polyvecl *v, const uint8_t seed[CRHBYTES], uint16_t *nonce)
{
  unsigned int i;

  for (i = 0; i < L; ++i)
    poly_uniform_eta(&v->vec[i], seed, (*nonce)++, 1);
}

void polyveck_add(polyveck *w, const polyveck *u, const polyveck *v)
{
  unsigned int i;

  for (i = 0; i < K; ++i)
    poly_add(&w->vec[i], &u->vec[i], &v->vec[i]);
}

void polyveck_sub(polyveck *w, const polyveck *u, const polyveck *v)
{
  unsigned int i;

  for (i = 0; i < K; ++i)
    poly_sub(&w->vec[i], &u->vec[i], &v->vec[i]);
}

void polyveck_ntt(polyveck *v)
{
  unsigned int i;

  for (i = 0; i < K; ++i)
    poly_ntt(&v->vec[i]);
}

void polyveck_invntt_tomont(polyveck *v)
{
  unsigned int i;

  for (i = 0; i < K; ++i)
    poly_invntt_tomont(&v->vec[i]);
}

void polyveck_pointwise_poly_montgomery(polyveck *r, const poly *a, const polyveck *v)
{
  unsigned int i;

  for (i = 0; i < K; ++i)
    poly_pointwise_montgomery(&r->vec[i], a, &v->vec[i]);
}

/*************************************************
 * Name:        polyveck_decompose
 *
 * Description: For all coefficients a of polynomials in vector of length K,
 *              compute high and low bits a0, a1 such a mod^+ Q = a1*ALPHA + a0
 *              with -ALPHA/2 < a0 <= ALPHA/2 except a1 = (Q-1)/ALPHA where we
 *              set a1 = 0 and -ALPHA/2 <= a0 = a mod Q - Q < 0.
 *              Assumes coefficients to be standard representatives.
 *
 * Arguments:   - polyveck *v1: pointer to output vector of polynomials with
 *                              coefficients a1
 *              - polyveck *v0: pointer to output vector of polynomials with
 *                              coefficients a0
 *              - const polyveck *v: pointer to input vector
 **************************************************/
void polyveck_decompose(polyveck *v1, polyveck *v0, const polyveck *v)
{
  unsigned int i;

  for (i = 0; i < K; ++i)
    poly_decompose(&v1->vec[i], &v0->vec[i], &v->vec[i]);
}

void polyveck_pack_P(uint8_t r[K * NBYTES * LOGQ], polyveck *P)
{
  unsigned int i;

  for (i = 0; i < K; ++i)
    poly_pack_ETA(r + i * NBYTES * LOGQ, &P->vec[i], LOGQ);
}

void polyveck_unpack_P(polyveck *P, uint8_t r[K * NBYTES * LOGQ])
{
  unsigned int i;

  for (i = 0; i < K; ++i)
    poly_unpack_ETA(&P->vec[i], r + i * NBYTES * LOGQ, LOGQ);
}

/*************************************************
 * Name:        polyvecl_chknorm
 *
 * Description: Check infinity norm of polynomials in vector of length L.
 *              Assumes input polyvecl to be reduced by polyvecl_reduce().
 *
 * Arguments:   - const polyvecl *v: pointer to vector
 *              - int32_t B: norm bound
 *
 * Returns 0 if norm of all polynomials is strictly smaller than B <= (Q-1)/8
 * and 1 otherwise.
 **************************************************/
int polyvecl_chknorm(const polyvecl *v, int32_t B)
{
  unsigned int i;

  for (i = 0; i < L; ++i)
    if (poly_chknorm(&v->vec[i], B))
      return -1;

  return 0;
}

/*************************************************
 * Name:        polyveck_chknorm
 *
 * Description: Check infinity norm of polynomials in vector of length K.
 *              Assumes input polyvecl to be reduced by polyvecl_reduce().
 *
 * Arguments:   - const polyvecl *v: pointer to vector
 *              - int32_t B: norm bound
 *
 * Returns 0 if norm of all polynomials is strictly smaller than B <= (Q-1)/8
 * and 1 otherwise.
 **************************************************/
int polyveck_chknorm(const polyveck *v, int32_t B)
{
  unsigned int i;

  for (i = 0; i < K; ++i)
    if (poly_chknorm(&v->vec[i], B))
      return -1;

  return 0;
}