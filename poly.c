/*
 * Implementors: EagleSign Team
 * This implementation is highly inspired from Dilithium and
 * Falcon Signatures' implementations
 */

#include <stdint.h>
#include <stdio.h>
#include "params.h"
#include "poly.h"
#include "ntt.h"
#include "reduce.h"
#include "symmetric.h"
#include "polymatrix.h"

void poly_add(poly *c, const poly *a, const poly *b)
{
  unsigned int i;

  for (i = 0; i < N; ++i)
    c->coeffs[i] = (S_Q_SIZE)add(a->coeffs[i], b->coeffs[i]);
}

void poly_sub(poly *c, const poly *a, const poly *b)
{
  unsigned int i;

  for (i = 0; i < N; ++i)
    c->coeffs[i] = (S_Q_SIZE)sub(a->coeffs[i], b->coeffs[i]);
}

void poly_ntt(poly *a)
{

  ntt(a->coeffs, LOG_N);
}

void poly_invntt_tomont(poly *a)
{

  invntt_tomont(a->coeffs, LOG_N);
}

void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b)
{
  unsigned int i;

  for (i = 0; i < N; ++i)
    c->coeffs[i] = (S_Q_SIZE)montymul2(a->coeffs[i], b->coeffs[i]);
}

/*************************************************
 * Name:        poly_decompose
 *
 * Description: For all coefficients c of the input polynomial,
 *              compute high and low bits c0, c1 such c mod Q = c1*ALPHA + c0
 *              with -ALPHA/2 < c0 <= ALPHA/2 except c1 = (Q-1)/ALPHA where we
 *              set c1 = 0 and -ALPHA/2 <= c0 = c mod Q - Q < 0.
 *              Assumes coefficients to be standard representatives.
 *
 * Arguments:   - poly *a1: pointer to output polynomial with coefficients c1
 *              - poly *a0: pointer to output polynomial with coefficients c0
 *              - const poly *a: pointer to input polynomial
 **************************************************/
void poly_decompose(poly *a1, poly *a0, const poly *a)
{
  unsigned int i;

  for (i = 0; i < N; ++i)
    a1->coeffs[i] = decompose(&a0->coeffs[i], a->coeffs[i]);
}

/*************************************************
 * Name:        poly_chknorm
 *
 * Description: Check infinity norm of polynomial against given bound.
 *
 * Arguments:   - const poly *a: pointer to polynomial
 *              - int32_t B: norm bound
 *
 * Returns 0 if norm is strictly smaller than B <= (Q-1)/8 and 1 otherwise.
 **************************************************/
int poly_chknorm(const poly *a, int32_t B)
{
  unsigned int i;
  int32_t t;

  for (i = 0; i < N; ++i)
  {
    /* Absolute value */
    t = a->coeffs[i] >> (Q_BIT_SIZE - 1);
    t = a->coeffs[i] - (t & 2 * a->coeffs[i]);

    if (t >= B)
      return -1;
  }

  return 0;
}

static unsigned int rej_uniform(S_Q_SIZE *a,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  unsigned int ctr, pos;
  Q_SIZE t;

  ctr = pos = 0;
  while (ctr < len && pos + 2 <= buflen)
  {
    t = buf[pos++];
    t |= (Q_SIZE)buf[pos++] << 8;
    t &= 0x3FFF;

    if (t < Q)
      a[ctr++] = t;
  }

  return ctr;
}

#define POLY_UNIFORM_NBLOCKS ((768 + STREAM128_BLOCKBYTES - 1) / STREAM128_BLOCKBYTES)
void poly_uniform(poly *a,
                  const uint8_t seed[SEEDBYTES],
                  uint16_t nonce)
{
  unsigned int i, ctr, off;
  unsigned int buflen = POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES;
  uint8_t buf[POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES + 1];
  stream128_state state;

  stream128_init(&state, seed, nonce);
  stream128_squeezeblocks(buf, POLY_UNIFORM_NBLOCKS, &state);

  ctr = rej_uniform(a->coeffs, N, buf, buflen);

  while (ctr < N)
  {
    off = buflen % 2;
    for (i = 0; i < off; ++i)
      buf[i] = buf[buflen - off + i];

    stream128_squeezeblocks(buf + off, 1, &state);
    buflen = STREAM128_BLOCKBYTES + off;
    ctr += rej_uniform(a->coeffs + ctr, N - ctr, buf, buflen);
  }
}

static unsigned int rej_challenge(S_Q_SIZE *c,
                                  unsigned int *ctr,
                                  const uint8_t *buf,
                                  unsigned int buflen,
                                  uint64_t signs,
                                  unsigned int *h,
                                  int T_param)
{
  Q_SIZE t;
  unsigned int pos = 0;
  uint64_t v;

  while ((*ctr) < N && (*h) < T_param && pos + 3 <= buflen && signs != 0)
  {
    t = buf[pos++];
    t |= (Q_SIZE)buf[pos++] << 8;
    t &= N - 1;

    if (t < (*ctr))
    {
      c[*ctr] = c[t];

      v = 1 - 2 * (signs & 1);
      v += ((-Q) & ~-((v - (Q >> 1)) >> 63)) | (Q & -((v + (Q >> 1)) >> 63));

      c[t] = (S_Q_SIZE)v;
      signs >>= 1;
      *ctr += 1;
      *h += 1;
    }
  }
}

void poly_challenge(poly *c, const uint8_t seed[SEEDBYTES],
                    uint16_t nonce, int param)
{
  unsigned int i, ctr, off, b, pos, h = 0, T_param;
  unsigned int buflen = POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES;
  uint8_t buf[POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES + 2];
  stream128_state state;

  stream128_init(&state, seed, nonce);
  stream128_squeezeblocks(buf, POLY_UNIFORM_NBLOCKS, &state);

  uint64_t signs;

  for (i = 0; i < N; ++i)
    c->coeffs[i] = 0;

  T_param = ((2 - param) * (1 - param) * TG + 2 * param * (2 - param) * TD - param * (1 - param) * TAU) >> 1;
  ctr = N - T_param;

  while (ctr < N)
  {
    off = buflen % 3;
    for (i = 0; i < off; ++i)
      buf[i] = buf[buflen - off + i];

    stream128_squeezeblocks(buf + off, 1, &state);
    buflen = STREAM128_BLOCKBYTES + off;

    signs = 0;
    for (i = 0; i < 8; ++i)
      signs |= (uint64_t)buf[i] << 8 * i;

    rej_challenge(c->coeffs, &ctr, buf, buflen, signs, &h, T_param);
  }
  poly_ntt(c);
}

static unsigned int rej_eta(S_Q_SIZE *a,
                            unsigned int len,
                            const uint8_t *buf,
                            unsigned int buflen,
                            unsigned int eta)
{
  unsigned int ctr, pos;
  Q_SIZE t0, t1, t2, t3;

  ctr = pos = 0;

  while (ctr < len && pos < buflen - 4)
  {
    if (eta == 1)
    {
      t0 = buf[pos] & 0x03;
      t1 = (buf[pos] >> 2) & 0x03;
      t2 = (buf[pos] >> 4) & 0x03;
      t3 = buf[pos++] >> 6;
      if (t0 < 3)
        a[ctr++] = 1 - t0;
      if (t1 < 3 && ctr < len)
        a[ctr++] = 1 - t1;
      if (t2 < 3 && ctr < len)
        a[ctr++] = 1 - t2;
      if (t3 < 3 && ctr < len)
        a[ctr++] = 1 - t3;
    }

    else if (eta == 2)
    {
      t0 = buf[pos] & 0x07;
      t1 = buf[pos++] >> 5;
      if (t0 < 5)
        a[ctr++] = 2 - t0;
      if (t1 < 5 && ctr < len)
        a[ctr++] = 2 - t1;
    }

    else if (eta == 3)
    {
      t0 = buf[pos] & 0x07;
      t1 = buf[pos++] >> 5;
      if (t0 < 7)
        a[ctr++] = 3 - t0;
      if (t1 < 7 && ctr < len)
        a[ctr++] = 3 - t1;
    }

    else if (eta == 32)
    {
      t0 = buf[pos++] & 0x7f;
      if (t0 < 65)
        a[ctr++] = 32 - t0;
    }

    else if (eta == 64)
    {
      t0 = buf[pos++];
      if (t0 < 129)
        a[ctr++] = 64 - t0;
    }

    else if (eta == GAMMA1)
    {
      t0 = 0;
      for (int al = 0; al < (LOGGAMMA1 + 2); al = al + 8)
      {
        t0 |= buf[pos++] << al;
      }
      t0 &= 4 * GAMMA1 - 1;
      if (t0 < (2 * GAMMA1 + 1))
        a[ctr++] = GAMMA1 - t0;
    }
  }

  return ctr;
}

void poly_uniform_eta(poly *a,
                      const uint8_t seed[CRHBYTES],
                      uint16_t nonce, int param)
{
  unsigned int ctr, POLY_UNIFORM_ETA_NBLOCKS = ((227 + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES);
  unsigned int buflen = POLY_UNIFORM_ETA_NBLOCKS * STREAM256_BLOCKBYTES;
  uint8_t buf[POLY_UNIFORM_ETA_NBLOCKS * STREAM256_BLOCKBYTES];
  stream256_state state;
  int eta_param;

  stream256_init(&state, seed, nonce);
  stream256_squeezeblocks(buf, POLY_UNIFORM_ETA_NBLOCKS, &state);

  eta_param = GAMMA1 * param + ETAF * (1 - param);

  ctr = rej_eta(a->coeffs, N, buf, buflen, eta_param);

  while (ctr < N)
  {
    stream256_squeezeblocks(buf, 1, &state);
    ctr += rej_eta(a->coeffs + ctr, N - ctr, buf, STREAM256_BLOCKBYTES, eta_param);
  }
  poly_ntt(a);
}

void poly_pack_ETA(uint8_t *r, const poly *a, unsigned int logeta)
{
  unsigned int i;

  if (logeta == 2)
  {
    Q_SIZE t[4];
    for (i = 0; i < N / 4; ++i)
    {
      t[0] = (1 << (logeta - 1)) - a->coeffs[4 * i + 0];
      t[1] = (1 << (logeta - 1)) - a->coeffs[4 * i + 1];
      t[2] = (1 << (logeta - 1)) - a->coeffs[4 * i + 2];
      t[3] = (1 << (logeta - 1)) - a->coeffs[4 * i + 3];

      r[i] = t[0];
      r[i] |= t[1] << 2;
      r[i] |= t[2] << 4;
      r[i] |= t[3] << 6;
    }
  }

  else if (logeta == 3)
  {
    Q_SIZE t[8];
    for (i = 0; i < N / 8; ++i)
    {
      t[0] = (1 << (logeta - 1)) - a->coeffs[8 * i + 0];
      t[1] = (1 << (logeta - 1)) - a->coeffs[8 * i + 1];
      t[2] = (1 << (logeta - 1)) - a->coeffs[8 * i + 2];
      t[3] = (1 << (logeta - 1)) - a->coeffs[8 * i + 3];
      t[4] = (1 << (logeta - 1)) - a->coeffs[8 * i + 4];
      t[5] = (1 << (logeta - 1)) - a->coeffs[8 * i + 5];
      t[6] = (1 << (logeta - 1)) - a->coeffs[8 * i + 6];
      t[7] = (1 << (logeta - 1)) - a->coeffs[8 * i + 7];

      r[3 * i + 0] = t[0];
      r[3 * i + 0] |= t[1] << 3;
      r[3 * i + 0] |= t[2] << 6;
      r[3 * i + 1] = t[2] >> 2;
      r[3 * i + 1] |= t[3] << 1;
      r[3 * i + 1] |= t[4] << 4;
      r[3 * i + 1] |= t[5] << 7;
      r[3 * i + 2] = t[5] >> 1;
      r[3 * i + 2] |= t[6] << 2;
      r[3 * i + 2] |= t[7] << 5;
    }
  }
  else if (logeta == 18)
  {
    Q_SIZE t[4];
    for (i = 0; i < N / 4; ++i)
    {

      t[0] = (1 << (logeta - 1)) - a->coeffs[4 * i + 0];
      t[1] = (1 << (logeta - 1)) - a->coeffs[4 * i + 1];
      t[2] = (1 << (logeta - 1)) - a->coeffs[4 * i + 2];
      t[3] = (1 << (logeta - 1)) - a->coeffs[4 * i + 3];

      r[9 * i + 0] = t[0];
      r[9 * i + 1] = t[0] >> 8;
      r[9 * i + 2] = t[0] >> 16;
      r[9 * i + 2] |= t[1] << 2;
      r[9 * i + 3] = t[1] >> 6;
      r[9 * i + 4] = t[1] >> 14;
      r[9 * i + 4] |= t[2] << 4;
      r[9 * i + 5] = t[2] >> 4;
      r[9 * i + 6] = t[2] >> 12;
      r[9 * i + 6] |= t[3] << 6;
      r[9 * i + 7] = t[3] >> 2;
      r[9 * i + 8] = t[3] >> 10;
    }
  }
  else if (logeta == 21)
  {
    Q_SIZE t[8];
    for (i = 0; i < N / 8; ++i)
    {

      t[0] = (1 << (logeta - 1)) - a->coeffs[8 * i + 0];
      t[1] = (1 << (logeta - 1)) - a->coeffs[8 * i + 1];
      t[2] = (1 << (logeta - 1)) - a->coeffs[8 * i + 2];
      t[3] = (1 << (logeta - 1)) - a->coeffs[8 * i + 3];
      t[4] = (1 << (logeta - 1)) - a->coeffs[8 * i + 4];
      t[5] = (1 << (logeta - 1)) - a->coeffs[8 * i + 5];
      t[6] = (1 << (logeta - 1)) - a->coeffs[8 * i + 6];
      t[7] = (1 << (logeta - 1)) - a->coeffs[8 * i + 7];

      r[21 * i + 0] = t[0];
      r[21 * i + 1] = t[0] >> 8;
      r[21 * i + 2] = t[0] >> 16;
      r[21 * i + 2] |= t[1] << 5;
      r[21 * i + 3] = t[1] >> 3;
      r[21 * i + 4] = t[1] >> 11;
      r[21 * i + 5] = t[1] >> 19;
      r[21 * i + 5] |= t[2] << 2;
      r[21 * i + 6] = t[2] >> 6;
      r[21 * i + 7] = t[2] >> 14;
      r[21 * i + 7] |= t[3] << 7;
      r[21 * i + 8] = t[3] >> 1;
      r[21 * i + 9] = t[3] >> 9;
      r[21 * i + 10] = t[3] >> 17;
      r[21 * i + 10] |= t[4] << 4;
      r[21 * i + 11] = t[4] >> 4;
      r[21 * i + 12] = t[4] >> 12;
      r[21 * i + 13] = t[4] >> 20;
      r[21 * i + 13] |= t[5] << 1;
      r[21 * i + 14] = t[5] >> 7;
      r[21 * i + 15] = t[5] >> 15;
      r[21 * i + 15] |= t[6] << 6;
      r[21 * i + 16] = t[6] >> 2;
      r[21 * i + 17] = t[6] >> 10;
      r[21 * i + 18] = t[6] >> 18;
      r[21 * i + 18] |= t[7] << 3;
      r[21 * i + 19] = t[7] >> 5;
      r[21 * i + 20] = t[7] >> 13;
    }
  }
  else if (logeta == 22)
  {
    Q_SIZE t[4];
    for (i = 0; i < N / 4; ++i)
    {

      t[0] = (1 << (logeta - 1)) - a->coeffs[4 * i + 0];
      t[1] = (1 << (logeta - 1)) - a->coeffs[4 * i + 1];
      t[2] = (1 << (logeta - 1)) - a->coeffs[4 * i + 2];
      t[3] = (1 << (logeta - 1)) - a->coeffs[4 * i + 3];

      r[11 * i + 0] = t[0];
      r[11 * i + 1] = t[0] >> 8;
      r[11 * i + 2] = t[0] >> 16;
      r[11 * i + 2] |= t[1] << 6;
      r[11 * i + 3] = t[1] >> 2;
      r[11 * i + 4] = t[1] >> 10;
      r[11 * i + 5] = t[1] >> 18;
      r[11 * i + 5] |= t[2] << 4;
      r[11 * i + 6] = t[2] >> 4;
      r[11 * i + 7] = t[2] >> 12;
      r[11 * i + 8] = t[2] >> 20;
      r[11 * i + 8] |= t[3] << 2;
      r[11 * i + 9] = t[3] >> 6;
      r[11 * i + 10] = t[3] >> 14;
    }
  }
  else if (logeta == 23)
  {
    Q_SIZE t[8];
    for (i = 0; i < N / 8; ++i)
    {

      t[0] = (1 << (logeta - 1)) - a->coeffs[8 * i + 0];
      t[1] = (1 << (logeta - 1)) - a->coeffs[8 * i + 1];
      t[2] = (1 << (logeta - 1)) - a->coeffs[8 * i + 2];
      t[3] = (1 << (logeta - 1)) - a->coeffs[8 * i + 3];
      t[4] = (1 << (logeta - 1)) - a->coeffs[8 * i + 4];
      t[5] = (1 << (logeta - 1)) - a->coeffs[8 * i + 5];
      t[6] = (1 << (logeta - 1)) - a->coeffs[8 * i + 6];
      t[7] = (1 << (logeta - 1)) - a->coeffs[8 * i + 7];

      r[23 * i + 0] = t[0];
      r[23 * i + 1] = t[0] >> 8;
      r[23 * i + 2] = t[0] >> 16;
      r[23 * i + 2] |= t[1] << 7;
      r[23 * i + 3] = t[1] >> 1;
      r[23 * i + 4] = t[1] >> 9;
      r[23 * i + 5] = t[1] >> 17;
      r[23 * i + 5] |= t[2] << 6;
      r[23 * i + 6] = t[2] >> 2;
      r[23 * i + 7] = t[2] >> 10;
      r[23 * i + 8] = t[2] >> 18;
      r[23 * i + 8] |= t[3] << 5;
      r[23 * i + 9] = t[3] >> 3;
      r[23 * i + 10] = t[3] >> 11;
      r[23 * i + 11] = t[3] >> 19;
      r[23 * i + 11] |= t[4] << 4;
      r[23 * i + 12] = t[4] >> 4;
      r[23 * i + 13] = t[4] >> 12;
      r[23 * i + 14] = t[4] >> 20;
      r[23 * i + 14] |= t[5] << 3;
      r[23 * i + 15] = t[5] >> 5;
      r[23 * i + 16] = t[5] >> 13;
      r[23 * i + 17] = t[5] >> 21;
      r[23 * i + 17] |= t[6] << 2;
      r[23 * i + 18] = t[6] >> 6;
      r[23 * i + 19] = t[6] >> 14;
      r[23 * i + 20] = t[6] >> 22;
      r[23 * i + 20] |= t[7] << 1;
      r[23 * i + 21] = t[7] >> 7;
      r[23 * i + 22] = t[7] >> 15;
    }
  }
  else if (logeta == 24)
  {
    Q_SIZE t;
    for (i = 0; i < N / 1; ++i)
    {

      t = (1 << (logeta - 1)) - a->coeffs[1 * i];

      r[3 * i + 0] = t;
      r[3 * i + 1] = t >> 8;
      r[3 * i + 2] = t >> 16;
    }
  }
  else if (logeta == 25)
  {
    Q_SIZE t[8];
    for (i = 0; i < N / 8; ++i)
    {

      t[0] = (1 << (logeta - 1)) - a->coeffs[8 * i + 0];
      t[1] = (1 << (logeta - 1)) - a->coeffs[8 * i + 1];
      t[2] = (1 << (logeta - 1)) - a->coeffs[8 * i + 2];
      t[3] = (1 << (logeta - 1)) - a->coeffs[8 * i + 3];
      t[4] = (1 << (logeta - 1)) - a->coeffs[8 * i + 4];
      t[5] = (1 << (logeta - 1)) - a->coeffs[8 * i + 5];
      t[6] = (1 << (logeta - 1)) - a->coeffs[8 * i + 6];
      t[7] = (1 << (logeta - 1)) - a->coeffs[8 * i + 7];

      r[25 * i + 0] = t[0];
      r[25 * i + 1] = t[0] >> 8;
      r[25 * i + 2] = t[0] >> 16;
      r[25 * i + 3] = t[0] >> 24;
      r[25 * i + 3] |= t[1] << 1;
      r[25 * i + 4] = t[1] >> 7;
      r[25 * i + 5] = t[1] >> 15;
      r[25 * i + 6] = t[1] >> 23;
      r[25 * i + 6] |= t[2] << 2;
      r[25 * i + 7] = t[2] >> 6;
      r[25 * i + 8] = t[2] >> 14;
      r[25 * i + 9] = t[2] >> 22;
      r[25 * i + 9] |= t[3] << 3;
      r[25 * i + 10] = t[3] >> 5;
      r[25 * i + 11] = t[3] >> 13;
      r[25 * i + 12] = t[3] >> 21;
      r[25 * i + 12] |= t[4] << 4;
      r[25 * i + 13] = t[4] >> 4;
      r[25 * i + 14] = t[4] >> 12;
      r[25 * i + 15] = t[4] >> 20;
      r[25 * i + 15] |= t[5] << 5;
      r[25 * i + 16] = t[5] >> 3;
      r[25 * i + 17] = t[5] >> 11;
      r[25 * i + 18] = t[5] >> 19;
      r[25 * i + 18] |= t[6] << 6;
      r[25 * i + 19] = t[6] >> 2;
      r[25 * i + 20] = t[6] >> 10;
      r[25 * i + 21] = t[6] >> 18;
      r[25 * i + 21] |= t[7] << 7;
      r[25 * i + 22] = t[7] >> 1;
      r[25 * i + 23] = t[7] >> 9;
      r[25 * i + 24] = t[7] >> 17;
    }
  }
}

void poly_unpack_ETA(poly *r, const uint8_t *a, unsigned int logeta)
{
  unsigned int i;
  if (logeta == 2)
  {
    for (i = 0; i < N / 4; ++i)
    {
      r->coeffs[4 * i + 0] = a[i];
      r->coeffs[4 * i + 0] &= 0x3;

      r->coeffs[4 * i + 1] = a[i] >> 2;
      r->coeffs[4 * i + 1] &= 0x3;

      r->coeffs[4 * i + 2] = a[i] >> 4;
      r->coeffs[4 * i + 2] &= 0x3;

      r->coeffs[4 * i + 3] = a[i] >> 6;
      r->coeffs[4 * i + 3] &= 0x3;

      r->coeffs[4 * i + 0] = (1 << (logeta - 1)) - r->coeffs[4 * i + 0];
      r->coeffs[4 * i + 1] = (1 << (logeta - 1)) - r->coeffs[4 * i + 1];
      r->coeffs[4 * i + 2] = (1 << (logeta - 1)) - r->coeffs[4 * i + 2];
      r->coeffs[4 * i + 3] = (1 << (logeta - 1)) - r->coeffs[4 * i + 3];
    }
  }
  else if (logeta == 3)
  {
    for (i = 0; i < N / 8; ++i)
    {
      r->coeffs[8 * i + 0] = a[3 * i + 0];
      r->coeffs[8 * i + 0] &= 0x7;

      r->coeffs[8 * i + 1] = a[3 * i + 0] >> 3;
      r->coeffs[8 * i + 1] &= 0x7;

      r->coeffs[8 * i + 2] = a[3 * i + 0] >> 6;
      r->coeffs[8 * i + 2] |= (Q_SIZE)a[3 * i + 1] << 2;
      r->coeffs[8 * i + 2] &= 0x7;

      r->coeffs[8 * i + 3] = a[3 * i + 1] >> 1;
      r->coeffs[8 * i + 3] &= 0x7;

      r->coeffs[8 * i + 4] = a[3 * i + 1] >> 4;
      r->coeffs[8 * i + 4] &= 0x7;

      r->coeffs[8 * i + 5] = a[3 * i + 1] >> 7;
      r->coeffs[8 * i + 5] |= (Q_SIZE)a[3 * i + 2] << 1;
      r->coeffs[8 * i + 5] &= 0x7;

      r->coeffs[8 * i + 6] = a[3 * i + 2] >> 2;
      r->coeffs[8 * i + 6] &= 0x7;

      r->coeffs[8 * i + 7] = a[3 * i + 2] >> 5;
      r->coeffs[8 * i + 7] &= 0x7;

      r->coeffs[8 * i + 0] = (1 << (logeta - 1)) - r->coeffs[8 * i + 0];
      r->coeffs[8 * i + 1] = (1 << (logeta - 1)) - r->coeffs[8 * i + 1];
      r->coeffs[8 * i + 2] = (1 << (logeta - 1)) - r->coeffs[8 * i + 2];
      r->coeffs[8 * i + 3] = (1 << (logeta - 1)) - r->coeffs[8 * i + 3];
      r->coeffs[8 * i + 4] = (1 << (logeta - 1)) - r->coeffs[8 * i + 4];
      r->coeffs[8 * i + 5] = (1 << (logeta - 1)) - r->coeffs[8 * i + 5];
      r->coeffs[8 * i + 6] = (1 << (logeta - 1)) - r->coeffs[8 * i + 6];
      r->coeffs[8 * i + 7] = (1 << (logeta - 1)) - r->coeffs[8 * i + 7];
    }
  }
  else if (logeta == 18)
  {
    for (i = 0; i < N / 4; ++i)
    {

      r->coeffs[4 * i + 0] = a[9 * i + 0];
      r->coeffs[4 * i + 0] |= (Q_SIZE)a[9 * i + 1] << 8;
      r->coeffs[4 * i + 0] |= (Q_SIZE)a[9 * i + 2] << 16;
      r->coeffs[4 * i + 0] &= 0x3ffff;

      r->coeffs[4 * i + 1] = a[9 * i + 2] >> 2;
      r->coeffs[4 * i + 1] |= (Q_SIZE)a[9 * i + 3] << 6;
      r->coeffs[4 * i + 1] |= (Q_SIZE)a[9 * i + 4] << 14;
      r->coeffs[4 * i + 1] &= 0x3ffff;

      r->coeffs[4 * i + 2] = a[9 * i + 4] >> 4;
      r->coeffs[4 * i + 2] |= (Q_SIZE)a[9 * i + 5] << 4;
      r->coeffs[4 * i + 2] |= (Q_SIZE)a[9 * i + 6] << 12;
      r->coeffs[4 * i + 2] &= 0x3ffff;

      r->coeffs[4 * i + 3] = a[9 * i + 6] >> 6;
      r->coeffs[4 * i + 3] |= (Q_SIZE)a[9 * i + 7] << 2;
      r->coeffs[4 * i + 3] |= (Q_SIZE)a[9 * i + 8] << 10;
      r->coeffs[4 * i + 3] &= 0x3ffff;

      r->coeffs[4 * i + 0] = (1 << (logeta - 1)) - r->coeffs[4 * i + 0];
      r->coeffs[4 * i + 1] = (1 << (logeta - 1)) - r->coeffs[4 * i + 1];
      r->coeffs[4 * i + 2] = (1 << (logeta - 1)) - r->coeffs[4 * i + 2];
      r->coeffs[4 * i + 3] = (1 << (logeta - 1)) - r->coeffs[4 * i + 3];
    }
  }
  else if (logeta == 21)
  {
    for (i = 0; i < N / 8; ++i)
    {

      r->coeffs[8 * i + 0] = a[21 * i + 0];
      r->coeffs[8 * i + 0] |= (Q_SIZE)a[21 * i + 1] << 8;
      r->coeffs[8 * i + 0] |= (Q_SIZE)a[21 * i + 2] << 16;
      r->coeffs[8 * i + 0] &= 0x1fffff;

      r->coeffs[8 * i + 1] = a[21 * i + 2] >> 5;
      r->coeffs[8 * i + 1] |= (Q_SIZE)a[21 * i + 3] << 3;
      r->coeffs[8 * i + 1] |= (Q_SIZE)a[21 * i + 4] << 11;
      r->coeffs[8 * i + 1] |= (Q_SIZE)a[21 * i + 5] << 19;
      r->coeffs[8 * i + 1] &= 0x1fffff;

      r->coeffs[8 * i + 2] = a[21 * i + 5] >> 2;
      r->coeffs[8 * i + 2] |= (Q_SIZE)a[21 * i + 6] << 6;
      r->coeffs[8 * i + 2] |= (Q_SIZE)a[21 * i + 7] << 14;
      r->coeffs[8 * i + 2] &= 0x1fffff;

      r->coeffs[8 * i + 3] = a[21 * i + 7] >> 7;
      r->coeffs[8 * i + 3] |= (Q_SIZE)a[21 * i + 8] << 1;
      r->coeffs[8 * i + 3] |= (Q_SIZE)a[21 * i + 9] << 9;
      r->coeffs[8 * i + 3] |= (Q_SIZE)a[21 * i + 10] << 17;
      r->coeffs[8 * i + 3] &= 0x1fffff;

      r->coeffs[8 * i + 4] = a[21 * i + 10] >> 4;
      r->coeffs[8 * i + 4] |= (Q_SIZE)a[21 * i + 11] << 4;
      r->coeffs[8 * i + 4] |= (Q_SIZE)a[21 * i + 12] << 12;
      r->coeffs[8 * i + 4] |= (Q_SIZE)a[21 * i + 13] << 20;
      r->coeffs[8 * i + 4] &= 0x1fffff;

      r->coeffs[8 * i + 5] = a[21 * i + 13] >> 1;
      r->coeffs[8 * i + 5] |= (Q_SIZE)a[21 * i + 14] << 7;
      r->coeffs[8 * i + 5] |= (Q_SIZE)a[21 * i + 15] << 15;
      r->coeffs[8 * i + 5] &= 0x1fffff;

      r->coeffs[8 * i + 6] = a[21 * i + 15] >> 6;
      r->coeffs[8 * i + 6] |= (Q_SIZE)a[21 * i + 16] << 2;
      r->coeffs[8 * i + 6] |= (Q_SIZE)a[21 * i + 17] << 10;
      r->coeffs[8 * i + 6] |= (Q_SIZE)a[21 * i + 18] << 18;
      r->coeffs[8 * i + 6] &= 0x1fffff;

      r->coeffs[8 * i + 7] = a[21 * i + 18] >> 3;
      r->coeffs[8 * i + 7] |= (Q_SIZE)a[21 * i + 19] << 5;
      r->coeffs[8 * i + 7] |= (Q_SIZE)a[21 * i + 20] << 13;
      r->coeffs[8 * i + 7] &= 0x1fffff;

      r->coeffs[8 * i + 0] = (1 << (logeta - 1)) - r->coeffs[8 * i + 0];
      r->coeffs[8 * i + 1] = (1 << (logeta - 1)) - r->coeffs[8 * i + 1];
      r->coeffs[8 * i + 2] = (1 << (logeta - 1)) - r->coeffs[8 * i + 2];
      r->coeffs[8 * i + 3] = (1 << (logeta - 1)) - r->coeffs[8 * i + 3];
      r->coeffs[8 * i + 4] = (1 << (logeta - 1)) - r->coeffs[8 * i + 4];
      r->coeffs[8 * i + 5] = (1 << (logeta - 1)) - r->coeffs[8 * i + 5];
      r->coeffs[8 * i + 6] = (1 << (logeta - 1)) - r->coeffs[8 * i + 6];
      r->coeffs[8 * i + 7] = (1 << (logeta - 1)) - r->coeffs[8 * i + 7];
    }
  }
  else if (logeta == 22)
  {
    for (i = 0; i < N / 4; ++i)
    {

      r->coeffs[4 * i + 0] = a[11 * i + 0];
      r->coeffs[4 * i + 0] |= (Q_SIZE)a[11 * i + 1] << 8;
      r->coeffs[4 * i + 0] |= (Q_SIZE)a[11 * i + 2] << 16;
      r->coeffs[4 * i + 0] &= 0x3fffff;

      r->coeffs[4 * i + 1] = a[11 * i + 2] >> 6;
      r->coeffs[4 * i + 1] |= (Q_SIZE)a[11 * i + 3] << 2;
      r->coeffs[4 * i + 1] |= (Q_SIZE)a[11 * i + 4] << 10;
      r->coeffs[4 * i + 1] |= (Q_SIZE)a[11 * i + 5] << 18;
      r->coeffs[4 * i + 1] &= 0x3fffff;

      r->coeffs[4 * i + 2] = a[11 * i + 5] >> 4;
      r->coeffs[4 * i + 2] |= (Q_SIZE)a[11 * i + 6] << 4;
      r->coeffs[4 * i + 2] |= (Q_SIZE)a[11 * i + 7] << 12;
      r->coeffs[4 * i + 2] |= (Q_SIZE)a[11 * i + 8] << 20;
      r->coeffs[4 * i + 2] &= 0x3fffff;

      r->coeffs[4 * i + 3] = a[11 * i + 8] >> 2;
      r->coeffs[4 * i + 3] |= (Q_SIZE)a[11 * i + 9] << 6;
      r->coeffs[4 * i + 3] |= (Q_SIZE)a[11 * i + 10] << 14;
      r->coeffs[4 * i + 3] &= 0x3fffff;

      r->coeffs[4 * i + 0] = (1 << (logeta - 1)) - r->coeffs[4 * i + 0];
      r->coeffs[4 * i + 1] = (1 << (logeta - 1)) - r->coeffs[4 * i + 1];
      r->coeffs[4 * i + 2] = (1 << (logeta - 1)) - r->coeffs[4 * i + 2];
      r->coeffs[4 * i + 3] = (1 << (logeta - 1)) - r->coeffs[4 * i + 3];
    }
  }
  else if (logeta == 23)
  {
    for (i = 0; i < N / 8; ++i)
    {

      r->coeffs[8 * i + 0] = a[23 * i + 0];
      r->coeffs[8 * i + 0] |= (Q_SIZE)a[23 * i + 1] << 8;
      r->coeffs[8 * i + 0] |= (Q_SIZE)a[23 * i + 2] << 16;
      r->coeffs[8 * i + 0] &= 0x7fffff;

      r->coeffs[8 * i + 1] = a[23 * i + 2] >> 7;
      r->coeffs[8 * i + 1] |= (Q_SIZE)a[23 * i + 3] << 1;
      r->coeffs[8 * i + 1] |= (Q_SIZE)a[23 * i + 4] << 9;
      r->coeffs[8 * i + 1] |= (Q_SIZE)a[23 * i + 5] << 17;
      r->coeffs[8 * i + 1] &= 0x7fffff;

      r->coeffs[8 * i + 2] = a[23 * i + 5] >> 6;
      r->coeffs[8 * i + 2] |= (Q_SIZE)a[23 * i + 6] << 2;
      r->coeffs[8 * i + 2] |= (Q_SIZE)a[23 * i + 7] << 10;
      r->coeffs[8 * i + 2] |= (Q_SIZE)a[23 * i + 8] << 18;
      r->coeffs[8 * i + 2] &= 0x7fffff;

      r->coeffs[8 * i + 3] = a[23 * i + 8] >> 5;
      r->coeffs[8 * i + 3] |= (Q_SIZE)a[23 * i + 9] << 3;
      r->coeffs[8 * i + 3] |= (Q_SIZE)a[23 * i + 10] << 11;
      r->coeffs[8 * i + 3] |= (Q_SIZE)a[23 * i + 11] << 19;
      r->coeffs[8 * i + 3] &= 0x7fffff;

      r->coeffs[8 * i + 4] = a[23 * i + 11] >> 4;
      r->coeffs[8 * i + 4] |= (Q_SIZE)a[23 * i + 12] << 4;
      r->coeffs[8 * i + 4] |= (Q_SIZE)a[23 * i + 13] << 12;
      r->coeffs[8 * i + 4] |= (Q_SIZE)a[23 * i + 14] << 20;
      r->coeffs[8 * i + 4] &= 0x7fffff;

      r->coeffs[8 * i + 5] = a[23 * i + 14] >> 3;
      r->coeffs[8 * i + 5] |= (Q_SIZE)a[23 * i + 15] << 5;
      r->coeffs[8 * i + 5] |= (Q_SIZE)a[23 * i + 16] << 13;
      r->coeffs[8 * i + 5] |= (Q_SIZE)a[23 * i + 17] << 21;
      r->coeffs[8 * i + 5] &= 0x7fffff;

      r->coeffs[8 * i + 6] = a[23 * i + 17] >> 2;
      r->coeffs[8 * i + 6] |= (Q_SIZE)a[23 * i + 18] << 6;
      r->coeffs[8 * i + 6] |= (Q_SIZE)a[23 * i + 19] << 14;
      r->coeffs[8 * i + 6] |= (Q_SIZE)a[23 * i + 20] << 22;
      r->coeffs[8 * i + 6] &= 0x7fffff;

      r->coeffs[8 * i + 7] = a[23 * i + 20] >> 1;
      r->coeffs[8 * i + 7] |= (Q_SIZE)a[23 * i + 21] << 7;
      r->coeffs[8 * i + 7] |= (Q_SIZE)a[23 * i + 22] << 15;
      r->coeffs[8 * i + 7] &= 0x7fffff;

      r->coeffs[8 * i + 0] = (1 << (logeta - 1)) - r->coeffs[8 * i + 0];
      r->coeffs[8 * i + 1] = (1 << (logeta - 1)) - r->coeffs[8 * i + 1];
      r->coeffs[8 * i + 2] = (1 << (logeta - 1)) - r->coeffs[8 * i + 2];
      r->coeffs[8 * i + 3] = (1 << (logeta - 1)) - r->coeffs[8 * i + 3];
      r->coeffs[8 * i + 4] = (1 << (logeta - 1)) - r->coeffs[8 * i + 4];
      r->coeffs[8 * i + 5] = (1 << (logeta - 1)) - r->coeffs[8 * i + 5];
      r->coeffs[8 * i + 6] = (1 << (logeta - 1)) - r->coeffs[8 * i + 6];
      r->coeffs[8 * i + 7] = (1 << (logeta - 1)) - r->coeffs[8 * i + 7];
    }
  }
  else if (logeta == 24)
  {
    for (i = 0; i < N / 1; ++i)
    {

      r->coeffs[i] = a[3 * i + 0];
      r->coeffs[i] |= (Q_SIZE)a[3 * i + 1] << 8;
      r->coeffs[i] |= (Q_SIZE)a[3 * i + 2] << 16;
      r->coeffs[i] &= 0xffffff;

      r->coeffs[1 * i] = (1 << (logeta - 1)) - r->coeffs[1 * i];
    }
  }
  else if (logeta == 25)
  {
    for (i = 0; i < N / 8; ++i)
    {

      r->coeffs[8 * i + 0] = a[25 * i + 0];
      r->coeffs[8 * i + 0] |= (Q_SIZE)a[25 * i + 1] << 8;
      r->coeffs[8 * i + 0] |= (Q_SIZE)a[25 * i + 2] << 16;
      r->coeffs[8 * i + 0] |= (Q_SIZE)a[25 * i + 3] << 24;
      r->coeffs[8 * i + 0] &= 0x1ffffff;

      r->coeffs[8 * i + 1] = a[25 * i + 3] >> 1;
      r->coeffs[8 * i + 1] |= (Q_SIZE)a[25 * i + 4] << 7;
      r->coeffs[8 * i + 1] |= (Q_SIZE)a[25 * i + 5] << 15;
      r->coeffs[8 * i + 1] |= (Q_SIZE)a[25 * i + 6] << 23;
      r->coeffs[8 * i + 1] &= 0x1ffffff;

      r->coeffs[8 * i + 2] = a[25 * i + 6] >> 2;
      r->coeffs[8 * i + 2] |= (Q_SIZE)a[25 * i + 7] << 6;
      r->coeffs[8 * i + 2] |= (Q_SIZE)a[25 * i + 8] << 14;
      r->coeffs[8 * i + 2] |= (Q_SIZE)a[25 * i + 9] << 22;
      r->coeffs[8 * i + 2] &= 0x1ffffff;

      r->coeffs[8 * i + 3] = a[25 * i + 9] >> 3;
      r->coeffs[8 * i + 3] |= (Q_SIZE)a[25 * i + 10] << 5;
      r->coeffs[8 * i + 3] |= (Q_SIZE)a[25 * i + 11] << 13;
      r->coeffs[8 * i + 3] |= (Q_SIZE)a[25 * i + 12] << 21;
      r->coeffs[8 * i + 3] &= 0x1ffffff;

      r->coeffs[8 * i + 4] = a[25 * i + 12] >> 4;
      r->coeffs[8 * i + 4] |= (Q_SIZE)a[25 * i + 13] << 4;
      r->coeffs[8 * i + 4] |= (Q_SIZE)a[25 * i + 14] << 12;
      r->coeffs[8 * i + 4] |= (Q_SIZE)a[25 * i + 15] << 20;
      r->coeffs[8 * i + 4] &= 0x1ffffff;

      r->coeffs[8 * i + 5] = a[25 * i + 15] >> 5;
      r->coeffs[8 * i + 5] |= (Q_SIZE)a[25 * i + 16] << 3;
      r->coeffs[8 * i + 5] |= (Q_SIZE)a[25 * i + 17] << 11;
      r->coeffs[8 * i + 5] |= (Q_SIZE)a[25 * i + 18] << 19;
      r->coeffs[8 * i + 5] &= 0x1ffffff;

      r->coeffs[8 * i + 6] = a[25 * i + 18] >> 6;
      r->coeffs[8 * i + 6] |= (Q_SIZE)a[25 * i + 19] << 2;
      r->coeffs[8 * i + 6] |= (Q_SIZE)a[25 * i + 20] << 10;
      r->coeffs[8 * i + 6] |= (Q_SIZE)a[25 * i + 21] << 18;
      r->coeffs[8 * i + 6] &= 0x1ffffff;

      r->coeffs[8 * i + 7] = a[25 * i + 21] >> 7;
      r->coeffs[8 * i + 7] |= (Q_SIZE)a[25 * i + 22] << 1;
      r->coeffs[8 * i + 7] |= (Q_SIZE)a[25 * i + 23] << 9;
      r->coeffs[8 * i + 7] |= (Q_SIZE)a[25 * i + 24] << 17;
      r->coeffs[8 * i + 7] &= 0x1ffffff;

      r->coeffs[8 * i + 0] = (1 << (logeta - 1)) - r->coeffs[8 * i + 0];
      r->coeffs[8 * i + 1] = (1 << (logeta - 1)) - r->coeffs[8 * i + 1];
      r->coeffs[8 * i + 2] = (1 << (logeta - 1)) - r->coeffs[8 * i + 2];
      r->coeffs[8 * i + 3] = (1 << (logeta - 1)) - r->coeffs[8 * i + 3];
      r->coeffs[8 * i + 4] = (1 << (logeta - 1)) - r->coeffs[8 * i + 4];
      r->coeffs[8 * i + 5] = (1 << (logeta - 1)) - r->coeffs[8 * i + 5];
      r->coeffs[8 * i + 6] = (1 << (logeta - 1)) - r->coeffs[8 * i + 6];
      r->coeffs[8 * i + 7] = (1 << (logeta - 1)) - r->coeffs[8 * i + 7];
    }
  }
}