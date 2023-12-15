/*
 * Implementors: EagleSign Team
 * This implementation is highly inspired from Dilithium and
 * Falcon Signatures' implementations
 */

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "../randombytes.h"
#include "../fips202.h"
#include "../params.h"
#include "../sign.h"
#include "../poly.h"
#include "../polyvec.h"
#include "../polymatrix.h"
#include "../packing.h"

#define MLEN 32
#define NVECTORS 100

void randombytes(uint8_t *out, size_t outlen)
{
  unsigned int i;
  uint8_t buf[8];
  static uint64_t ctr = 0;

  for (i = 0; i < 8; ++i)
    buf[i] = ctr >> 8 * i;

  ctr++;
  shake128(out, outlen, buf, 8);
}

int main(void)
{
  unsigned int i, j, k, l;
  uint8_t pk[CRYPTO_EAGLESIGN_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_EAGLESIGN_SECRETKEYBYTES];
  uint8_t sig[CRYPTO_EAGLESIGN_BYTES];
  uint8_t m[MLEN];
  uint8_t seed[CRHBYTES];
  uint8_t buf[CRYPTO_EAGLESIGN_SECRETKEYBYTES];
  size_t siglen;
  poly c, tmp, g_inv, g;
  polyvecl A[K], F[L], D[K], F_INV[L], Y1;
  polyveck Y2;
  uint16_t nonce = 0;

  for (i = 0; i < NVECTORS; ++i)
  {
    printf("count = %u\n", i);

    randombytes(m, MLEN);
    printf("m = ");
    for (j = 0; j < MLEN; ++j)
      printf("%02x", m[j]);
    printf("\n");

    crypto_sign_keypair(pk, sk);
    shake256(buf, 32, pk, CRYPTO_EAGLESIGN_PUBLICKEYBYTES);
    printf("pk = ");
    for (j = 0; j < 32; ++j)
      printf("%02x", buf[j]);
    printf("\n");
    shake256(buf, 32, sk, CRYPTO_EAGLESIGN_SECRETKEYBYTES);
    printf("sk = ");
    for (j = 0; j < 32; ++j)
      printf("%02x", buf[j]);
    printf("\n");

    crypto_sign_signature(sig, &siglen, m, MLEN, sk);
    shake256(buf, 32, sig, CRYPTO_EAGLESIGN_BYTES);
    printf("sig = ");
    for (j = 0; j < 32; ++j)
      printf("%02x", buf[j]);
    printf("\n");

    if (crypto_sign_verify(sig, siglen, m, MLEN, pk))
      fprintf(stderr, "Signature verification failed!\n");

    randombytes(seed, sizeof(seed));
    printf("seed = ");
    for (j = 0; j < sizeof(seed); ++j)
      printf("%02X", seed[j]);
    printf("\n");

    polymatrix_expand(A, seed);
    printf("A = ([");
    for (j = 0; j < K; ++j)
    {
      for (k = 0; k < L; ++k)
      {
        for (l = 0; l < N; ++l)
        {
          printf("%8d", A[j].vec[k].coeffs[l]);
          if (l < N - 1)
            printf(", ");
          else if (k < L - 1)
            printf("], [");
          else if (j < K - 1)
            printf("];\n     [");
          else
            printf("])\n");
        }
      }
    }

    poly_expand_g_invertible(&g_inv, &g, seed);
    poly_invntt_tomont(&g);

    if (poly_chknorm(&g, 2))
      fprintf(stderr, "ERROR in poly_expand_g_invertible!\n");

    printf("g = [");
    for (l = 0; l < N; ++l)
    {
      printf("%8d", g.coeffs[l]);
      if (l < N - 1)
        printf(", ");
      else
        printf("]\n");
    }

    polymatrix_k_l_expand_d(D, seed);
    polymatrix_invntt_tomont_k_l(D);

    if (polyvecl_chknorm(&D[0], 2))
      fprintf(stderr, "ERROR in polymatrix_invntt_tomont_k_l!\n");

    printf("D = ([");
    for (j = 0; j < K; ++j)
    {
      for (k = 0; k < L; ++k)
      {
        for (l = 0; l < N; ++l)
        {
          printf("%8d", D[j].vec[k].coeffs[l]);
          if (l < N - 1)
            printf(", ");
          else if (k < L - 1)
            printf("], [");
          else if (j < K - 1)
            printf("];\n     [");
          else
            printf("])\n");
        }
      }
    }

    polymatrix_l_expand_f_invertible(F_INV, F, seed);
    polymatrix_invntt_tomont_l_l(F);

    if (polyvecl_chknorm(&F[0], LOGETAF))
      fprintf(stderr, "ERROR in polymatrix_invntt_tomont_l_l!\n");

    poly_pack_ETA(buf, &F[0].vec[0], LOGETAF);
    poly_unpack_ETA(&tmp, buf, LOGETAF);

    for (j = 0; j < N; ++j)
      if (tmp.coeffs[j] != F[0].vec[0].coeffs[j])
        fprintf(stderr, "ERROR in poly_(un)pack_ETA!\n");

    printf("F = ([");
    for (j = 0; j < L; ++j)
    {
      for (k = 0; k < L; ++k)
      {
        for (l = 0; l < N; ++l)
        {
          printf("%8d", F[j].vec[k].coeffs[l]);
          if (l < N - 1)
            printf(", ");
          else if (k < L - 1)
            printf("], [");
          else if (j < K - 1)
            printf("];\n     [");
          else
            printf("])\n");
        }
      }
    }
  }
  return 0;
}
