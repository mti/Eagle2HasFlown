/*
 * Implementors: EagleSign Team
 * This implementation is highly inspired from Dilithium and
 * Falcon Signatures' implementations
 */

#ifndef POLYMATRIX_H
#define POLYMATRIX_H

#include <stdint.h>
#include "params.h"
#include "polyvec.h"
#include "poly.h"

#define polymatrix_expand EAGLESIGN_NAMESPACE(polymatrix_expand)
void polymatrix_expand(polyvecl mat[K], const uint8_t rho[SEEDBYTES]);

#define polymatrix_pointwise_montgomery EAGLESIGN_NAMESPACE(polymatrix_pointwise_montgomery)
void polymatrix_pointwise_montgomery(polyveck *t, const polyvecl mat[K], const polyvecl *v);

#define polymatrix_pointwise_montgomery_l_l EAGLESIGN_NAMESPACE(polymatrix_pointwise_montgomery_l_l)
void polymatrix_pointwise_montgomery_l_l(polyvecl *t, const polyvecl mat[L], const polyvecl *v);

#define polymatrix_pointwise_product EAGLESIGN_NAMESPACE(polymatrix_pointwise_product)
void polymatrix_pointwise_product(polyveck c[L], const polyvecl a[K], const polyvecl b[L]);

#define polymatrix_pointwise_product_l_l EAGLESIGN_NAMESPACE(polymatrix_pointwise_product_l_l)
void polymatrix_pointwise_product_l_l(polyvecl c[L], const polyvecl a[L], const polyvecl b[L]);

#define poly_pointwise_matrix_product EAGLESIGN_NAMESPACE(poly_pointwise_matrix_product)
void poly_pointwise_matrix_product(polyvecl c[K], const poly a, const polyvecl b[K]);

#define poly_pointwise_matrix_product_l_l EAGLESIGN_NAMESPACE(poly_pointwise_matrix_product_l_l)
void poly_pointwise_matrix_product_l_l(polyvecl c[L], const poly a, const polyvecl b[L]);

#define polymatrix_reformat EAGLESIGN_NAMESPACE(polymatrix_reformat)
void polymatrix_reformat(polyvecl b[K], const polyveck a[L]);

#define polymatrix_reformat_l_l EAGLESIGN_NAMESPACE(polymatrix_reformat_l_l)
void polymatrix_reformat_l_l(polyvecl b[L], const polyvecl a[L]);

#define polymatrix_l_expand_f EAGLESIGN_NAMESPACE(polymatrix_l_expand_f)
void polymatrix_l_expand_f(polyvecl v[L], const uint8_t seed[CRHBYTES]);

#define polymatrix_k_l_expand_d EAGLESIGN_NAMESPACE(polymatrix_k_l_expand_d)
void polymatrix_k_l_expand_d(polyvecl v[K], const uint8_t seed[SEEDBYTES]);

#define polymatrix_ntt_k_l EAGLESIGN_NAMESPACE(polymatrix_ntt_k_l)
void polymatrix_ntt_k_l(polyvecl a[K]);

#define polymatrix_ntt_l_l EAGLESIGN_NAMESPACE(polymatrix_ntt_l_l)
void polymatrix_ntt_l_l(polyvecl a[L]);

#define polymatrix_invntt_tomont_k_l EAGLESIGN_NAMESPACE(polymatrix_invntt_tomont_k_l)
void polymatrix_invntt_tomont_k_l(polyvecl a[K]);

#define polymatrix_invntt_tomont_l_l EAGLESIGN_NAMESPACE(polymatrix_invntt_tomont_l_l)
void polymatrix_invntt_tomont_l_l(polyvecl a[L]);

#define polymatrix_l_is_invertible EAGLESIGN_NAMESPACE(polymatrix_l_is_invertible)
int polymatrix_l_is_invertible(poly *d, const polyvecl v[L]);

#define polymatrix_l_expand_f_invertible EAGLESIGN_NAMESPACE(polymatrix_l_expand_f_invertible)
int polymatrix_l_expand_f_invertible(polyvecl v[L], polyvecl f[L], const uint8_t seed[CRHBYTES]);

#define poly_is_invertible EAGLESIGN_NAMESPACE(poly_is_invertible)
int poly_is_invertible(poly *d);

#define poly_expand_g_invertible EAGLESIGN_NAMESPACE(poly_expand_g_invertible)
int poly_expand_g_invertible(poly *v, poly *g, const uint8_t seed[SEEDBYTES]);

#define polymatrix_pointwise_add EAGLESIGN_NAMESPACE(polymatrix_pointwise_add)
void polymatrix_pointwise_add(polyvecl c[K], const polyvecl a[K], const polyvecl b[K]);

#define print_byte_array EAGLESIGN_NAMESPACE(print_byte_array)
void print_byte_array(const uint8_t *pk, int size);

#define print_poly EAGLESIGN_NAMESPACE(print_poly)
void print_poly(const poly *s);

#define print_polyveck EAGLESIGN_NAMESPACE(print_polyveck)
void print_polyveck(const polyveck s, int l);

#define print_polyvecl EAGLESIGN_NAMESPACE(print_polyvecl)
void print_polyvecl(const polyvecl s, int l);

#define print_polymatrix EAGLESIGN_NAMESPACE(print_polymatrix)
void print_polymatrix(const polyvecl *s, int m, int l);

#endif
