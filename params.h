/*
 * Implementors: EagleSign Team
 * This implementation is highly inspired from Dilithium and
 * Falcon Signatures' implementations
 */

#ifndef PARAMS_H
#define PARAMS_H

#include "config.h"

/*
 * Computation of Constants for NTT.
 *
 *   phi = X^n + 1
 *   Q = 1 mod 2n
 *   Q0I = -1/Q mod 2^32
 *   R = 2^32 mod Q
 *   R2 = 2^64 mod Q
 */

/* ---------First Level Parameters--------- */
#if EAGLESIGN_MODE == 2
// Fixed Parameters
#define K 1
#define L 1
#define N 1024
#define Q 2021377
#define GAMMA1 0x4000
#define GAMMA2 0x20000
#define ETAF 1
#define TD 7
#define TG 7
#define TAU 16

// Dependencies
#define Q_EXACT_SIZE_HIGHEST 0x1FFFFF
#define Q0I 2849953791
#define R 1562548
#define R2 1679445
#define LOG_N 10
#define NBYTES 128
#define LOGQ 21
#define LOGETAF 2
#define LOGGAMMA1 14
#define LOGGAMMA2 17
#define LOGTG 3
#define BETA 16

#elif EAGLESIGN_MODE == 5
// Fixed Parameters
#define K 1
#define L 1
#define N 2048
#define Q 33292289
#define GAMMA1 0x10000
#define GAMMA2 0x100000
#define ETAF 1
#define TD 13
#define TG 14
#define TAU 30

// Dependencies
#define Q_EXACT_SIZE_HIGHEST 0x1FFFFFF
#define Q0I 33292287
#define R 262015
#define R2 3160307
#define LOG_N 11
#define NBYTES 256
#define LOGQ 25
#define LOGETAF 2
#define LOGGAMMA1 16
#define LOGGAMMA2 20
#define LOGTG 4
#define BETA 30

#elif EAGLESIGN_MODE == 52
// Fixed Parameters
#define K 3
#define L 2
#define N 1024
#define Q 7340033
#define GAMMA1 0x10000
#define GAMMA2 0x100000
#define ETAF 1
#define TD 3
#define TG 16
#define TAU 16

// Dependencies
#define Q_EXACT_SIZE_HIGHEST 0x7FFFFF
#define Q0I 7340031
#define R 1047991
#define R2 3338324
#define LOG_N 10
#define NBYTES 128
#define LOGQ 23
#define LOGETAF 2
#define LOGGAMMA1 16
#define LOGGAMMA2 20
#define LOGTG 4
#define BETA 32

#endif
// ------------------------------------------

/* ---------Common Parameters-------------- */
#define Q_SIZE uint32_t
#define S_Q_SIZE int32_t
#define Q_BIT_SIZE 32
#define DOUBLE_Q_BIT_SIZE 64
#define DOUBLE_Q_SIZE uint64_t
#define S_DOUBLE_Q_SIZE int64_t
#define TWO_POWER_SIZE_Q_MINUS_ONE 0xFFFFFFFF
#define SEEDBYTES 32
#define CRHBYTES 48
#define Z_SIZE (LOGGAMMA1 + LOGTG + 1)

#define T_(s) T_##s

#define CRYPTO_EAGLESIGN_PUBLICKEYBYTES (SEEDBYTES + NBYTES * L * K * LOGQ)
#define CRYPTO_EAGLESIGN_SECRETKEYBYTES (2 * SEEDBYTES + NBYTES * (2 + K * L * 2 + L * L * LOGETAF + K * L * LOGQ))
#define CRYPTO_EAGLESIGN_BYTES (SEEDBYTES + NBYTES * L * Z_SIZE)
// ------------------------------------------
#endif
