/*
 * Implementors: EagleSign Team
 * This implementation is highly inspired from Dilithium and
 * Falcon Signatures' implementations
 */

#ifndef NTT_H
#define NTT_H

#include <stdint.h>
#include "params.h"

#define ntt EAGLESIGN_NAMESPACE(ntt)
void ntt(S_Q_SIZE a[N], unsigned logn);

#define invntt_tomont EAGLESIGN_NAMESPACE(invntt_tomont)
void invntt_tomont(S_Q_SIZE a[N], unsigned logn);

#define GMb EAGLESIGN_NAMESPACE(GMb)
extern const S_Q_SIZE GMb[N];

#define iGMb EAGLESIGN_NAMESPACE(iGMb)
extern const S_Q_SIZE iGMb[N];

#define basemul EAGLESIGN_NAMESPACE(basemul)
void basemul(S_Q_SIZE r[2],
             const S_Q_SIZE a[2],
             const S_Q_SIZE b[2],
             S_Q_SIZE zeta);

#endif
