#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include "omp.h"

#include "../fips202.h"
#include "../sign.h"
#include "../packing.h"
#include "../poly.h"
#include "../polyvec.h"
#include "../polymatrix.h"
#include "../params.h"
#include "../api.h"

S_DOUBLE_Q_SIZE modInverse(S_DOUBLE_Q_SIZE a);

#define MSGLEN ( 20 )

#define NUM_THREADS 96

void print_progress(size_t count, size_t max, double u[], S_Q_SIZE v[]) {
    const int bar_width = 50;

    double progress = (double) count / max;
    int bar_length = progress * bar_width;

    double dotp = 0., normu = 0., normv = 0., corr;
    int reccoeffs = 0;
    for(int j=0; j<N; j++) {
        normu += u[j] * u[j];
        normv += v[j] * v[j];
        dotp  += u[j] * v[j];
    }
    corr = dotp / sqrt(normu*normv);
    if(normu > 0) {
        double factor = sqrt(ETAF*(ETAF+1)*N/3./normu);
        //double factor = sqrt(normv/normu);
        for(int j=0; j<N; j++) {
			S_Q_SIZE recuj = lround(u[j]*factor);
			if(recuj > ETAF)
				recuj = ETAF;
			else if(recuj < -ETAF)
				recuj = -ETAF;
			reccoeffs += (recuj == v[j]);
		}
    }

    printf("\rSigs: [");
    for (int i = 0; i < bar_width; ++i) {
        printf("%c", (i<bar_length)?'#':' ');
    }
    printf("] %.2f%% [%.2f; %d]", progress * 100, corr, reccoeffs);
    fflush(stdout);
}

void print_F(const uint8_t *private_key)
{
	uint8_t rho[SEEDBYTES], tr[SEEDBYTES];
	poly g;
	polyvecl D[K], F[L], Es[K];
	unpack_sk(rho, tr, &g, D, F, Es, private_key);

	// ================================================================================
	// Print F
	// ================================================================================
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			printf("F[%d][%d] = (", i, j);
			for (int k = 0; k < N; k++) {
				if (k % 64 == 0) printf("\n");
				printf("%2d", F[i].vec[j].coeffs[k]);
				if (k != N-1) printf(",");
			}
			printf(")\n");
		}
	}
}

// Public key:
uint8_t public_key[CRYPTO_EAGLESIGN_PUBLICKEYBYTES], rho[SEEDBYTES];
uint8_t private_key[CRYPTO_EAGLESIGN_SECRETKEYBYTES];
polyvecl E[K];
polyvecl A[K]; // matrix A is even more 'unpacked'

void determine_c(polyvecl *c, uint8_t *m, size_t mlen, uint8_t *r)
{
	uint8_t cp[SEEDBYTES], mu[CRHBYTES];
	keccak_state state;
	uint16_t nonce_c = 0;

	/* Compute mu = CRH(H(pk), msg) */
	shake256(mu, SEEDBYTES, public_key, CRYPTO_EAGLESIGN_PUBLICKEYBYTES);
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
	polyvecl_challenge(c, cp, &nonce_c, 2);

	polyvecl_invntt_tomont(c);
}

void private_key_recovery(polyvecl *predicted_F, long long num_sigs, void (*gensig)(uint8_t*, size_t*, uint8_t*))
{
	// ================================================================================

	// Unpacked signature:
	uint8_t tr[SEEDBYTES];

	// Statistics on the signatures:
	long long num_traces_par[NUM_THREADS][L] = {},
		 trace_F_par[NUM_THREADS][L][L][N] = {};
	long long num_traces[L] = {},
		 trace_F[L][L][N] = {};
	double trace_F_d[N];
	double rec_F[N];

	poly g, ginv;
	polyvecl D[K], F[L], Es[K];
	unpack_sk(rho, tr, &g, D, F, Es, private_key);

	poly_ntt(&g);
	for(int k=0; k<N; k++)
		ginv.coeffs[k] = (S_Q_SIZE)modInverse(g.coeffs[k]);

	// ================================================================================

#pragma omp parallel for num_threads(NUM_THREADS) schedule(static)
	for (long long i = 0; i < num_sigs; i++) {
		uint8_t msg[MSGLEN], sig[CRYPTO_EAGLESIGN_BYTES];
		size_t msglen;
		uint8_t sig_r[SEEDBYTES];
		polyvecl sig_z, sig_c, sig_u;
		int thread_id = omp_get_thread_num();

		gensig(msg, &msglen, sig);

		unpack_sig(sig_r, &sig_z, sig);

		determine_c(&sig_c, msg, msglen, sig_r);

		polyvecl_ntt(&sig_z);
		polyvecl_pointwise_poly_montgomery(&sig_u, &ginv, &sig_z);
		polyvecl_invntt_tomont(&sig_u);

		for (int lu = 0; lu < L; lu++) {
			for (int idx_c = 0; idx_c < N; idx_c++) {
				if (sig_c.vec[lu].coeffs[idx_c] != 0) {
					long long c0 = sig_c.vec[lu].coeffs[idx_c];
					num_traces_par[thread_id][lu]++;

					// Gather statistics on F:
					for (int lz = 0; lz < L; lz++) {
						for (int j = 0; j < N; j++) {
							int idx_z = idx_c + j, sign = 1;
							if (idx_z >= N) idx_z -= N, sign = -sign;
							long long zj = sign * (long long)sig_u.vec[lz].coeffs[idx_z];
							trace_F_par[thread_id][lz][lu][j] += c0 * zj;
						}
					}
				}
			}
		}
		if(thread_id == 0 && i%256 == 0) {
#pragma omp critical
			{
				for(int j=0; j<N; j++) {
					trace_F_d[j] = 0;
					for(int t=0; t<NUM_THREADS; t++) 
						trace_F_d[j] += (double)trace_F_par[t][0][0][j];
				}
				print_progress(NUM_THREADS*(i+1), num_sigs, trace_F_d, F[0].vec[0].coeffs);
			}
		}
	}

	printf("\n\n");

	for(int lu=0; lu<L; lu++) {
		num_traces[lu] = 0;
		for(int t=0; t<NUM_THREADS; t++)
			num_traces[lu] += num_traces_par[t][lu];
		for(int lz=0; lz<L; lz++) {
			for(int j=0; j<N; j++) {
				trace_F[lz][lu][j] = 0;
				for(int t=0; t<NUM_THREADS; t++) 
					trace_F[lz][lu][j] += trace_F_par[t][lz][lu][j];
			}
		}
	}

	for (int lu = 0; lu < L; lu++) {
		// Recover F:
		for (int lz = 0; lz < L; lz++) {
			double norm_rec_F = 0;
			for (int i = 0; i < N; i++) {
				rec_F[i] = ((double)trace_F[lz][lu][i]) / num_traces[lu];
				norm_rec_F += rec_F[i] * rec_F[i];
			}
			for (int i = 0; i < N; i++) {
				S_Q_SIZE pFi;

				rec_F[i] *= sqrt(ETAF * (ETAF+1) * N / (3 * norm_rec_F));
				pFi = lround(rec_F[i]);
				if(pFi > ETAF)
					pFi = ETAF;
				else if(pFi < -ETAF)
					pFi = -ETAF;
				predicted_F[lz].vec[lu].coeffs[i] = pFi;
			}
		}
	}
}

void compare_keys(const uint8_t *private_key, const polyvecl *predicted_F)
{
	uint8_t rho[SEEDBYTES], tr[SEEDBYTES];
	poly g, ginv;
	polyvecl D[K], F[L], Es[K];
	unpack_sk(rho, tr, &g, D, F, Es, private_key);

	for (int lu = 0; lu < L; lu++) {

		// Show results for F:
		for (int lz = 0; lz < L; lz++) {
			int n_correct = 0;
			double corrcoeff = 0., norm_F = 0., norm_pF = 0., dotp = 0.;

			for (int i = 0; i < N; i++) {
				int64_t aGi = F[lz].vec[lu].coeffs[i];
				int64_t pGi = predicted_F[lz].vec[lu].coeffs[i];
				n_correct += (aGi == pGi);
				norm_F += (double)aGi * aGi;
				norm_pF+= (double)pGi * pGi;
				dotp    += (double)aGi * pGi;
			}
			corrcoeff = dotp / sqrt(norm_F*norm_pF);
			printf("%d out of %d coefficients of F[%d,%d] are recovered succesfully. Correlation: %f.\n", n_correct, N, lz, lu, corrcoeff);
		}
	}

}




void generate_msg_sig(uint8_t *msg, size_t *msglen, uint8_t sig[CRYPTO_EAGLESIGN_BYTES]) {
	size_t siglen;

	*msglen = MSGLEN;
	for (int i = 0; i < MSGLEN; i++) {
		msg[i] = rand() & 0xFF;
	}

	assert(crypto_sign_signature(sig, &siglen, msg, MSGLEN, private_key) == 0);
	assert(siglen == CRYPTO_EAGLESIGN_BYTES);
}

int main(int argc, char **argv) {
	long long num_sigs = argc > 1 ? atoll(argv[1]) : 100000000LL;
	polyvecl predicted_F[L];

	srand(time(NULL));

	// Key Generation:
	crypto_sign_keypair(public_key, private_key);

	// Recover G:
	private_key_recovery(predicted_F, num_sigs, generate_msg_sig);

	// Print the actual private key.
	print_F(private_key);

	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			printf("predicted_F[%d][%d] = (", i, j);
			for (int k = 0; k < N; k++) {
				if (k % 64 == 0) printf("\n");
				printf("%2d", predicted_F[i].vec[j].coeffs[k]);
				if (k != N-1) printf(",");
			}
			printf(")\n");
		}
	}

	// Compare G and predicted_G:
	compare_keys(private_key, predicted_F);

	return 0;
}

// vim: ts=4
