#include "api.h"

#include "randombytes.h"
#include "crypto_hash_sha512.h"

#include "snowshoe.h"

// This is supposed to be an implementation of EdDSA using Snowshoe

int crypto_sign_snowshoe_ref_keypair(
	unsigned char *pk,
	unsigned char *sk
) {
	unsigned char extsk[64];
	int i;

	// sk = hi,lo = random
	randombytes(sk, 64);

	// A = lo * G
	snowshoe_secret_gen((char*)sk);
	if (!snowshoe_mul_gen((char*)sk, false, (char*)pk)) {
		return -1;
	}

	return 0;
}

int crypto_sign_snowshoe_ref(
	unsigned char *sm,unsigned long long *smlen,
	const unsigned char *m,unsigned long long mlen,
	const unsigned char *sk
) {
	char ger[64], r[32], s[32];
	unsigned char extsk[64];
	unsigned long long i;

	// Copy message into slot 3
	*smlen = mlen + 96;
	for (i = 0; i < mlen; ++i) {
		sm[96 + i] = m[i];
	}

	// Load hi part of secret key hash into slot 2
	for (i = 0; i < 32; ++i) {
		sm[64 + i] = sk[32 + i];
	}

	// r = H(hi, M) (mod q)
	crypto_hash_sha512(extsk, sm + 64, mlen + 32);
	snowshoe_mod_q((char*)extsk, r);

	// ger = r * 4 * G
	if (!snowshoe_mul_gen(r, true, ger)) {
		return -1;
	}

	// Load point R in slot 1
	for (i = 0; i < 64; ++i) {
		sm[i] = ger[i];
	}

	// Copy X-coordinate of point A into slot 2
	for (i = 0; i < 32; ++i) {
		sm[64 + i] = sk[64 + i];
	}

	// s = H(R, Ax, M) (mod q)
	crypto_hash_sha512(extsk, sm, mlen + 96);
	snowshoe_mod_q((char*)extsk, s);

	// s = r + s*lo (mod q)
	snowshoe_mul_mod_q((char*)sk, s, r, s);

	// Copy s into slot 2
	for (i = 0; i < 32; ++i) {
		sm[64 + i] = s[i];
	}

	// Zero stack
	for (i = 0; i < 32; ++i) {
		r[i] = 0;
	}

	return 0;
}

int crypto_sign_snowshoe_ref_open(
	unsigned char *m,unsigned long long *mlen,
	const unsigned char *sm,unsigned long long smlen,
	const unsigned char *pk
) {
	int i;
	char u[32], pkn[64], rtest[64];
	unsigned char extsk[64];

	// m = R(64) s(32) message(mlen - 96)

	*mlen = ~(unsigned long long)0;
	if (smlen < 96) {
		return -1;
	}

	// Load R into slot 1
	for (i = 0; i < 64; ++i) {
		m[i] = sm[i];
	}

	// Load Ax into slot 2
	for (i = 0; i < 32; ++i) {
		m[64 + i] = pk[i];
	}

	// Load M into slot 3
	for (i = 0; i < smlen - 96; ++i) {
		m[96 + i] = sm[96 + i];
	}

	// u = H(R, Ax, M) (mod q)
	crypto_hash_sha512(extsk, m, smlen);
	snowshoe_mod_q((char*)extsk, u);

	// pkn = -pk
	snowshoe_neg((char*)pk, pkn);

	// rtest = s * 4 * G + u * 4 * A
	if (!snowshoe_simul_gen((char*)(m + 64), u, pkn, rtest)) {
		return -1;
	}

	// Verify that rtest == r (does not need to happen in constant time)
	for (i = 0; i < 64; ++i) {
		// If verification failed,
		if (rtest[i] != sm[i]) {
			// Zero message
			for (i = 0; i < smlen - 96; ++i) {
				m[i] = 0;
			}

			return -2;
		}
	}

	// Copy message over
	for (i = 0; i < smlen - 96; ++i) {
		m[i] = sm[i + 96];
	}

	*mlen = smlen - 96;

	return 0;
}

