#include "api.h"
#include "randombytes.h"

#include "snowshoe.h"

int crypto_dh_snowshoe_ref_keypair(
	unsigned char *pk,
	unsigned char *sk
) {
	randombytes(sk, 32);
	snowshoe_secret_gen((char*)sk);
	if (!snowshoe_mul_gen((char*)sk, false, (char*)pk)) {
		return -1;
	}
	return 0;
}

int crypto_dh_snowshoe_ref(
	unsigned char *out,
	const unsigned char *pk,
	const unsigned char *sk
) {
	if (!snowshoe_mul((const char*)sk, (const char*)pk, (char*)out)) {
		return -1;
	}
	return 0;
}

