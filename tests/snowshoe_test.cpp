#include <iostream>
#include <iomanip>
#include <cassert>
using namespace std;

#include "Clock.hpp"
using namespace cat;

// Math library
#include "../snowshoe/Snowshoe.h"

//// Test Driver

static void generate_k(char kb[32]) {
	unsigned int k[16] = {0};

	for (int ii = 0; ii < 16; ++ii) {
		for (int jj = 0; jj < 30; ++jj) {
			k[ii] ^= (k[ii] << 3) | (rand() >> 2);
		}
	}

	memcpy(kb, k, 32);
}

/*
 * EC-DH:
 *
 * client:
 * Generate sk_c, pp_c
 *
 * server:
 * Generate sk_s, pp_s
 *
 * client -> server : pp_c
 * server -> client : pp_s
 *
 * client:
 * Calculate sp_c = sk_c * pp_s
 *
 * server:
 * Calculate sp_s = sk_s * pp_c
 *
 * Validate sp_c == sp_s
 */

bool ec_dh_test() {
	char sk_c[32], sk_s[32];
	char pp_c[64], pp_s[64];
	char sp_c[64], sp_s[64];

	generate_k(sk_c);
	snowshoe_secret_gen(sk_c);

	generate_k(sk_s);
	snowshoe_secret_gen(sk_s);

	assert(snowshoe_mul_gen(sk_c, pp_c));

	assert(snowshoe_mul_gen(sk_s, pp_s));

	assert(snowshoe_mul(sk_c, pp_s, sp_c));

	assert(snowshoe_mul(sk_s, pp_c, sp_s));

	for (int ii = 0; ii < 64; ++ii) {
		if (sp_c[ii] != sp_s[ii]) {
			return false;
		}
	}

	return true;
}

/*
 * server side ec_mul:
 * h = H(stuff) > 1000
 * d = (long-term server private key) + h * (ephemeral private key) (mod q) > 1000
 * (private point) = d * 4 * (client public point)
 *
 * client side ec_simul:
 * h = H(stuff) > 1000
 * a = h * (client private) (mod q)
 * (private point) = a * 4 * (ephemeral public) + (client private) * 4 * (server public)
 */

bool ec_dh_fs_test() {
	char h[32], d[32], a[32];
	char sk_c[32], sk_s[32], sk_e[32];
	char pp_c[64], pp_s[64], pp_e[64];
	char sp_c[64], sp_s[64];

	// Offline precomputation:

	generate_k(sk_s);
	snowshoe_secret_gen(sk_s);
	assert(snowshoe_mul_gen(sk_s, pp_s));

	generate_k(sk_e);
	snowshoe_secret_gen(sk_e);
	assert(snowshoe_mul_gen(sk_e, pp_e));

	// Online: Client setup

	generate_k(sk_c);
	snowshoe_secret_gen(sk_c);

	u32 t0 = Clock::cycles();

	assert(snowshoe_mul_gen(sk_c, pp_c));

	u32 t1 = Clock::cycles();

	cout << (t1 - t0) << endl;

	generate_k(h);

	// Online: Server handles client request

	// d = h * sk_e + sk_s (mod q)
	assert(snowshoe_mul_mod_q(h, sk_e, sk_s, d));
	assert(snowshoe_mul(d, pp_c, sp_s));

	// Online: Client handles server response

	// a = h * sk_c (mod q)
	assert(snowshoe_mul_mod_q(h, sk_c, 0, a));
	assert(snowshoe_simul(a, pp_e, sk_c, pp_s, sp_c));

	for (int ii = 0; ii < 64; ++ii) {
		if (sp_c[ii] != sp_s[ii]) {
			return false;
		}
	}

	return true;
}


//// Entrypoint

int main() {
	cout << "Snowshoe Unit Tester" << endl;

	srand(0);

	for (int ii = 0; ii < 10000; ++ii) {
		assert(ec_dh_fs_test());
	}

	for (int ii = 0; ii < 10000; ++ii) {
		assert(ec_dh_test());
	}

	cout << "All tests passed successfully." << endl;

	return 0;
}

