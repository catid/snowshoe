#include <iostream>
#include <iomanip>
#include <cassert>
using namespace std;

// Math library
#include "../snowshoe/Snowshoe.hpp"

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
 *
 * Required math operations:
 * H()
 * d = S + H*T (mod q)
 * ec_mul()
 * ec_simul()
 * a = H*T (mod q)
 */


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

	snowshoe_mul_gen(sk_c, pp_c);

	snowshoe_mul_gen(sk_s, pp_s);

	snowshoe_mul(sk_c, pp_s, sp_c);

	snowshoe_mul(sk_s, pp_c, sp_s);

	for (int ii = 0; ii < 64; ++ii) {
		if (sp_c[ii] != sp_s[ii]) {
			return false;
		}
	}

	return true;
}


//// Entrypoint

int main() {
	cout << "Snowshoe Unit Tester: EC Scalar Multiplication" << endl;

	srand(0);

	for (int ii = 0; ii < 100000; ++ii) {
		assert(ec_dh_test());
	}

	cout << "All tests passed successfully." << endl;

	return 0;
}


