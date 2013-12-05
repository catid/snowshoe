#include <iostream>
#include <iomanip>
#include <cassert>
using namespace std;

#include "Clock.hpp"
using namespace cat;

// Math library
#include "snowshoe.h"

static Clock m_clock;


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

	if (snowshoe_mul_gen(sk_c, 0, pp_c)) {
		return false;
	}

	if (snowshoe_mul_gen(sk_s, 0, pp_s)) {
		return false;
	}

	double s0 = m_clock.usec();
	u32 t0 = Clock::cycles();

	if (snowshoe_mul(sk_c, pp_s, sp_c)) {
		return false;
	}

	u32 t1 = Clock::cycles();
	double s1 = m_clock.usec();

	cout << "EC-DH client: " << (t1 - t0) << " cycles " << (s1 - s0) << " usec" << endl;

	s0 = m_clock.usec();
	t0 = Clock::cycles();

	if (snowshoe_mul(sk_s, pp_c, sp_s)) {
		return false;
	}

	t1 = Clock::cycles();
	s1 = m_clock.usec();

	cout << "EC-DH server: " << (t1 - t0) << " cycles " << (s1 - s0) << " usec" << endl;

	for (int ii = 0; ii < 64; ++ii) {
		if (sp_c[ii] != sp_s[ii]) {
			return false;
		}
	}

	return true;
}

/*
 * EC-DH-FS:
 *
 * An improvement on EC-DH is the concept of "forward secrecy", where if the server's long-term
 * secret key is revealed it is not any easier to decrypt past logged communication.
 *
 * This example sketches how to implement a protocol that can achieve this goal.  It requires
 * that the client generate a new public key periodically to delete the past sessions, and
 * the server should also regenerate an ephemeral key periodically.
 *
 * After the key agreement completes, the server should provide a proof that it knows the
 * secret key in its response to the client.  This message can also carry the first encrypted
 * packet.
 *
 * server side ec_mul:
 * h = H(stuff) > 1000
 * d = (long-term server private key) * h + (ephemeral private key) (mod q) > 1000
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

	// Offline: Server long-term public key generation

	generate_k(sk_s);
	snowshoe_secret_gen(sk_s);
	if (snowshoe_mul_gen(sk_s, 0, pp_s)) {
		return false;
	}

	// Online: Server ephemeral public key (changes periodically)

	generate_k(sk_e);
	snowshoe_secret_gen(sk_e);
	if (snowshoe_mul_gen(sk_e, 0, pp_e)) {
		return false;
	}

	// Online: Client ephemeral public key (changes periodically)

	generate_k(sk_c);
	snowshoe_secret_gen(sk_c);
	if (snowshoe_mul_gen(sk_c, 0, pp_c)) {
		return false;
	}

	// h = H(pp_s, pp_e, pp_c, client_nonce, server_nonce)
	generate_k(h);

	// Online: Server handles client request

	double s0 = m_clock.usec();
	u32 t0 = Clock::cycles();

	// d = sk_e + h * sk_s (mod q)
	snowshoe_mul_mod_q(h, sk_s, sk_e, d);
	if (snowshoe_mul(d, pp_c, sp_s)) {
		return false;
	}

	u32 t1 = Clock::cycles();
	double s1 = m_clock.usec();

	cout << "EC-DH-FS server: " << (t1 - t0) << " cycles " << (s1 - s0) << " usec" << endl;

	// Online: Client handles server response

	s0 = m_clock.usec();
	t0 = Clock::cycles();

	// a = h * sk_c (mod q)
	snowshoe_mul_mod_q(h, sk_c, 0, a);
	if (snowshoe_simul(sk_c, pp_e, a, pp_s, sp_c)) {
		return false;
	}

	t1 = Clock::cycles();
	s1 = m_clock.usec();

	cout << "EC-DH-FS client: " << (t1 - t0) << " cycles " << (s1 - s0) << " usec" << endl;

	for (int ii = 0; ii < 64; ++ii) {
		if (sp_c[ii] != sp_s[ii]) {
			return false;
		}
	}

	return true;
}

/*
 * EdDSA from http://ed25519.cr.yp.to/ed25519-20110926.pdf
 * Summarized pretty well here: http://blog.mozilla.org/warner/2011/11/29/ed25519-keys/
 *
 * Note that the "key massaging" would be different for my group order.
 * Instead of Ed25519 key masking, use snowshoe_secret_gen.
 *
 * Key generation:
 * 	Generate a random number k < 2^256
 * 	hi,lo = H(k)
 * 	a = snowshoe_secret_gen(lo)
 * 	A = a*G
 *
 * Sign message M:
 * 	r = H(hi,M) (mod q)
 * 	t = H(R,A,M) (mod q)
 * 	R = r*4*G
 * 	s = r + t*a (mod q)
 * 	Produce: R, s
 *
 * Verify:
 * 	u = H(R,A,M) (mod q)
 * 	nA = -A
 *	R =?= s*G + u*nA
 */

bool ec_dsa_test() {
	char a[32], h_hi_m[64], h_r_a_m[64], r[32], t[32], s[32], u[32];
	char pp_A[64], pp_R[64], pp_Rtest[64];

	// Fake hashes to avoid implementing Skein-512 and Skein-256 just for testing
	generate_k(a);
	generate_k(h_hi_m);
	generate_k(h_hi_m+32);
	generate_k(h_r_a_m);
	generate_k(h_r_a_m+32);

	// Offline precomputation:

	snowshoe_secret_gen(a);
	if (snowshoe_mul_gen(a, MULGEN_VARTIME, pp_A)) {
		return false;
	}

	// Sign:

	double s0 = m_clock.usec();
	u32 t0 = Clock::cycles();

	snowshoe_mod_q(h_hi_m, r);
	snowshoe_mod_q(h_r_a_m, t);
	if (snowshoe_mul_gen(r, MULGEN_COFACTOR, pp_R)) {
		return false;
	}
	snowshoe_mul_mod_q(a, t, r, s); // s = a * t + r (mod q)

	u32 t1 = Clock::cycles();
	double s1 = m_clock.usec();

	cout << "EdDSA sign: " << (t1 - t0) << " cycles " << (s1 - s0) << " usec" << endl;

	// Verify:

	s0 = m_clock.usec();
	t0 = Clock::cycles();

	snowshoe_mod_q(h_r_a_m, u);
	snowshoe_neg(pp_A, pp_A);
	if (snowshoe_simul_gen(s, u, pp_A, pp_Rtest)) {
		return false;
	}

	for (int ii = 0; ii < 64; ++ii) {
		if (pp_Rtest[ii] != pp_R[ii]) {
			return false;
		}
	}

	t1 = Clock::cycles();
	s1 = m_clock.usec();

	cout << "EdDSA verify: " << (t1 - t0) << " cycles " << (s1 - s0) << " usec" << endl;

	return true;
}


//// Entrypoint

static void tscTime() {
	const u32 c0 = Clock::cycles();
	const double t0 = m_clock.usec();
	const double t_end = t0 + 1000000.0;

	double t;
	u32 c;
	do {
		c = Clock::cycles();
		t = m_clock.usec();
	} while (t < t_end);

	cout << "RDTSC instruction runs at " << (c - c0)/(t - t0)/1000.0 << " GHz" << endl;
}

int main() {
	cout << "Snowshoe Unit Tester" << endl;

	// Note that assert() should not be used for crypto code since it is often compiled
	// out in release mode.  It is only used here for testing.

	m_clock.OnInitialize();

	tscTime();

	srand(0);

	// Example of verifying API level on startup:

	if (snowshoe_init()) {
		throw "Wrong snowshoe static library is linked";
	}

	for (int ii = 0; ii < 10000; ++ii) {
		assert(ec_dh_test());
	}

	for (int ii = 0; ii < 10000; ++ii) {
		assert(ec_dh_fs_test());
	}

	for (int ii = 0; ii < 10000; ++ii) {
		assert(ec_dsa_test());
	}

	cout << "All tests passed successfully." << endl;

	m_clock.OnFinalize();

	return 0;
}

