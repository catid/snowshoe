#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <cstdlib>
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
	This Quickselect routine is based on the algorithm described in
	"Numerical recipes in C", Second Edition,
	Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
	This code by Nicolas Devillard - 1998. Public domain.
*/
#define ELEM_SWAP(a,b) { register u32 t=(a);(a)=(b);(b)=t; }
static u32 quick_select(u32 arr[], int n)
{
	int low, high ;
	int median;
	int middle, ll, hh;
	low = 0 ; high = n-1 ; median = (low + high) / 2;
	for (;;) {
		if (high <= low) /* One element only */
			return arr[median] ;
		if (high == low + 1) { /* Two elements only */
			if (arr[low] > arr[high])
				ELEM_SWAP(arr[low], arr[high]) ;
			return arr[median] ;
		}
		/* Find median of low, middle and high items; swap into position low */
		middle = (low + high) / 2;
		if (arr[middle] > arr[high]) ELEM_SWAP(arr[middle], arr[high]) ;
		if (arr[low] > arr[high]) ELEM_SWAP(arr[low], arr[high]) ;
		if (arr[middle] > arr[low]) ELEM_SWAP(arr[middle], arr[low]) ;
		/* Swap low item (now in position middle) into position (low+1) */
		ELEM_SWAP(arr[middle], arr[low+1]) ;
		/* Nibble from each end towards middle, swapping items when stuck */
		ll = low + 1;
		hh = high;
		for (;;) {
			do ll++; while (arr[low] > arr[ll]) ;
			do hh--; while (arr[hh] > arr[low]) ;
			if (hh < ll)
				break;
			ELEM_SWAP(arr[ll], arr[hh]) ;
		}
		/* Swap middle item (in position low) back into correct position */
		ELEM_SWAP(arr[low], arr[hh]) ;
		/* Re-set active partition */
		if (hh <= median)
			low = ll;
		if (hh >= median)
			high = hh - 1;
	}
}
#undef ELEM_SWAP

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

	vector<u32> tc, ts;
	double wc = 0, ws = 0;

	for (int iteration = 0; iteration < 10000; ++iteration) {
		generate_k(sk_c);
		snowshoe_secret_gen(sk_c);

		generate_k(sk_s);
		snowshoe_secret_gen(sk_s);

		if (snowshoe_mul_gen(sk_c, pp_c)) {
			return false;
		}

		if (snowshoe_mul_gen(sk_s, pp_s)) {
			return false;
		}

		double s0 = m_clock.usec();
		u32 t0 = Clock::cycles();

		if (snowshoe_mul(sk_c, pp_s, sp_c)) {
			return false;
		}

		u32 t1 = Clock::cycles();
		double s1 = m_clock.usec();

		tc.push_back(t1 - t0);
		wc += s1 - s0;

		s0 = m_clock.usec();
		t0 = Clock::cycles();

		if (snowshoe_mul(sk_s, pp_c, sp_s)) {
			return false;
		}

		t1 = Clock::cycles();
		s1 = m_clock.usec();

		ts.push_back(t1 - t0);
		ws += s1 - s0;

		for (int ii = 0; ii < 64; ++ii) {
			if (sp_c[ii] != sp_s[ii]) {
				return false;
			}
		}
	}

	u32 mc = quick_select(&tc[0], (int)tc.size());
	wc /= tc.size();
	u32 ms = quick_select(&ts[0], (int)ts.size());
	ws /= ts.size();

	cout << "+ EC-DH client: `" << dec << mc << "` median cycles, `" << wc << "` avg usec" << endl;
	cout << "+ EC-DH server: `" << dec << ms << "` median cycles, `" << ws << "` avg usec" << endl;

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

	vector<u32> tc, ts, tc1;
	double wc = 0, ws = 0, wc1 = 0;

	for (int iteration = 0; iteration < 10000; ++iteration) {
		// Offline: Server long-term public key generation

		generate_k(sk_s);
		snowshoe_secret_gen(sk_s);
		if (snowshoe_mul_gen(sk_s, pp_s)) {
			return false;
		}

		// Online: Server ephemeral public key (changes periodically)

		generate_k(sk_e);
		snowshoe_secret_gen(sk_e);
		if (snowshoe_mul_gen(sk_e, pp_e)) {
			return false;
		}

		// Online: Client ephemeral public key (changes periodically)

		double s0 = m_clock.usec();
		u32 t0 = Clock::cycles();

		generate_k(sk_c);
		snowshoe_secret_gen(sk_c);
		if (snowshoe_mul_gen(sk_c, pp_c)) {
			return false;
		}

		u32 t1 = Clock::cycles();
		double s1 = m_clock.usec();

		tc.push_back(t1 - t0);
		wc += s1 - s0;

		// h = H(pp_s, pp_e, pp_c, client_nonce, server_nonce)
		generate_k(h);

		// Online: Server handles client request

		s0 = m_clock.usec();
		t0 = Clock::cycles();

		// d = sk_e + h * sk_s (mod q)
		snowshoe_mul_mod_q(h, sk_s, sk_e, d);
		if (snowshoe_mul(d, pp_c, sp_s)) {
			return false;
		}

		t1 = Clock::cycles();
		s1 = m_clock.usec();

		ts.push_back(t1 - t0);
		ws += s1 - s0;

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

		tc1.push_back(t1 - t0);
		wc1 += s1 - s0;

		for (int ii = 0; ii < 64; ++ii) {
			if (sp_c[ii] != sp_s[ii]) {
				return false;
			}
		}
	}

	u32 mc = quick_select(&tc[0], (int)tc.size());
	wc /= tc.size();
	u32 mc1 = quick_select(&tc1[0], (int)tc1.size());
	wc1 /= tc1.size();
	u32 ms = quick_select(&ts[0], (int)ts.size());
	ws /= ts.size();

	cout << "+ EC-DH-FS client gen: `" << dec << mc << "` median cycles, `" << wc << "` avg usec" << endl;
	cout << "+ EC-DH-FS server proc: `" << dec << ms << "` median cycles, `" << ws << "` avg usec" << endl;
	cout << "+ EC-DH-FS client proc: `" << dec << mc1 << "` median cycles, `" << wc1 << "` avg usec" << endl;

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
 * 	R = r*G
 * 	t = H(R,A,M) (mod q)
 * 	s = r + t*a (mod q)
 * 	Produce: R, s
 *
 * Verify:
 * 	u = H(R,A,M) (mod q)
 * 	nA = -A
 *	4*R =?= s*4*G + u*4*nA
 */

bool ec_dsa_test() {
	char a[32], h_hi_m[64], h_r_a_m[64], r[32], t[32], s[32], u[32];
	char pp_A[64], pp_R[64], pp_Rtest[64];

	vector<u32> tc, ts;
	double wc = 0, ws = 0;

	for (int iteration = 0; iteration < 10000; ++iteration) {
		// Fake hashes to avoid implementing Skein-512 and Skein-256 just for testing
		generate_k(a);
		generate_k(h_hi_m);
		generate_k(h_hi_m+32);
		generate_k(h_r_a_m);
		generate_k(h_r_a_m+32);

		// Offline precomputation:

		snowshoe_secret_gen(a);
		if (snowshoe_mul_gen(a, pp_A)) {
			return false;
		}

		// Sign:

		double s0 = m_clock.usec();
		u32 t0 = Clock::cycles();

		snowshoe_mod_q(h_hi_m, r);
		snowshoe_mod_q(h_r_a_m, t);
		if (snowshoe_mul_gen(r, pp_R)) {
			return false;
		}
		snowshoe_mul_mod_q(a, t, r, s); // s = a * t + r (mod q)

		u32 t1 = Clock::cycles();
		double s1 = m_clock.usec();

		ts.push_back(t1 - t0);
		ws += s1 - s0;

		// Verify:

		s0 = m_clock.usec();
		t0 = Clock::cycles();

		snowshoe_mod_q(h_r_a_m, u);
		snowshoe_neg(pp_A, pp_A);
		if (snowshoe_simul_gen(s, u, pp_A, pp_Rtest)) {
			return false;
		}

		if (snowshoe_equals4(pp_Rtest, pp_R)) {
			return false;
		}

		t1 = Clock::cycles();
		s1 = m_clock.usec();

		tc.push_back(t1 - t0);
		wc += s1 - s0;
	}

	u32 mc = quick_select(&tc[0], (int)tc.size());
	wc /= tc.size();
	u32 ms = quick_select(&ts[0], (int)ts.size());
	ws /= ts.size();

	cout << "+ EdDSA sign: `" << dec << ms << "` median cycles, `" << ws << "` avg usec" << endl;
	cout << "+ EdDSA verify: `" << dec << mc << "` median cycles, `" << wc << "` avg usec" << endl;

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

#ifdef CAT_HAS_VECTOR_EXTENSIONS
	cout << "Using vector extensions for table lookups! <3" << endl;
#endif

	// Note that assert() should not be used for crypto code since it is often compiled
	// out in release mode.  It is only used here for testing.

	m_clock.OnInitialize();

	tscTime();

	srand(0);

	// Example of verifying API level on startup:

	if (snowshoe_init()) {
		throw "Wrong snowshoe static library is linked";
	}

	assert(ec_dh_test());
	assert(ec_dh_fs_test());
	assert(ec_dsa_test());

	cout << "All tests passed successfully." << endl;

	m_clock.OnFinalize();

	return 0;
}

