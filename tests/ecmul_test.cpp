#include <iostream>
#include <iomanip>
#include <cassert>
using namespace std;

// Math library
#include "../snowshoe/ecmul.cpp"

// TODO:
// ec_gen_table_4
// ec_recode_scalars_4
// ec_table_select_4


static void fe_print(const ufe &x) {
	cout << "Real(H:L) = " << hex << x.a.i[1] << " : " << x.a.i[0] << endl;
	cout << "Imag(H:L) = " << hex << x.b.i[1] << " : " << x.b.i[0] << endl;
}

static void ec_print(const ecpt &p) {
	cout << "Point = " << endl;
	cout << "X : = " << setw(16) << hex << p.x.a.i[1] << "," << setw(16) << p.x.a.i[0] << " + i * " << setw(16) << p.x.b.i[1] << "," << setw(16) << p.x.b.i[0] << endl;
	cout << "Y : = " << setw(16) << hex << p.y.a.i[1] << "," << setw(16) << p.y.a.i[0] << " + i * " << setw(16) << p.y.b.i[1] << "," << setw(16) << p.y.b.i[0] << endl;
	cout << "Z : = " << setw(16) << hex << p.z.a.i[1] << "," << setw(16) << p.z.a.i[0] << " + i * " << setw(16) << p.z.b.i[1] << "," << setw(16) << p.z.b.i[0] << endl;
	cout << "T : = " << setw(16) << hex << p.t.a.i[1] << "," << setw(16) << p.t.a.i[0] << " + i * " << setw(16) << p.t.b.i[1] << "," << setw(16) << p.t.b.i[0] << endl;
}

static bool ec_gen_table_2_test() {
	ecpt a, b;

	ec_set(EC_G, a);
	ec_set(EC_EG, b);

	ecpt table[8];

	ec_gen_table_2(a, b, table);

	ufe t2b;

	// Add all table points together, which should sum to 16a + 8b
	ecpt p; 
	ec_set(table[0], p);
	ec_add(p, table[1], p, false, true, false, t2b);
	ec_add(p, table[2], p, false, false, false, t2b);
	ec_add(p, table[3], p, false, false, false, t2b);
	ec_add(p, table[4], p, false, false, false, t2b);
	ec_add(p, table[5], p, false, false, false, t2b);
	ec_add(p, table[6], p, false, false, false, t2b);
	ec_add(p, table[7], p, false, false, true, t2b);
	ec_neg(p, p);

	ecpt a1, b1;
	ec_dbl(b, b1, true, t2b);
	ec_dbl(b1, b1, false, t2b);
	ec_dbl(b1, b1, false, t2b);
	ec_add(b1, p, b1, false, false, true, t2b);

	ec_dbl(a, a1, true, t2b);
	ec_dbl(a1, a1, false, t2b);
	ec_dbl(a1, a1, false, t2b);
	ec_dbl(a1, a1, false, t2b);
	ec_add(a1, b1, a1, false, false, true, t2b);

	ec_affine(a1, true, a1);

	ufe one;
	fe_set_smallk(1, one);

	if (!fe_iszero(a1.x) || !fe_isequal(a1.y, one)) {
		return false;
	}

	return true;
}

static bool ec_recode_scalars_2_test(ufp a, ufp b) {
	ufp a1, b1;

	a1 = a;
	b1 = b;

	u32 lsb = ec_recode_scalars_2(a1, b1);

	// Follow the recoded bits to reconstruct the original scalars
	ufp a2, b2;
	a2.w = 0;
	b2.w = 0;

	for (int ii = 127; ii >=0; ii--) {
		u32 u = (a1.w >> ii) & 1;
		u32 v = (b1.w >> ii) & 1;

		a2.w <<= 1;

		if (u) {
			a2.w++;
		} else {
			a2.w--;
		}

		b2.w <<= 1;

		if (v) {
			if (u) {
				b2.w++;
			} else {
				b2.w--;
			}
		}
	}

	if (lsb == 1) {
		a2.w++;
	}

	if (a.w != a2.w) {
		cout << "Recoding a failed" << endl;
		return false;
	}

	if (b.w != b2.w) {
		cout << "Recoding b failed" << endl;
		return false;
	}

	return true;
}

static CAT_INLINE bool ec_isequal(const ecpt &a, const ecpt &b) {
	if (!fe_isequal(a.x, b.x)) {
		return false;
	}
	if (!fe_isequal(a.y, b.y)) {
		return false;
	}
	if (!fe_isequal(a.t, b.t)) {
		return false;
	}
	if (!fe_isequal(a.z, b.z)) {
		return false;
	}

	return true;
}

static bool ec_table_select_2_test_try(ecpt *table, u32 a, u32 b, int expected) {
	ufp a1, b1;
	ecpt r;
	a1.i[0] = a;
	b1.i[0] = b;

	ecpt c;
	if ((a & 2) == 0) {
		ec_neg(table[expected], c);
	} else {
		ec_set(table[expected], c);
	}

	ec_table_select_2(table, a1, b1, 0, r);
	if (!ec_isequal(r, c)) {
		cout << a << ", " << b << endl;
		return false;
	}
	return true;
}

static bool ec_table_select_2_test() {
	ecpt p1, p2, r;
	ec_set(EC_G, p1);
	ec_set(EC_EG, p2);

	ecpt table[8];
	ec_gen_table_2(p1, p2, table);

	bool test = !ec_table_select_2_test_try(table, 0, 0, 0);
	test |= !ec_table_select_2_test_try(table, 0, 1, 1);
	test |= !ec_table_select_2_test_try(table, 0, 2, 2);
	test |= !ec_table_select_2_test_try(table, 0, 3, 3);
	test |= !ec_table_select_2_test_try(table, 1, 0, 4);
	test |= !ec_table_select_2_test_try(table, 1, 1, 5);
	test |= !ec_table_select_2_test_try(table, 1, 2, 6);
	test |= !ec_table_select_2_test_try(table, 1, 3, 7);
	test |= !ec_table_select_2_test_try(table, 2, 0, 4);
	test |= !ec_table_select_2_test_try(table, 2, 1, 5);
	test |= !ec_table_select_2_test_try(table, 2, 2, 6);
	test |= !ec_table_select_2_test_try(table, 2, 3, 7);
	test |= !ec_table_select_2_test_try(table, 3, 0, 0);
	test |= !ec_table_select_2_test_try(table, 3, 1, 1);
	test |= !ec_table_select_2_test_try(table, 3, 2, 2);
	test |= !ec_table_select_2_test_try(table, 3, 3, 3);

	return !test;
}

/*
 * server side ecmul:
 * h = H(stuff) > 1000
 * d = (long-term server private key) + h * (ephemeral private key) (mod q) > 1000
 * d = (A0 + A1v) + h * (B0 + B1v) (mod q)
 * (private point) = d * (client public point)
 *
 * client side ecsimul:
 * h = H(stuff) > 1000
 * (private point) = h * (client private) * (ephemeral public) + (client private) * (server public)
 */

//// Reference implementations for comparison

static bool ec_mul_ref(const u64 k[4], const ecpt &P0, ecpt &R) {
	ufe t2b;

	// Verify that curve has the specified order

	bool seen = false;

	ecpt p, g;
	ec_set(P0, g);

	for (int ii = 255; ii >= 0; --ii) {
		if (seen) {
			if (fe_iszero(p.x)) {
				cout << "Zero at dbl " << ii << endl;
				return false;
			}
			ec_dbl(p, p, false, t2b);
		}

		if ((k[ii / 64] >> (ii % 64)) & 1) {
			if (seen && fe_iszero(p.x)) {
				cout << "Zero at add " << ii << endl;
				return false;
			}
			if (seen) {
				ec_add(p, g, p, true, false, false, t2b);
			} else {
				ec_set(P0, p);
				seen = true;
			}
		}
	}

	ec_affine(p, false, R);

	return true;
}

static bool ec_simul_ref(const u64 k1[4], const ecpt &P0, const u64 k2[4], const ecpt &Q0, ecpt &R) {
	ecpt Pr, Qr;

	ec_mul_ref(k1, P0, Pr);
	ec_mul_ref(k2, Q0, Qr);

	ec_expand(Pr);
	ec_expand(Qr);

	ufe t2b;
	ec_add(Pr, Qr, R, false, true, false, t2b);

	ec_affine(R, false, R);

	return true;
}


//// Test Driver

bool ec_mul_gen_test() {
	u64 k[4] = {0};
	ecpt R1, R2;
	u8 a1[64], a2[64];

	for (int jj = 0; jj < 10000; ++jj) {
		for (int ii = 0; ii < 4; ++ii) {
			for (int jj = 0; jj < 30; ++jj) {
				k[ii] ^= (k[ii] << 3) | (rand() >> 2);
			}
		}
		ec_mask_scalar(k);

		ec_mul_ref(k, EC_G, R1);
		ec_mul_gen(k, R2);

		ec_save_xy(R1, a1);
		ec_save_xy(R2, a2);

		for (int ii = 0; ii < 64; ++ii) {
			if (a1[ii] != a2[ii]) {
				return false;
			}
		}
	}

	return true;
}

bool ec_mul_test() {
	u64 k[4] = {0};
	ecpt R1, R2;
	u8 a1[64], a2[64];

	for (int jj = 0; jj < 10000; ++jj) {
		for (int ii = 0; ii < 4; ++ii) {
			for (int jj = 0; jj < 30; ++jj) {
				k[ii] ^= (k[ii] << 3) | (rand() >> 2);
			}
		}
		ec_mask_scalar(k);

		ec_mul_ref(k, EC_G, R1);
		ec_mul(k, EC_G, R2);

		ec_save_xy(R1, a1);
		ec_save_xy(R2, a2);

		for (int ii = 0; ii < 64; ++ii) {
			if (a1[ii] != a2[ii]) {
				return false;
			}
		}
	}

	return true;
}

bool ec_simul_test() {
	u64 k1[4] = {0};
	u64 k2[4] = {0};
	ecpt R1, R2;
	u8 a1[64], a2[64];

	for (int jj = 0; jj < 10000; ++jj) {
		for (int ii = 0; ii < 4; ++ii) {
			for (int jj = 0; jj < 30; ++jj) {
				k1[ii] ^= (k1[ii] << 3) | (rand() >> 2);
			}
		}
		for (int ii = 0; ii < 4; ++ii) {
			for (int jj = 0; jj < 30; ++jj) {
				k2[ii] ^= (k2[ii] << 3) | (rand() >> 2);
			}
		}
		ec_mask_scalar(k1);
		ec_mask_scalar(k2);

		ec_simul_ref(k1, EC_G, k2, EC_EG, R1);
		ec_simul(k1, EC_G, k2, EC_EG, R2);

		ec_save_xy(R1, a1);
		ec_save_xy(R2, a2);

		for (int ii = 0; ii < 64; ++ii) {
			if (a1[ii] != a2[ii]) {
				return false;
			}
		}
	}

	return true;
}


//// Entrypoint

int main() {
	cout << "Snowshoe Unit Tester: EC Scalar Multiplication" << endl;

	srand(0);

	assert(ec_gen_table_2_test());

	ufp a, b;
	a.i[1] = 0x1af9f9557b981a24ULL;
	a.i[0] = 0xb25a5d1c138484e7ULL;
	b.i[1] = 0x13b714e78886c7d5ULL;
	b.i[0] = 0x585c40764421b75fULL;
	assert(ec_recode_scalars_2_test(a, b));
	a.i[1] = 0x018c0a3ded0f112eULL;
	a.i[0] = 0x7907e0549ac3793eULL;
	b.i[1] = 0x15b63bfe365757d5ULL;
	b.i[0] = 0xabf9db0384d24c26ULL;
	assert(ec_recode_scalars_2_test(a, b));

	assert(ec_table_select_2_test());

	assert(ec_mul_gen_test());

	assert(ec_mul_test());

	assert(ec_simul_test());

	cout << "All tests passed successfully." << endl;

	return 0;
}

