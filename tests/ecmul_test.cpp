#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
using namespace std;

#include "Clock.hpp"
using namespace cat;

static Clock m_clock;

// Math library
#include "../src/ecmul.inc"

static void random_k(u64 k[4]) {
	for (int ii = 0; ii < 4; ++ii) {
		for (int jj = 0; jj < 30; ++jj) {
			k[ii] ^= (k[ii] << 3) | (rand() >> 2);
		}
	}
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

static const ecpt_affine EC_G_AFFINE = {
	EC_GX, EC_GY
};

static const ecpt_affine EC_EG_AFFINE = {
	EC_EGX, EC_EGY
};

static CAT_INLINE bool ec_isequal_xy(const ecpt_affine &a, const ecpt_affine &b) {
	if (!fe_isequal(a.x, b.x)) {
		return false;
	}
	if (!fe_isequal(a.y, b.y)) {
		return false;
	}

	return true;
}

static void fe_print(const ufe &x) {
	cout << "Real(H:L) = " << hex << x.a.i[1] << " : " << x.a.i[0] << endl;
	cout << "Imag(H:L) = " << hex << x.b.i[1] << " : " << x.b.i[0] << endl;
}

static void ec_print(const ecpt &p) {
	cout << "Point = " << endl;
	cout << "X : = " << setw(16) << hex << p.x.a.i[1] << "," << setw(16) << p.x.a.i[0] << " + i * " << setw(16) << p.x.b.i[1] << "," << setw(16) << p.x.b.i[0] << endl;
	cout << "Y : = " << setw(16) << hex << p.y.a.i[1] << "," << setw(16) << p.y.a.i[0] << " + i * " << setw(16) << p.y.b.i[1] << "," << setw(16) << p.y.b.i[0] << endl;
	cout << "T : = " << setw(16) << hex << p.t.a.i[1] << "," << setw(16) << p.t.a.i[0] << " + i * " << setw(16) << p.t.b.i[1] << "," << setw(16) << p.t.b.i[0] << endl;
	cout << "Z : = " << setw(16) << hex << p.z.a.i[1] << "," << setw(16) << p.z.a.i[0] << " + i * " << setw(16) << p.z.b.i[1] << "," << setw(16) << p.z.b.i[0] << endl;
}

static void ec_print_xy(const ecpt_affine &p) {
	cout << "Point = " << endl;
	cout << "X : = " << setw(16) << hex << p.x.a.i[1] << "," << setw(16) << p.x.a.i[0] << " + i * " << setw(16) << p.x.b.i[1] << "," << setw(16) << p.x.b.i[0] << endl;
	cout << "Y : = " << setw(16) << hex << p.y.a.i[1] << "," << setw(16) << p.y.a.i[0] << " + i * " << setw(16) << p.y.b.i[1] << "," << setw(16) << p.y.b.i[0] << endl;
}

#define EC_GEN_PRINT_TABLES

/*
 * 000 = -7
 * 001 = -5
 * 010 = -3
 * 011 = -1
 * 100 = 1
 * 101 = 3
 * 110 = 5
 * 111 = 7
 *
 * high bit = sign
 * low bits xor high bit = table entry
 */

// Verify generator multiplication tables are correct
static bool ec_gen_tables_comb_test() {
	ecpt_affine table[MG_n][MG_e];

	double s0 = m_clock.usec();
	u32 t0 = Clock::cycles();

	ufe t2b;
	ecpt base;
	ec_set(EC_G, base);

	for (int comb = 0; comb < MG_n; ++comb) {
		// Shift base to next table
		if (comb > 0) {
			for (int ii = 0; ii < MG_w; ++ii) {
				ec_dbl(base, base, false, t2b);
			}
		}

		ec_affine(base, table[comb][0]);

		ecpt two;
		ufe ignore;
		ec_dbl(base, two, false, ignore);
		fe_mul(t2b, two.t, two.t);

		ecpt q;
		ec_set(base, q);

		for (int u = 1; u < MG_e; ++u) {
			ec_add(q, two, q, false, false, false, t2b);

			ec_affine(q, table[comb][u]);
		}
	}

	u32 t1 = Clock::cycles();
	double s1 = m_clock.usec();

	cout << "Regenerated comb tables in " << (t1 - t0) << " cycles and " << (s1 - s0) << " usec" << endl;

#ifdef EC_GEN_PRINT_TABLES

	cout << dec << "static const u64 PRECOMP_TABLE_0[" << MG_n << "][8 * " << MG_e << "] = {";
	for (int jj = 0; jj < MG_n; ++jj) {
		cout << "{" << endl;
		ecpt_affine *ptr = &table[jj][0];
		for (int ii = 0; ii < MG_e; ++ii) {
			cout << "0x" << hex << ptr->x.a.i[0] << "ULL, 0x" << ptr->x.a.i[1] << "ULL, 0x" << ptr->x.b.i[0] << "ULL, 0x" << ptr->x.b.i[1] << "ULL," << endl;
			cout << "0x" << hex << ptr->y.a.i[0] << "ULL, 0x" << ptr->y.a.i[1] << "ULL, 0x" << ptr->y.b.i[0] << "ULL, 0x" << ptr->y.b.i[1] << "ULL," << endl;
			ptr++;
		}
		cout << "},";
	}
	cout << "};" << endl;

#endif

	for (int jj = 0; jj < MG_n; ++jj) {
		if (0 != memcmp(table[jj], GEN_TABLE[jj], sizeof(table[jj]))) {
			return false;
		}
	}

	return true;
}

static bool ec_gen_tables3_comb_test() {
	//const int t = 252;
	const int w = 8;
	//const int e = 32; // ceil(t / wv)
	const int d = 32; // e*v, v = 1

	ecpt_z1 table[128];

	const int ul = 1 << (w - 1);
	for (int u = 0; u < ul; ++u) {
		// P[u][v'] = 2^(ev') * (1 + u0*2^d + ... + u_(w-2)*2^((w-1)*d)) * P

		// q = u * P
		ufe t2b;
		ecpt q, s;

		ec_set(EC_G, q);

		for (int ii = 0; ii < (w - 1); ++ii) {
			if (u & (1 << ii)) {
				ec_set(EC_G, s);
				for (int jj = 0; jj < (d * (ii + 1)); ++jj) {
					ec_add(s, s, s, false, true, true, t2b);
				}
				ec_add(q, s, q, false, true, true, t2b);
			}
		}

		// Get affine version
		ecpt_affine qa;
		ec_affine(q, qa);

		// Store in table entry without Z=1
		table[u].x = qa.x;
		table[u].y = qa.y;
		fe_mul(qa.x, qa.y, table[u].t);
	}

#if 0

	ecpt_z1 *ptr = table;
	cout << "static const u64 PRECOMP_TABLE_3[12*128] = {" << endl;
	for (int ii = 0; ii < 128; ++ii) {
		cout << "0x" << hex << ptr->x.a.i[0] << "ULL, 0x" << ptr->x.a.i[1] << "ULL, 0x" << ptr->x.b.i[0] << "ULL, 0x" << ptr->x.b.i[1] << "ULL," << endl;
		cout << "0x" << hex << ptr->y.a.i[0] << "ULL, 0x" << ptr->y.a.i[1] << "ULL, 0x" << ptr->y.b.i[0] << "ULL, 0x" << ptr->y.b.i[1] << "ULL," << endl;
		cout << "0x" << hex << ptr->t.a.i[0] << "ULL, 0x" << ptr->t.a.i[1] << "ULL, 0x" << ptr->t.b.i[0] << "ULL, 0x" << ptr->t.b.i[1] << "ULL," << endl;
		ptr++;
	}
	cout << "};" << endl;

#endif

	if (0 != memcmp(table, SIMUL_GEN_TABLE, sizeof(table))) {
		return false;
	}

	return true;
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

	ecpt_affine q;
	ec_affine(a1, q);
	ec_expand(q, a1);

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

	u32 lsb = ec_recode_scalars_2(a1, b1, 128);

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

	ec_table_select_2(table, a1, b1, 0, true, r);
	if (!ec_isequal(r, c)) {
		cout << a << ", " << b << endl;
		return false;
	}
	return true;
}

static bool ec_table_select_2_test() {
	ecpt p1, p2;
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

//// Reference implementations for comparison

static bool ec_mul_ref(const u64 k[4], const ecpt_affine &P0, ecpt_affine &R) {
	ufe t2b;

	// Verify that curve has the specified order

	bool seen = false;

	ecpt p, g;
	ec_expand(P0, g);

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
				ec_expand(P0, p);
				seen = true;
			}
		}
	}

	ec_dbl(p, p, false, t2b);
	ec_dbl(p, p, false, t2b);

	ec_affine(p, R);

	return true;
}

static bool ec_simul_ref(const u64 k1[4], const ecpt_affine &P0, const u64 k2[4], const ecpt_affine &Q0, ecpt_affine &R) {
	ecpt_affine Pr, Qr;

	ec_mul_ref(k1, P0, Pr);
	ec_mul_ref(k2, Q0, Qr);

	ecpt P1, Q1;
	ec_expand(Pr, P1);
	ec_expand(Qr, Q1);

	ufe t2b;
	ecpt r;
	ec_add(P1, Q1, r, true, true, false, t2b);

	ec_affine(r, R);

	return true;
}


//// Test Driver

bool ec_mul_gen_test() {
	u64 k[4] = {0};
	ecpt_affine R1, R2;
	u8 a1[64], a2[64];

	vector<u32> t;
	double wall = 0;

	for (int jj = 0; jj < 10000; ++jj) {
		random_k(k);
		ec_mask_scalar(k);

		ec_mul_ref(k, EC_G_AFFINE, R1);

		double s0 = m_clock.usec();
		u32 t0 = Clock::cycles();

		ec_mul_gen(k, R2);

		u32 t1 = Clock::cycles();
		double s1 = m_clock.usec();

		t.push_back(t1 - t0);
		wall += s1 - s0;

		ecpt temp;
		ec_expand(R2, temp);
		ufe t2b;
		ec_dbl(temp, temp, true, t2b);
		ec_dbl(temp, temp, false, t2b);
		ec_affine(temp, R2);

		ec_save_xy(R1, a1);
		ec_save_xy(R2, a2);

		for (int ii = 0; ii < 64; ++ii) {
			if (a1[ii] != a2[ii]) {
				return false;
			}
		}
	}

	u32 median = quick_select(&t[0], (int)t.size());
	wall /= t.size();

	cout << "+ ec_mul_gen: `" << dec << median << "` median cycles, `" << wall << "` avg usec" << endl;

	return true;
}

bool ec_mul_test() {
	u64 k[4] = {0};
	ecpt_affine R1, R2;
	ecpt_affine BP;
	u8 a1[64], a2[64];

	random_k(k);
	ec_mul_ref(k, EC_G_AFFINE, BP);

	vector<u32> t;
	double wall = 0;

	for (int jj = 0; jj < 10000; ++jj) {
		random_k(k);
		ec_mask_scalar(k);

		ec_mul_ref(k, BP, R1);

		double s0 = m_clock.usec();
		u32 t0 = Clock::cycles();

		ec_mul(k, BP, R2);

		u32 t1 = Clock::cycles();
		double s1 = m_clock.usec();

		t.push_back(t1 - t0);
		wall += s1 - s0;

		ec_save_xy(R1, a1);
		ec_save_xy(R2, a2);

		for (int ii = 0; ii < 64; ++ii) {
			if (a1[ii] != a2[ii]) {
				return false;
			}
		}
	}

	u32 median = quick_select(&t[0], (int)t.size());
	wall /= t.size();

	cout << "+ ec_mul: `" << dec << median << "` median cycles, `" << wall << "` avg usec" << endl;

	return true;
}

bool ec_simul_test() {
	u64 k1[4] = {0};
	u64 k2[4] = {0};
	ecpt_affine B1, B2;
	ecpt_affine R1, R2;
	u8 a1[64], a2[64];

	random_k(k1);
	random_k(k2);
	ec_mul_ref(k1, EC_G_AFFINE, B1);
	ec_mul_ref(k2, EC_G_AFFINE, B2);

	vector<u32> t;
	double wall = 0;

	for (int jj = 0; jj < 10000; ++jj) {
		random_k(k1);
		random_k(k2);
		ec_mask_scalar(k1);
		ec_mask_scalar(k2);

		ec_simul_ref(k1, B1, k2, B2, R1);

		double s0 = m_clock.usec();
		u32 t0 = Clock::cycles();

		ec_simul(k1, B1, k2, B2, R2);

		u32 t1 = Clock::cycles();
		double s1 = m_clock.usec();

		t.push_back(t1 - t0);
		wall += s1 - s0;

		ec_save_xy(R1, a1);
		ec_save_xy(R2, a2);

		for (int ii = 0; ii < 64; ++ii) {
			if (a1[ii] != a2[ii]) {
				return false;
			}
		}
	}

	u32 median = quick_select(&t[0], (int)t.size());
	wall /= t.size();

	cout << "+ ec_simul: `" << dec << median << "` median cycles, `" << wall << "` avg usec" << endl;

	return true;
}

bool ec_simul_gen_test() {
	u64 k1[4] = {0};
	u64 k2[4] = {0};
	ecpt_affine R1, R2;
	ecpt_affine BP;
	u8 a1[64], a2[64];

	random_k(k2);
	ec_mul_ref(k2, EC_G_AFFINE, BP);

	vector<u32> t;
	double wall = 0;

	for (int jj = 0; jj < 10000; ++jj) {
		random_k(k1);
		random_k(k2);
		ec_mask_scalar(k1);
		ec_mask_scalar(k2);

		ec_simul_ref(k1, EC_G_AFFINE, k2, BP, R1);

		double s0 = m_clock.usec();
		u32 t0 = Clock::cycles();

		ec_simul_gen(k1, k2, BP, R2);

		u32 t1 = Clock::cycles();
		double s1 = m_clock.usec();

		t.push_back(t1 - t0);
		wall += s1 - s0;

		ec_save_xy(R1, a1);
		ec_save_xy(R2, a2);

		for (int ii = 0; ii < 64; ++ii) {
			if (a1[ii] != a2[ii]) {
				return false;
			}
		}
	}

	u32 median = quick_select(&t[0], (int)t.size());
	wall /= t.size();

	cout << "+ ec_simul_gen: `" << dec << median << "` median cycles, `" << wall << "` avg usec" << endl;

	return true;
}

bool mod_q_test() {
	u64 x[8], r[4];

	x[0] = 0xffffffffffffffffULL;
	x[1] = 0xffffffffffffffffULL;
	x[2] = 0xffffffffffffffffULL;
	x[3] = 0xffffffffffffffffULL;
	x[4] = 0xffffffffffffffffULL;
	x[5] = 0xffffffffffffffffULL;
	x[6] = 0xffffffffffffffffULL;
	x[7] = 0xffffffffffffffffULL;

	mod_q(x, r);

	if (r[0] != 0x72A7E6A3F7A11C27ULL ||
		r[1] != 0xA52B0BE10884939EULL ||
		r[2] != 0x95EB7B0E1A988566ULL ||
		r[3] != 0x093F8B602171C88EULL) {
		return false;
	}

	return true;
}

bool mul_mod_q_test() {
	u64 x[4], y[4], z[4], r[4];

	x[0] = 0xFB8A86C9E6022515ULL;
	x[1] = 0xD97FE1124FD8CC92ULL;
	x[2] = 0x782777E7572BA130ULL;
	x[3] = 0x0A64E21CF80B9B64ULL;
	y[0] = 0xEC7442A2DDA82CE0ULL;
	y[1] = 0x85F16DA062E80241ULL;
	y[2] = 0x21309454C67D3636ULL;
	y[3] = 0xE9296E5F048E01CCULL;
	z[0] = 0x140A07B4AD54B996ULL;
	z[1] = 0x5B73600FD51C45CDULL;
	z[2] = 0xC83C13EF9A0A3AC3ULL;
	z[3] = 0x003445C52BC607CFULL;

	mul_mod_q(x, y, z, r);

	if (r[0] != 0x9A5FC58C4E29F36EULL ||
		r[1] != 0x0A03DAB8CF16D699ULL ||
		r[2] != 0x6F161E3B5D31BBCEULL ||
		r[3] != 0x063D680741CBB9A1ULL) {
		return false;
	}

	x[0] = 0xffffffffffffffffULL;
	x[1] = 0xffffffffffffffffULL;
	x[2] = 0xffffffffffffffffULL;
	x[3] = 0xffffffffffffffffULL;
	y[0] = EC_Q[0] - 1;
	y[1] = EC_Q[1];
	y[2] = EC_Q[2];
	y[3] = EC_Q[3];
	z[0] = EC_Q[0] - 1;
	z[1] = EC_Q[1];
	z[2] = EC_Q[2];
	z[3] = EC_Q[3];

	mul_mod_q(x, y, z, r);

	if (r[0] != 0xB851F71EBA7E1BF5ULL ||
		r[1] != 0x08875560CEA50510ULL ||
		r[2] != 0xFFFFFFFFFFFFFFFAULL ||
		r[3] != 0x0FFFFFFFFFFFFFFFULL) {
		return false;
	}

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

	// Note: this fails for GHz > 4
	cout << "RDTSC instruction runs at " << (c - c0)/(t - t0)/1000.0 << " GHz" << endl;
}

static bool ec_recode_comb_test(const u64 x[4]) {
	u64 y[4];
	ec_recode_comb(x, y);

	// Follow the recoded bits to reconstruct the original scalars
	u64 z[4] = {0};

	for (int ii = 251; ii >= 0; ii--) {
		u32 u = (y[ii >> 6] >> (ii & 63)) & 1;

		z[3] = (z[2] >> 63) | (z[3] << 1);
		z[2] = (z[1] >> 63) | (z[2] << 1);
		z[1] = (z[0] >> 63) | (z[1] << 1);
		z[0] = z[0] << 1;

		if (u) {
			u128 sum = (u128)z[0] + 1;
			z[0] = (u64)sum;
			sum = (u128)z[1] + (u64)(sum >> 64);
			z[1] = (u64)sum;
			sum = (u128)z[2] + (u64)(sum >> 64);
			z[2] = (u64)sum;
			sum = (u128)z[3] + (u64)(sum >> 64);
			z[3] = (u64)sum;
		} else {
			s128 diff = (s128)z[0] - 1;
			z[0] = (u64)diff;
			diff = (diff >> 64) + z[1];
			z[1] = (u64)diff;
			diff = (diff >> 64) + z[2];
			z[2] = (u64)diff;
			diff = (diff >> 64) + z[3];
			z[3] = (u64)diff;
		}
	}

	u128 sum = (u128)z[0] + EC_Q[0];
	z[0] = (u64)sum;
	sum = (u128)z[1] + EC_Q[1] + (u64)(sum >> 64);
	z[1] = (u64)sum;
	sum = (u128)z[2] + EC_Q[2] + (u64)(sum >> 64);
	z[2] = (u64)sum;
	sum = (u128)z[3] + EC_Q[3] + (u64)(sum >> 64);
	z[3] = (u64)sum;

	z[3] &= 0x0fffffffffffffffULL;

	if ((x[0] & 1) == 0) {
		if (z[0] != x[0] || z[1] != x[1] ||
			z[2] != x[2] || z[3] != x[3]) {
			cout << "fails in even case";
			return false;
		}
		return true;
	}

	sum = (u128)z[0] + EC_Q[0];
	z[0] = (u64)sum;
	sum = (u128)z[1] + EC_Q[1] + (u64)(sum >> 64);
	z[1] = (u64)sum;
	sum = (u128)z[2] + EC_Q[2] + (u64)(sum >> 64);
	z[2] = (u64)sum;
	sum = (u128)z[3] + EC_Q[3] + (u64)(sum >> 64);
	z[3] = (u64)sum;

	z[3] &= 0x0fffffffffffffffULL;

	if (z[0] != x[0] || z[1] != x[1] ||
		z[2] != x[2] || z[3] != x[3]) {
		cout << "fails in odd case";
		return false;
	}
	return true;
}

int main() {
	cout << "Snowshoe Unit Tester: EC Scalar Multiplication" << endl;

#ifdef CAT_HAS_VECTOR_EXTENSIONS
	cout << "Using vector extensions for table lookups! <3" << endl;
#endif

	srand(0);

	m_clock.OnInitialize();

	tscTime();

	// Verify tables have not been tampered with
	assert(ec_gen_tables3_comb_test());
	assert(ec_gen_tables_comb_test());

	for (int ii = 0; ii < 1000; ++ii) {
		u64 k[4];
		random_k(k);
		ec_mask_scalar(k);
		assert(ec_recode_comb_test(k));
	}

	assert(mod_q_test());

	assert(mul_mod_q_test());

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

	assert(ec_mul_test());
	assert(ec_mul_gen_test());
	assert(ec_simul_test());
	assert(ec_simul_gen_test());

	cout << "All tests passed successfully." << endl;

	m_clock.OnFinalize();

	return 0;
}

