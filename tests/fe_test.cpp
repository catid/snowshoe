#include <iostream>
#include <cassert>
using namespace std;

// Math library
#include "../src/fe.inc"

static const ufe C1 = {
	{{1, 0}}, // a (real part)
	{{0, 0}} // b (imaginary part)
};

/*
 * Using [6], c = a^b:
 *
 * a = 0x36E49794B5C0A0DE8BCD499E62DCFBB6*i + 0x5532CCBF3865DA90C05DFCCF2B1921F
 * b = 0x5A51F0A5F2B7F0D137E03ACCACC9D62C
 * c = 0x32F8EF809041D20EE31B7E4CFD07F57D*i + 0x4AC7C152C6E8D0D970472D154641FCA7
 */

/*
 * a * 0x7fdead41
 * = 47323413677645719676687429128767419487*i +
 *   29331494648517259664313532906541209288
 * = 0x239A286F3AF75086A445168838E0345F*i +
 *   0x16110A56B7F85C2AC5DBE085BC08D6C8
 */

static const ufe CRA = {
	{{0x0C05DFCCF2B1921FULL, 0x05532CCBF3865DA9ULL}}, // a (real part)
	{{0x8BCD499E62DCFBB6ULL, 0x36E49794B5C0A0DEULL}} // b (imaginary part)
};

static const u32 CSK = 0x7fdead41;

static const ufe CKP = {
	{{0xC5DBE085BC08D6C8ULL, 0x16110A56B7F85C2AULL}},
	{{0xA445168838E0345FULL, 0x239A286F3AF75086ULL}}
};

static const ufp CN2 = {
	{0xfffffffffffffffdULL, 0x7fffffffffffffffULL}
};

static const ufp CN1 = {
	{0xfffffffffffffffeULL, 0x7fffffffffffffffULL}
};

static const ufp CRB = {
	{0x37E03ACCACC9D62CULL, 0x5A51F0A5F2B7F0D1ULL}
};

static const ufe CXC = {
	{{0x70472D154641FCA7ULL, 0x4AC7C152C6E8D0D9ULL}}, // a (real part)
	{{0xE31B7E4CFD07F57DULL, 0x32F8EF809041D20EULL}} // b (imaginary part)
};

static const ufe CI = {
	{{0x29BE6E4878ADC1DAULL, 0x094B07FBA966C9DBULL}},
	{{0x3259A6EF6AAB1954ULL, 0x237C4B0DE91C3A16ULL}}
};

static const ufe CSR = {
	{{0xF2902B89B829926AULL, 0x4EDF9690D65F8F91ULL}},
	{{0xAF1FD5388F1BC28FULL, 0x078F791696C0DCCBULL}}
};

static const ufe COF1 = {
	{{0x8C05DFCCF2B1921FULL, 0x85532CCBF3865DA9ULL}}, // a (real part)
	{{0x8BCD499E62DCFBB6ULL, 0x36E49794B5C0A0DEULL}} // b (imaginary part)
};

static const ufe COF2 = {
	{{0x8C05DFCCF2B1921FULL, 0x05532CCBF3865DA9ULL}}, // a (real part)
	{{0x8BCD499E62DCFBB6ULL, 0x86E49794B5C0A0DEULL}} // b (imaginary part)
};

static const ufe COF3 = {
	{{0x8C05DFCCF2B1921FULL, 0x05532CCBF3865DA9ULL}}, // a (real part)
	{{0xffffffffffffffffULL, 0x7fffffffffffffffULL}} // b (imaginary part)
};

static const ufe COF4 = {
	{{0xffffffffffffffffULL, 0x7fffffffffffffffULL}},
	{{0x8C05DFCCF2B1921FULL, 0x05532CCBF3865DA9ULL}}
};

static const ufe CIF1 = {
	{{0xefffffffffffffffULL, 0x7fffffffffffffffULL}},
	{{0x8C05DFCCF2B1921FULL, 0x05532CCBF3865DA9ULL}}
};

static const ufe CIF2 = {
	{{0x8C05DFCCF2B1921FULL, 0x05532CCBF3865DA9ULL}},
	{{0xefffffffffffffffULL, 0x7fffffffffffffffULL}}
};

static const ufe COF5 = {
	{{0xffffffffffffffffULL, 0x7fffffffffffffffULL}},
	{{0xffffffffffffffffULL, 0x7fffffffffffffffULL}}
};

static const ufe C0 = {
	{{0, 0}},
	{{0, 0}}
};

static const ufe CU = {
	{{2, 0}},
	{{1, 0}}
};

static bool fe_isequal_test(const ufe &x, const ufe &y) {
	bool ct = fe_isequal_ct(x, y);
	bool vt = fe_isequal_vartime(x, y);

	assert(ct == vt);

	return ct;
}

void fe_print(const ufe &x) {
	cout << "Real(H:L) = " << hex << x.a.i[1] << " : " << x.a.i[0] << endl;
	cout << "Imag(H:L) = " << hex << x.b.i[1] << " : " << x.b.i[0] << endl;
}

bool fe_exp_test(const ufe &x, const ufp &e, const ufe &expected) {
	ufe r;
	fe_set_smallk(1, r);

	bool seen = false;

	for (int ii = 126; ii >= 0; --ii) {
		if (seen) {
			fe_sqr(r, r);
		}

		if ((e.i[ii/64] >> (ii%64)) & 1) {
			fe_mul(r, x, r);
			seen = true;
		}
	}

	fe_complete_reduce(r);
	return fe_isequal_test(r, expected);
}

bool fe_inv_test(const ufe &a, const ufe &expected) {
	ufe x, z;

	fe_set(a, x);
	fe_inv(x, z);
	if (!fe_isequal_test(z, expected)) {
		return false;
	}

	fe_set(a, x);
	fe_inv(x, x);
	fe_complete_reduce(x);
	if (!fe_isequal_test(x, expected)) {
		return false;
	}

	return true;
}

bool fe_mul_sqr_test(const ufe &a) {
	ufe x, y;

	fe_set(a, x);
	fe_mul(x, x, x);
	fe_set(a, y);
	fe_sqr(y, y);
	return fe_isequal_test(x, y);
}

bool fe_mul_inv_test(const ufe &a) {
	ufe x, y;

	fe_set(a, x);
	fe_inv(x, x);
	fe_set(a, y);
	fe_mul(x, y, x);

	fe_complete_reduce(x);
	return fe_isequal_test(x, C1);
}

bool fe_conj_test(const ufe &a) {
	ufe x;

	fe_set(a, x);
	fe_conj(x, x);
	if (fe_isequal_test(x, a)) {
		return false;
	}
	fe_conj(x, x);
	return fe_isequal_test(x, a);
}

bool fe_neg_test(const ufe &a) {
	ufe x;

	fe_set(a, x);
	fe_neg(x, x);
	if (fe_isequal_test(x, a)) {
		return false;
	}
	fe_neg(x, x);
	if (!fe_isequal_test(x, a)) {
		return false;
	}

	fe_set(a, x);
	fe_neg(x, x);
	fe_add(a, x, x);
	return fe_iszero_vartime(x);
}

bool fe_save_load_test(const ufe &a) {
	ufe x, y;
	u8 buffer[33] = {0};
	fe_set(a, x);
	fe_save(x, buffer);
	fe_load(buffer, y);
	if (!fe_isequal_test(x, y)) {
		return false;
	}
	if (buffer[32] != 0) {
		return false;
	}
	return true;
}

bool fe_complete_reduce_test() {
	ufe x;
	fe_set(COF5, x);
	fe_complete_reduce(x);
	return fe_isequal_test(x, C0);
}

bool fe_mulu_test(const ufe &a) {
	ufe x, y, z;
	fe_set(a, x);
	fe_set(CU, y);
	fe_mul(x, y, x);
	fe_mul_u(a, z);
	return fe_isequal_test(x, z);
}

bool fe_set_mask_test() {
	ufe a, r;

	static const u64 zero = 0, neg = ~(u64)0;

	fe_set(CXC, a);
	fe_set(CSR, r);

	// Verify that passing a -1 mask will set r,
	// and passing 0 will leave r unaffected

	fe_set_mask(a, zero, r);

	if (!fe_isequal_test(r, CSR)) {
		return false;
	}

	fe_set_mask(a, neg, r);

	if (!fe_isequal_test(r, CXC)) {
		return false;
	}

	return true;
}

bool fe_xor_mask_test() {
	ufe a, r;

	static const u64 zero = 0, neg = ~(u64)0;

	fe_set(CXC, a);
	fe_zero(r);

	// Verify that passing a -1 mask will xor r,
	// and passing 0 will leave r unaffected

	fe_xor_mask(a, zero, r);

	if (!fe_iszero_vartime(r)) {
		return false;
	}

	fe_xor_mask(a, neg, r);

	if (!fe_isequal_test(r, CXC)) {
		return false;
	}

	return true;
}

bool fe_neg_mask_test() {
	ufe a, r;

	static const u64 zero = 0, neg = ~(u64)0;

	// Verify that passing a -1 mask will leave r unaffected,
	// and passing 0 will negate r in place

	fe_neg_mask(neg, CXC, r);

	fe_set(CXC, a);
	fe_neg(a, a);

	if (fe_isequal_test(r, CXC)) {
		return false;
	}

	if (!fe_isequal_test(r, a)) {
		return false;
	}

	fe_neg_mask(neg, r, r);

	if (!fe_isequal_test(r, CXC)) {
		return false;
	}

	fe_neg_mask(zero, r, r);

	if (!fe_isequal_test(r, CXC)) {
		return false;
	}

	return true;
}

bool fe_mul_smallk_test() {
	ufe r;

	fe_set(CRA, r);
	fe_mul_smallk(r, CSK, r);

	return fe_isequal_test(r, CKP);
}

bool fe_add_set_smallk_test() {
	ufe a, b, r;

	fe_set_smallk(CSK, a);
	fe_zero(r);
	fe_add_smallk(r, CSK, b);
	fe_add(a, r, r);

	return fe_isequal_test(a, b) && fe_isequal_test(a, r);
}

static bool fe_sqrt_test(const ufe &a, bool expected_valid) {
	ufe a2, a1, a1n;

	fe_sqr(a, a2);
	bool valid = fe_sqrt(a2, a1, true);

	int chi = fe_chi(a2);
	if (chi == -1) {
		cout << "fe_chi falsely accused a quadratic residue of being a non-residue" << endl;
		return false;
	}

	fe_neg(a1, a1n);
	fe_complete_reduce(a1);
	fe_complete_reduce(a1n);

	if (expected_valid != valid) {
		cout << "fe_sqrt failed validation test" << endl;
		return false;
	}

	if (!valid) {
		return true;
	}

	// Result is +/- the input
	return fe_isequal_test(a1, a) || fe_isequal_test(a1n, a);
}


//// Entrypoint

int main() {
	cout << "Snowshoe Unit Tester: Fp^2 Optimal Extension Field Arithmetic" << endl;

	// fe_mul <-> fe_sqr:
	assert(fe_mul_sqr_test(CRA));
	assert(fe_mul_sqr_test(CXC));
	assert(fe_mul_sqr_test(CI));
	assert(fe_mul_sqr_test(CSR));

	// fe_mul <-> fe_inv:
	assert(fe_mul_inv_test(CRA));
	assert(fe_mul_inv_test(CI));
	assert(fe_mul_inv_test(CXC));
	assert(fe_mul_inv_test(CSR));

	// fe_mul, fe_sqr:
	assert(fe_exp_test(CRA, CRB, CXC));

	// fe_inv:
	assert(fe_inv_test(CRA, CI));

	// fe_conj:
	assert(fe_conj_test(CRA));
	assert(fe_conj_test(CXC));
	assert(fe_conj_test(CI));
	assert(fe_conj_test(CSR));

	// fe_infield_vartime:
	assert(!fe_infield_vartime(COF1));
	assert(!fe_infield_vartime(COF2));
	assert(!fe_infield_vartime(COF3));
	assert(!fe_infield_vartime(COF4));
	assert(fe_infield_vartime(CIF1));
	assert(fe_infield_vartime(CIF2));

	// fe_neg:
	assert(fe_neg_test(CRA));
	assert(fe_neg_test(CXC));
	assert(fe_neg_test(CI));
	assert(fe_neg_test(CSR));

	// fe_save, fe_load:
	assert(fe_save_load_test(CIF1));
	assert(fe_save_load_test(CIF2));
	assert(fe_save_load_test(CRA));
	assert(fe_save_load_test(CXC));
	assert(fe_save_load_test(CI));
	assert(fe_save_load_test(CSR));

	// fe_mul_u:
	assert(fe_mulu_test(CIF1));
	assert(fe_mulu_test(CIF2));
	assert(fe_mulu_test(CRA));
	assert(fe_mulu_test(CXC));
	assert(fe_mulu_test(CI));
	assert(fe_mulu_test(CSR));

	// fe_complete_reduce:
	assert(fe_complete_reduce_test());

	// fe_set_mask:
	assert(fe_set_mask_test());

	// fe_xor_mask:
	assert(fe_xor_mask_test());

	// fe_neg_mask:
	assert(fe_neg_mask_test());

	// fe_mul_smallk:
	assert(fe_mul_smallk_test());

	// fe_set_smallk <-> fe_add_smallk:
	assert(fe_add_set_smallk_test());

	// fe_sqrt:
	assert(fe_sqrt_test(CIF1, true));
	assert(fe_sqrt_test(CIF2, true));
	assert(fe_sqrt_test(CRA, true));
	assert(fe_sqrt_test(CXC, true));
	assert(fe_sqrt_test(CI, true));
	assert(fe_sqrt_test(CSR, true));

	cout << "All tests passed successfully." << endl;

	return 0;
}

