#include <iostream>
#include <cassert>
using namespace std;

// Math library
#include "../src/fp.inc"


//// Test Driver

static const ufp C0 = {
	{0, 0}
};

static const ufp C1 = {
	{1, 0}
};

static const ufp C2 = {
	{2, 0}
};

static const ufp C0F = {
	{0xffffffffffffffffULL, 0}
};

static const ufp C64 = {
	{0, 1}
};

static const ufp C65 = {
	{0, 2}
};

static const ufp CN2 = {
	{0xfffffffffffffffdULL, 0x7fffffffffffffffULL}
};

static const ufp CN1 = {
	{0xfffffffffffffffeULL, 0x7fffffffffffffffULL}
};

static const ufp CP = {
	{0xffffffffffffffffULL, 0x7fffffffffffffffULL}
};

static const ufp CP1 = {
	{0, 0x8000000000000000ULL}
};

static const ufp CF0 = {
	{0, 0xffffffffffffffffULL}
};

static const ufp CFF = {
	{0xffffffffffffffffULL, 0xffffffffffffffffULL}
};

static const ufp CR1 = {
	{0x09744238EF199911ULL, 0x6541AA8FCD8C4C65ULL}
};

static const ufp CR2 = {
	{0xD204049593D4A1D1ULL, 0x5281A3886F35ED6FULL}
};

static const ufp CX3 = {
	{0xB766E7802FB7635FULL, 0x3F42AC9208EEFF87ULL}
};

void fp_print(const ufp &x) {
	cout << "Value(H:L) = " << hex << x.i[1] << " : " << x.i[0] << endl;
}

static bool fp_isequal_test(const ufp &a, const ufp &b) {
	bool ct = fp_isequal_ct(a, b);
	bool vt = fp_isequal_vartime(a, b);

	assert(ct == vt);

	return ct;
}

bool fp_complete_reduce_test(const ufp &a, const ufp &expected) {
	ufp x;
	fp_set(a, x);
	fp_complete_reduce(x);
	return fp_isequal_test(x, expected);
}

bool fp_add_test(const ufp &a, const ufp &b, const ufp &expected) {
	ufp x, y, z;

	fp_set(a, x);
	fp_set(b, y);
	fp_add(x, y, z);
	fp_complete_reduce(z);
	if (!fp_isequal_test(z, expected)) {
		return false;
	}

	fp_set(a, x);
	fp_set(b, y);
	fp_add(x, y, x);
	fp_complete_reduce(x);
	if (!fp_isequal_test(x, expected)) {
		return false;
	}

	fp_set(a, x);
	fp_set(b, y);
	fp_add(x, y, y);
	fp_complete_reduce(y);
	if (!fp_isequal_test(y, expected)) {
		return false;
	}

	return true;
}

bool fp_add_smallk_test(const ufp &a, const u32 b, const ufp &expected) {
	ufp x, z;

	fp_set(a, x);
	fp_add_smallk(x, b, z);
	fp_complete_reduce(z);
	if (!fp_isequal_test(z, expected)) {
		return false;
	}

	fp_set(a, x);
	fp_add_smallk(x, b, x);
	fp_complete_reduce(x);
	if (!fp_isequal_test(x, expected)) {
		return false;
	}

	return true;
}

bool fp_sub_test(const ufp &a, const ufp &b, const ufp &expected) {
	ufp x, y, z;

	fp_set(a, x);
	fp_set(b, y);
	fp_sub(x, y, z);
	fp_complete_reduce(z);
	if (!fp_isequal_test(z, expected)) {
		return false;
	}

	fp_set(a, x);
	fp_set(b, y);
	fp_sub(x, y, x);
	fp_complete_reduce(x);
	if (!fp_isequal_test(x, expected)) {
		return false;
	}

	fp_set(a, x);
	fp_set(b, y);
	fp_sub(x, y, y);
	fp_complete_reduce(y);
	if (!fp_isequal_test(y, expected)) {
		return false;
	}

	return true;
}

bool fp_neg_test(const ufp &a, const ufp &expected) {
	ufp x, z;

	fp_set(a, x);
	fp_neg(x, z);
	if (!fp_isequal_test(z, expected)) {
		return false;
	}

	fp_set(a, x);
	fp_neg(x, x);
	if (!fp_isequal_test(x, expected)) {
		return false;
	}

	return true;
}

bool fp_set_smallk_test(const u32 a, const ufp &expected) {
	ufp x;
	fp_set_smallk(a, x);
	return fp_isequal_test(x, expected);
}

bool fp_mul_small_test(const ufp &a, const u32 b) {
	ufp x, y, z, w;

	fp_set(a, x);
	fp_set_smallk(b, y);
	fp_mul(x, y, z);
	fp_mul_smallk(x, b, w);
	fp_complete_reduce(z);
	fp_complete_reduce(w);
	if (!fp_isequal_test(z, w)) {
		return false;
	}

	fp_set(a, x);
	fp_set_smallk(b, y);
	fp_mul_smallk(x, b, w);
	fp_mul(x, y, x);
	fp_complete_reduce(x);
	fp_complete_reduce(w);
	if (!fp_isequal_test(x, w)) {
		return false;
	}

	fp_set(a, x);
	fp_set_smallk(b, y);
	fp_mul_smallk(x, b, w);
	fp_mul(x, y, y);
	fp_complete_reduce(y);
	fp_complete_reduce(w);
	if (!fp_isequal_test(y, w)) {
		return false;
	}

	return true;
}

bool fp_mul_sqr_test(const ufp &a) {
	ufp x, y, z, w;

	fp_set(a, x);
	fp_set(a, y);
	fp_mul(x, y, z);
	fp_sqr(x, w);
	fp_complete_reduce(z);
	fp_complete_reduce(w);
	if (!fp_isequal_test(z, w)) {
		return false;
	}

	fp_set(a, x);
	fp_set(a, y);
	fp_sqr(x, w);
	fp_mul(x, y, x);
	fp_complete_reduce(x);
	fp_complete_reduce(w);
	if (!fp_isequal_test(x, w)) {
		return false;
	}

	fp_set(a, x);
	fp_set(a, y);
	fp_sqr(x, w);
	fp_mul(x, y, y);
	fp_complete_reduce(y);
	fp_complete_reduce(w);
	if (!fp_isequal_test(y, w)) {
		return false;
	}

	fp_set(a, x);
	fp_set(a, y);
	fp_mul(x, y, z);
	fp_sqr(x, x);
	fp_complete_reduce(z);
	fp_complete_reduce(x);
	if (!fp_isequal_test(z, x)) {
		return false;
	}

	return true;
}

bool fp_mul_inv_test(const ufp &a, const ufp &expected) {
	ufp x, y, z;

	fp_set(a, x);
	fp_inv(x, y);
	fp_mul(x, y, z);
	if (!fp_isequal_test(z, expected)) {
		return false;
	}

	fp_set(a, x);
	fp_inv(x, x);
	fp_set(a, y);
	fp_mul(x, y, z);
	if (!fp_isequal_test(z, expected)) {
		return false;
	}

	return true;
}

bool fp_save_load_test(const ufp &a) {
	ufp x, y;
	u8 buffer[17] = {0};
	fp_set(a, x);
	fp_save(x, buffer);
	fp_load(buffer, y);
	if (!fp_isequal_test(x, y)) {
		return false;
	}
	if (buffer[16] != 0) {
		return false;
	}
	return true;
}

bool fp_exp_test(const ufp &x, const ufp &e, const ufp &expected) {
	ufp r;
	fp_set_smallk(1, r);

	bool seen = false;

	for (int ii = 126; ii >= 0; --ii) {
		if (seen) {
			fp_sqr(r, r);
		}

		if ((e.w >> ii) & 1) {
			fp_mul(r, x, r);
			seen = true;
		}
	}

	return fp_isequal_test(r, expected);
}

bool fp_exp_inv_test(const ufp &a) {
	ufp x, y;

	fp_set(a, x);
	fp_inv(x, y);
	return fp_exp_test(x, CN2, y);
}

bool fp_mul_test(const ufp &a, const ufp &b, const ufp &expected) {
	ufp x, y, z;

	fp_set(a, x);
	fp_set(b, y);
	fp_mul(x, y, z);

	return fp_isequal_test(z, expected);
}

bool fp_set_mask_test() {
	ufp a, r;

	u128 zero = 0;
	u128 neg = ~(u128)0;

	fp_set(CR1, a);
	fp_set(CR2, r);

	// Verify that passing a -1 mask will set r,
	// and passing 0 will leave r unaffected

	fp_set_mask(a, zero, r);

	if (!fp_isequal_test(r, CR2)) {
		return false;
	}

	fp_set_mask(a, neg, r);

	if (!fp_isequal_test(r, CR1)) {
		return false;
	}

	return true;
}

bool fp_xor_mask_test() {
	ufp a, r;

	u128 zero = 0;
	u128 neg = ~(u128)0;

	fp_set(CR1, a);
	fp_zero(r);

	// Verify that passing a -1 mask will xor r,
	// and passing 0 will leave r unaffected

	fp_xor_mask(a, zero, r);

	if (!fp_iszero_vartime(r)) {
		return false;
	}

	fp_xor_mask(a, neg, r);

	if (!fp_isequal_test(r, CR1)) {
		return false;
	}

	return true;
}

bool fp_neg_mask_test() {
	ufp a, r;

	u128 zero = 0;
	u128 neg = ~(u128)0;

	// Verify that passing a -1 mask will leave r unaffected,
	// and passing 0 will negate r in place

	fp_neg_mask(neg, CR1, r);

	fp_set(CR1, a);
	fp_neg(a, a);

	if (fp_isequal_test(r, CR1)) {
		return false;
	}

	if (!fp_isequal_test(r, a)) {
		return false;
	}

	fp_neg_mask(neg, r, r);

	if (!fp_isequal_test(r, CR1)) {
		return false;
	}

	fp_neg_mask(zero, r, r);

	if (!fp_isequal_test(r, CR1)) {
		return false;
	}

	return true;
}

bool fp_div2_test(const ufp &a) {
	ufp a2, a1;

	fp_mul(a, C2, a2);
	fp_div2(a2, a1);
	fp_complete_reduce(a1);

	return fp_isequal_test(a1, a);
}

bool fp_chi_test() {
	ufp a;

	fp_set(CN1, a);
	if (fp_chi(a) != -1) {
		cout << "chi failed -1" << endl;
		return false;
	}
	fp_set(C0, a);
	if (fp_chi(a) != 0) {
		cout << "chi failed 0" << endl;
		return false;
	}
	fp_set(C1, a);
	if (fp_chi(a) != 1) {
		cout << "chi failed 1" << endl;
		return false;
	}
	fp_set(C2, a);
	if (fp_chi(a) != 1) {
		cout << "chi failed 2" << endl;
		return false;
	}
	fp_set(CR1, a);
	if (fp_chi(a) != 1) {
		cout << "chi failed CR1" << endl;
		return false;
	}
	fp_set(CR2, a);
	if (fp_chi(a) != -1) {
		cout << "chi failed CR2" << endl;
		return false;
	}
	fp_set(CX3, a);
	if (fp_chi(a) != 1) {
		cout << "chi failed CX3" << endl;
		return false;
	}

	return true;
}

bool fp_sqrt_test(const ufp &a) {
	ufp a2, a1;

	int chi = fp_chi(a);

	fp_mul(a, a, a2);
	fp_sqrt(a2, a1);

	if (chi == -1) {
		fp_neg(a1, a1);
	}

	fp_complete_reduce(a1);

	return fp_isequal_test(a1, a);
}


//// Entrypoint

int main() {
	cout << "Snowshoe Unit Tester: Fp base finite field arithmetic" << endl;

	// fp_iszero_vartime:
	assert(fp_iszero_vartime(C0));
	assert(fp_iszero_vartime(CP));
	assert(!fp_iszero_vartime(CP1));
	assert(!fp_iszero_vartime(CN1));
	assert(!fp_iszero_vartime(CFF));

	// fp_infield_vartime:
	assert(!fp_infield_vartime(CFF));
	assert(!fp_infield_vartime(CF0));
	assert(!fp_infield_vartime(CP1));
	assert(!fp_infield_vartime(CP));
	assert(fp_infield_vartime(CN1));
	assert(fp_infield_vartime(C64));
	assert(fp_infield_vartime(C65));
	assert(fp_infield_vartime(C0F));
	assert(fp_infield_vartime(C0));

	// fp_set, fp_neg:
	assert(fp_neg_test(C0, C0));
	assert(fp_neg_test(C1, CN1));
	assert(fp_neg_test(CP, C0));
	assert(fp_neg_test(CN1, C1));

	// fp_set, fp_add_smallk:
	assert(fp_add_smallk_test(C0, 0, C0));
	assert(fp_add_smallk_test(C0, 1, C1));
	assert(fp_add_smallk_test(C1, 1, C2));
	assert(fp_add_smallk_test(CN1, 1, C0));
	assert(fp_add_smallk_test(C0, 2, C2));
	assert(fp_add_smallk_test(CN1, 2, C1));
	assert(fp_add_smallk_test(CN1, 3, C2));
	assert(fp_add_smallk_test(C0F, 1, C64));

	// fp_set_smallk:
	assert(fp_set_smallk_test(0, C0));
	assert(fp_set_smallk_test(1, C1));
	assert(fp_set_smallk_test(2, C2));

	// fp_zero:
	ufp x;
	fp_set(CFF, x);
	fp_zero(x);
	assert(fp_isequal_test(x, C0));
	assert(fp_iszero_vartime(x));

	// fp_set, fp_complete_reduce:
	assert(fp_complete_reduce_test(C0, C0));
	assert(fp_complete_reduce_test(C1, C1));
	assert(fp_complete_reduce_test(C64, C64));
	assert(fp_complete_reduce_test(C65, C65));
	assert(fp_complete_reduce_test(CN1, CN1));
	assert(fp_complete_reduce_test(CP, C0));

	// fp_set, fp_add, fp_complete_reduce: (infield + infield ?= expected)
	assert(fp_add_test(C0, C1, C1));
	assert(fp_add_test(C1, C1, C2));
	assert(fp_add_test(C0F, C1, C64));
	assert(fp_add_test(C0, C64, C64));
	assert(fp_add_test(C64, C64, C65));
	assert(fp_add_test(CN1, C0, CN1));
	assert(fp_add_test(CN1, C1, C0));
	assert(fp_add_test(CN1, C2, C1));

	// fp_set, fp_sub, fp_complete_reduce: (infield - infield ?= expected)
	assert(fp_sub_test(C2, C1, C1));
	assert(fp_sub_test(C1, C2, CN1));
	assert(fp_sub_test(C0, C1, CN1));
	assert(fp_sub_test(C1, C1, C0));
	assert(fp_sub_test(C1, C0, C1));
	assert(fp_sub_test(C65, C64, C64));
	assert(fp_sub_test(C65, C65, C0));
	assert(fp_sub_test(C64, C0F, C1));
	assert(fp_sub_test(C0F, C64, CN1));

	// fp_mul
	assert(fp_mul_test(C64, C2, C65));
	assert(fp_mul_test(C0, C1, C0));
	assert(fp_mul_test(C1, C2, C2));
	assert(fp_mul_test(C1, C1, C1));

	// fp_mul_smallk <-> fp_mul:
	assert(fp_mul_small_test(C0, 0));
	assert(fp_mul_small_test(C0, 1));
	assert(fp_mul_small_test(C1, 0));
	assert(fp_mul_small_test(C1, 1));
	assert(fp_mul_small_test(C1, 2));
	assert(fp_mul_small_test(C2, 0));
	assert(fp_mul_small_test(C2, 1));
	assert(fp_mul_small_test(C0F, 109));
	assert(fp_mul_small_test(C64, 109));
	assert(fp_mul_small_test(C65, 109));
	assert(fp_mul_small_test(CN1, 109));
	assert(fp_mul_small_test(CP, 109));
	assert(fp_mul_small_test(C0, 0xffffffff));
	assert(fp_mul_small_test(C1, 0xffffffff));
	assert(fp_mul_small_test(C2, 0xffffffff));
	assert(fp_mul_small_test(C0F, 0xffffffff));
	assert(fp_mul_small_test(C64, 0xffffffff));
	assert(fp_mul_small_test(C65, 0xffffffff));
	assert(fp_mul_small_test(CN1, 0xffffffff));
	assert(fp_mul_small_test(CP, 0xffffffff));

	// fp_mul <-> fp_sqr:
	assert(fp_mul_sqr_test(C0));
	assert(fp_mul_sqr_test(C1));
	assert(fp_mul_sqr_test(C2));
	assert(fp_mul_sqr_test(C0F));
	assert(fp_mul_sqr_test(C64));
	assert(fp_mul_sqr_test(C65));
	assert(fp_mul_sqr_test(CN1));
	assert(fp_mul_sqr_test(CP));

	// fp_mul <-> fp_inv:
	assert(fp_mul_inv_test(C0, C0));
	assert(fp_mul_inv_test(C1, C1));
	assert(fp_mul_inv_test(C2, C1));
	assert(fp_mul_inv_test(C0F, C1));
	assert(fp_mul_inv_test(C64, C1));
	assert(fp_mul_inv_test(C65, C1));
	assert(fp_mul_inv_test(CN1, C1));
	assert(fp_mul_inv_test(CP, CP));

	// fp_save, fp_load:
	assert(fp_save_load_test(C0));
	assert(fp_save_load_test(C1));
	assert(fp_save_load_test(C2));
	assert(fp_save_load_test(C0F));
	assert(fp_save_load_test(C64));
	assert(fp_save_load_test(C65));
	assert(fp_save_load_test(CN1));

	// fp_mul, fp_sqr:
	assert(fp_exp_test(C0, C0, C1));
	assert(fp_exp_test(C64, C0, C1));
	assert(fp_exp_test(C65, C1, C65));
	assert(fp_exp_test(C1, C2, C1));
	assert(fp_exp_test(C0, C2, C0));
	assert(fp_exp_test(CR1, CR2, CX3));

	// fp_inv <-> fp_mul, fp_sqr:
	assert(fp_exp_inv_test(C1));
	assert(fp_exp_inv_test(C2));
	assert(fp_exp_inv_test(C0F));
	assert(fp_exp_inv_test(C64));
	assert(fp_exp_inv_test(C65));
	assert(fp_exp_inv_test(CN1));

	// fp_set_mask:
	assert(fp_set_mask_test());

	// fp_xor_mask:
	assert(fp_xor_mask_test());

	// fp_neg_mask:
	assert(fp_neg_mask_test());

	// fp_div2:
	assert(fp_div2_test(C0));
	assert(fp_div2_test(C1));
	assert(fp_div2_test(C2));
	assert(fp_div2_test(CN1));
	assert(fp_div2_test(C64));
	assert(fp_div2_test(C65));
	assert(fp_div2_test(C0F));
	assert(fp_div2_test(CR1));
	assert(fp_div2_test(CR2));
	assert(fp_div2_test(CX3));

	// fp_chi:
	assert(fp_chi_test());

	// fp_sqrt:
	assert(fp_sqrt_test(C0));
	assert(fp_sqrt_test(C1));
	assert(fp_sqrt_test(C2));
	assert(fp_sqrt_test(CN1));
	assert(fp_sqrt_test(C64));
	assert(fp_sqrt_test(C65));
	assert(fp_sqrt_test(C0F));
	assert(fp_sqrt_test(CR1));
	assert(fp_sqrt_test(CR2));
	assert(fp_sqrt_test(CX3));

	cout << "All tests passed successfully." << endl;

	return 0;
}

