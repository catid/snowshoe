#include <iostream>
#include <cassert>
using namespace std;

// Math library
#include "../snowshoe/endo.cpp"

// static void gls_decompose(const u64 k[4], s32 &k1sign, u64 k1[2], s32 &k2sign, u64 k2[2])

static const u64 TEST_K0[4] = {
	0x91BB95B26470B944ULL,
	0x186A2F1F33217F72ULL,
	0xA058974AD3C6F3CDULL,
	0x0399805098D7D56FULL
};

static const u64 TEST_K1[4] = {
	0xCE9B68E3B09E01A4ULL,
	0xA6261414C0DC87D3ULL,
	0xFFFFFFFFFFFFFFFFULL,
	0x0FFFFFFFFFFFFFFFULL
};

static const u64 TEST_K2[4] = {
	0x679DFE17D6AC412FULL,
	0x43F1C74EDC9DC196ULL,
	0xA8A8D98EDB18E410ULL,
	0x0985EE47C6F67E9EULL
};

static const u64 TEST_K3[4] = {
	0xCE3469C57A30173EULL,
	0x5F6CD48A0AFFA60FULL,
	0x2519EDB7B96F26B1ULL,
	0x0B4AD868CD1641ACULL
};

static void gls_decompose_test() {
	s32 k1sign, k2sign;
	ufp k1, k2;

	gls_decompose(TEST_K1, k1sign, k1, k2sign, k2);

	assert(k1sign == 1);
	assert(k1.i[0] == 1);
	assert(k1.i[1] == 0);
	assert(k2sign == 0);
	assert(k2.i[0] == 0);
	assert(k2.i[1] == 0);

	gls_decompose(TEST_K0, k1sign, k1, k2sign, k2);

	assert(k1sign == 0);
	assert(k1.i[0] == 0xC14AABE9E079D148ULL);
	assert(k1.i[1] == 0x1E3B0E8CE06C74E5ULL);
	assert(k2sign == 1);
	assert(k2.i[0] == 0x680445984E433D40ULL);
	assert(k2.i[1] == 0x170475967D197366ULL);

	gls_decompose(TEST_K2, k1sign, k1, k2sign, k2);

	assert(k1sign == 0);
	assert(k1.i[0] == 0xC7620B2B8C69B128ULL);
	assert(k1.i[1] == 0x1354C079D167C5BCULL);
	assert(k2sign == 1);
	assert(k2.i[0] == 0x132501035CC11F8EULL);
	assert(k2.i[1] == 0x12BCB74AF1B58892ULL);

	gls_decompose(TEST_K3, k1sign, k1, k2sign, k2);

	assert(k1sign == 1);
	assert(k1.i[0] == 0xD8236D4762C9CD88ULL);
	assert(k1.i[1] == 0x0E546BB9D4D29156ULL);
	assert(k2sign == 0);
	assert(k2.i[0] == 0xE6809F829E581646ULL);
	assert(k2.i[1] == 0x00A1D93A9F379601ULL);
}


//// Entrypoint

int main() {
	cout << "Snowshoe Unit Tester: Scalar Decomposition" << endl;

	gls_decompose_test();

	cout << "All tests passed successfully." << endl;

	return 0;
}


