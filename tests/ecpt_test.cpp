#include <iostream>
#include <iomanip>
#include <cassert>
using namespace std;

// Math library
#include "../snowshoe/ecpt.inc"

static const u64 PRIME_ORDER[4] = {
	0xCE9B68E3B09E01A5ULL,
	0xA6261414C0DC87D3ULL,
	0xFFFFFFFFFFFFFFFFULL,
	0x0FFFFFFFFFFFFFFFULL
};

static const u64 RS1[4] = {
	0x752a766243247047ULL,
	0x4038437a4730512aULL,
	0x472c4f4a284c5e41ULL,
	0x053f5676464c4772ULL
};
static const u64 RS2[4] = {
	0xCE9B68E3B09E01A5ULL,
	0xA6261414C0DC87D3ULL,
	0xFFFFFFFFFFFFFFFFULL,
	0x0FFFFFFFFFFFFFFFULL // same as group order
};
static const u64 RS3[4] = {
	0x64CDBE5FF1299E8AULL,
	0x285A57A4F8FC117BULL,
	0x266159ED9FDD9427ULL,
	0x6F3A44DADC9FFB37ULL
};
static const u64 RS4[4] = {
	0x8867CDE452E85B4BULL,
	0xDEDB6B963522785EULL,
	0x8314224BEF434315ULL,
	0xAB8729EF444E56E2ULL
};
static const u64 RS5[4] = {
	0xFB38A9818E5CD667ULL,
	0x331C599C24F3DC55ULL,
	0x46111D989ACCA9B6ULL,
	0x75185153C9324A17ULL
};
static const u64 RS6[4] = {
	0x2EFAF7BB2959BA76ULL,
	0xB6BF25156CD86952ULL,
	0x38F83FD44F95B353ULL,
	0x529745EDDB1CB76AULL
};
static const u64 RS7[4] = {
	0x43AD6643822BB3FFULL,
	0x58D3BFF2B16DC3B6ULL,
	0xCAA2FA22858DAE47ULL,
	0x8D3E31ABB8A18242ULL
};
static const u64 RS8[4] = {
	0xEFDEBDC39C4FC996ULL,
	0x8F43652138EA399fULL,
	0xD6E9AC329AB7527EULL,
	0xB180000000E922E1ULL
};

static const ufe Cnn = {
	{{0xffffffffffffffffULL, 0x7fffffffffffffffULL}},
	{{0xffffffffffffffffULL, 0x7fffffffffffffffULL}}
};

static const ufe C0n = {
	{{0, 0}},
	{{0xffffffffffffffffULL, 0x7fffffffffffffffULL}}
};

static const ufe Cn0 = {
	{{0xffffffffffffffffULL, 0x7fffffffffffffffULL}},
	{{0, 0}}
};

static const ufe C00 = {
	{{0, 0}},
	{{0, 0}}
};

static const ufe GXn = {
	{{15}},
	{{0xffffffffffffffffULL, 0x7fffffffffffffffULL}}
};

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
	cout << "Z : = " << setw(16) << hex << p.z.a.i[1] << "," << setw(16) << p.z.a.i[0] << " + i * " << setw(16) << p.z.b.i[1] << "," << setw(16) << p.z.b.i[0] << endl;
	cout << "T : = " << setw(16) << hex << p.t.a.i[1] << "," << setw(16) << p.t.a.i[0] << " + i * " << setw(16) << p.t.b.i[1] << "," << setw(16) << p.t.b.i[0] << endl;
}

static CAT_INLINE void ec_add_ref(const ecpt &a, const ecpt &b, ecpt &r, const bool z2_one, const bool in_precomp_t1, const bool out_precomp_t1, ufe &t2b) {
	ufe A,B,C,D,E,F,G,H;

	fe_mul(a.x, b.x, A);
	fe_mul(a.y, b.y, B);

	if (!in_precomp_t1) {
		fe_mul(a.t, t2b, C);
	} else {
		fe_set(a.t, C);
	}

	fe_mul(C, b.t, C);
	fe_mul_u(C, C);
	fe_mul_smallk(C, EC_D, C);

	if (z2_one) {
		fe_set(a.z, D);
	} else {
		fe_mul(a.z, b.z, D);
	}

	ufe E1, E2;
	fe_add(a.x, a.y, E1);
	fe_add(b.x, b.y, E2);
	fe_mul(E1, E2, E);
	fe_sub(E, A, E);
	fe_sub(E, B, E);

	fe_sub(D, C, F);
	fe_add(D, C, G);

	ufe H1;
	fe_mul_u(A, H1);
	fe_add(B, H1, H);

	fe_mul(E, F, r.x);
	fe_mul(G, H, r.y);

	if (out_precomp_t1) {
		fe_mul(E, H, r.t);
	} else {
		fe_set(H, t2b);
		fe_set(E, r.t);
	}

	fe_mul(F, G, r.z);
}

static CAT_INLINE void ec_dbl_ref(const ecpt &p, ecpt &r, const bool z_one, ufe &t2b) {
	ufe A, B, C, D, E, F, G, H;

	fe_sqr(p.x, A);
	fe_sqr(p.y, B);
	if (z_one) {
		fe_set_smallk(1, C);
	} else {
		fe_sqr(p.z, C);
	}
	fe_add(C, C, C);
	fe_neg(A, D);
	fe_mul_u(D, D);

	fe_add(p.x, p.y, E);
	fe_sqr(E, E);
	fe_sub(E, A, E);
	fe_sub(E, B, E);

	fe_add(D, B, G);
	fe_sub(G, C, F);
	fe_sub(D, B, H);

	fe_mul(E, F, r.x);
	fe_mul(G, H, r.y);
	fe_set(E, r.t);
	fe_set(H, t2b);
	fe_mul(F, G, r.z);
}

static bool ec_dbl_add_ref_test(const u64 *scalar, int mode) {
	ufe tp, tq;

	const bool z2_one = (mode & 1) != 0;
	const bool use_add = (mode & 2) != 0;
	const bool use_eg = (mode & 4) != 0;

	// Verify that dbl and add are the same as the simpler reference implementation

	bool seen = false;

	ecpt p, q;
	ec_identity(p);
	ec_identity(q);
	fe_zero(tp);
	fe_zero(tq);

	ecpt g;
	if (use_eg) {
		ec_set(EC_EG, g);
	} else {
		ec_set(EC_G, g);
	}

	for (int ii = 255; ii >= 0; --ii) {
		if (seen) {
			if (fe_iszero(p.x)) {
				cout << "Zero at dbl " << ii << endl;
				return false;
			}
			if (use_add) {
				// Verifies add is unified
				ec_add(p, p, p, false, false, false, tp);
			} else {
				ec_dbl(p, p, false, tp);
			}
			ec_dbl_ref(q, q, false, tq);

			ecpt_affine aa, ba;
			ecpt a, b;
			ec_affine(p, aa);
			ec_affine(q, ba);
			ec_expand(aa, a);
			ec_expand(ba, b);

			if (!ec_isequal(a, b)) {
				cout << "Double mismatch mode " << mode << endl;
				return false;
			}
		}

		if ((scalar[ii / 64] >> (ii % 64)) & 1) {
			if (seen && fe_iszero(p.x)) {
				cout << "Zero at add " << ii << endl;
				return false;
			}
			ec_add(p, g, p, z2_one, false, false, tp);
			ec_add_ref(q, g, q, z2_one, false, false, tq);
			seen = true;

			ecpt_affine aa, ba;
			ecpt a, b;
			ec_affine(p, aa);
			ec_affine(q, ba);
			ec_expand(aa, a);
			ec_expand(ba, b);

			if (!ec_isequal(a, b)) {
				cout << "Add mismatch mode " << mode << endl;
				return false;
			}
		}
	}

	return true;
}

static bool ec_curve_order_test(bool use_eg, bool z2_one) {
	ufe t2b;

	// Verify that curve has the specified order

	bool seen = false;

	ecpt p;
	ec_identity(p);
	fe_zero(t2b);

	ecpt g;
	if (use_eg) {
		ec_set(EC_EG, g);
	} else {
		ec_set(EC_G, g);
	}

	for (int ii = 255; ii >= 0; --ii) {
		if (seen) {
			if (fe_iszero(p.x)) {
				cout << "Zero at dbl " << ii << endl;
				return false;
			}
			ec_dbl(p, p, false, t2b);
		}

		if ((PRIME_ORDER[ii / 64] >> (ii % 64)) & 1) {
			if (seen && fe_iszero(p.x)) {
				cout << "Zero at add " << ii << endl;
				return false;
			}
			ec_add(p, g, p, z2_one, false, false, t2b);
			seen = true;
		}
	}

	ecpt_affine pa;
	ec_affine(p, pa);
	ec_expand(pa, p);
	//ec_print(p);

	ufe one;
	fe_set_smallk(1, one);

	// Verifies base points are in q torsion group

	if (!fe_iszero(p.x) || !fe_isequal(p.y, one)) {
		return false;
	}

	ec_dbl(p, p, false, t2b);
	ec_dbl(p, p, false, t2b);

	// Verifies that doubling works on identity element

	ec_affine(p, pa);
	ec_expand(pa, p);
	//ec_print(p);

	return fe_iszero(p.x) && fe_isequal(p.y, one);
}

static bool ec_base_point_test() {
	ecpt g, eg;
	ec_set(EC_G, g);
	ec_set(EC_EG, eg);

	// Verify that base point and its endomorphism are on the curve

	if (!ec_valid(g.x, g.y)) {
		cout << "Generator not on curve" << endl;
		return false;
	}

	if (!ec_valid(eg.x, eg.y)) {
		cout << "Generator endomorphism not on curve" << endl;
		return false;
	}

	// Verify that precomputed endomorphism for base point is correct

	ecpt eg2;
	gls_morph(EC_G.x, EC_G.y, eg2.x, eg2.y);

	if (!fe_isequal(eg2.x, EC_EG.x)) {
		return false;
	}

	if (!fe_isequal(eg2.y, EC_EG.y)) {
		return false;
	}

	// Verify that the T and Z coordinates for base point
	// and endomorphism of base point are all correct

	ufe r, one;
	fe_set_smallk(1, one);

	fe_mul(g.x, g.y, r);
	if (!fe_isequal(g.t, r)) {
		cout << "Generator T wrong" << endl;
		return false;
	}
	if (!fe_isequal(g.z, one)) {
		cout << "Generator Z wrong" << endl;
		return false;
	}

	fe_mul(eg.x, eg.y, r);
	if (!fe_isequal(eg.t, r)) {
		cout << "Generator endo T wrong" << endl;
		return false;
	}
	if (!fe_isequal(eg.z, one)) {
		cout << "Generator endo Z wrong" << endl;
		return false;
	}

	return true;
}

static bool ec_zero_test() {
	ecpt r;

	ec_set(EC_G, r);

	if (fe_iszero(r.x)) {
		return false;
	}
	if (fe_iszero(r.y)) {
		return false;
	}
	if (fe_iszero(r.t)) {
		return false;
	}
	if (fe_iszero(r.z)) {
		return false;
	}

	ec_zero(r);

	if (!fe_iszero(r.x)) {
		return false;
	}
	if (!fe_iszero(r.y)) {
		return false;
	}
	if (!fe_iszero(r.t)) {
		return false;
	}
	if (!fe_iszero(r.z)) {
		return false;
	}

	return true;
}

static bool ec_expand_test() {
	ecpt_affine s;
	fe_set(EC_GX, s.x);
	fe_set(EC_GY, s.y);

	ecpt t;
	ec_expand(s, t);

	if (!ec_isequal(EC_G, t)) {
		return false;
	}

	return true;
}

static bool ec_identity_test() {
	ecpt a, p;

	ec_identity(p);
	ec_set(p, a);

	ufe tp;
	ec_dbl(p, p, false, tp);

	ecpt_affine pa;
	ec_affine(p, pa);
	ec_expand(pa, p);

	if (!ec_isequal(a, p)) {
		return false;
	}

	ec_add(p, p, p, false, true, false, tp);

	ec_affine(p, pa);
	ec_expand(pa, p);

	if (!ec_isequal(a, p)) {
		return false;
	}

	return true;
}

static bool ec_save_load_test(const ecpt &P) {
	ecpt_affine a, r;

	fe_set(P.x, a.x);
	fe_set(P.y, a.y);

	u8 buffer[65] = {0};

	ec_save_xy(a, buffer);

	ec_load_xy(buffer, r);

	if (!ec_isequal_xy(r, a)) {
		return false;
	}

	if (buffer[64] != 0) {
		return false;
	}

	return true;
}

static bool ec_valid_test() {
	ecpt a;

	ec_set(EC_G, a);
	if (!ec_valid(a.x, a.y)) {
		return false;
	}

	ec_set(EC_EG, a);
	if (!ec_valid(a.x, a.y)) {
		return false;
	}

	ufe one;
	fe_set_smallk(1, one);

	if (ec_valid(Cnn, one)) {
		return false;
	}
	if (ec_valid(C0n, one)) {
		return false;
	}
	if (ec_valid(Cn0, one)) {
		return false;
	}
	if (ec_valid(C00, one)) {
		return false;
	}

	if (ec_valid(GXn, EC_GY)) {
		return false;
	}
	if (!ec_valid(EC_GX, EC_GY)) {
		return false;
	}

	return true;
}

static bool ec_neg_test(const ecpt &P) {
	ecpt a, d, n, p, i;
	ufe tp;

	ec_set(P, a);

	ec_identity(i);

	ec_neg(a, n);

	ec_add(a, n, p, false, true, false, tp);

	ecpt_affine pa;
	ec_affine(p, pa);
	ec_expand(pa, p);

	if (!ec_isequal(i, p)) {
		return false;
	}

	ecpt nd;
	ec_dbl(n, nd, false, tp);
	fe_mul(tp, nd.t, nd.t);

	ec_dbl(a, d, false, tp);

	ec_add(d, nd, p, false, false, false, tp);

	ec_affine(p, pa);
	ec_expand(pa, p);

	if (!ec_isequal(i, p)) {
		return false;
	}

	return true;
}

static bool ec_gen_mask_test(const ecpt &P) {
	int x, y;
	u128 mask;
	const u128 N = (s128)(-1);

	x = 0;
	y = 0;

	mask = ec_gen_mask(x, y);
	if (mask != N) {
		return false;
	}

	x = 1;
	y = 0;

	mask = ec_gen_mask(x, y);
	if (mask != 0) {
		return false;
	}

	x = 0x7fffffff;
	y = x;

	mask = ec_gen_mask(x, y);
	if (mask != N) {
		return false;
	}

	x = 0x7fffffff;
	y = 1;

	mask = ec_gen_mask(x, y);
	if (mask != 0) {
		return false;
	}

	return true;
}

static bool ec_xor_mask_test(const ecpt &P) {
	ecpt a, p;

	ec_set(P, a);

	ec_zero(p);

	u128 mask;
  
	ecpt z;
	ec_zero(z);

	mask = ec_gen_mask(1, 3);
	ec_xor_mask(a, mask, p);
	if (!ec_isequal(z, p)) {
		return false;
	}

	mask = ec_gen_mask(3, 3);
	ec_xor_mask(a, mask, p);
	if (!ec_isequal(P, p)) {
		return false;
	}

	return true;
}

static bool ec_set_mask_test() {
	ecpt p;

	ec_zero(p);

	u128 mask;
  
	ecpt z;
	ec_zero(z);

	mask = ec_gen_mask(1, 3);
	ec_set_mask(EC_G, mask, p);
	if (!ec_isequal(z, p)) {
		return false;
	}

	mask = ec_gen_mask(3, 3);
	ec_set_mask(EC_G, mask, p);
	if (!ec_isequal(EC_G, p)) {
		return false;
	}

	mask = ec_gen_mask(4, 4);
	ec_set_mask(EC_EG, mask, p);
	if (!ec_isequal(EC_EG, p)) {
		return false;
	}

	mask = ec_gen_mask(2, 4);
	ec_set_mask(EC_G, mask, p);
	if (!ec_isequal(EC_EG, p)) {
		return false;
	}

	return true;
}

static bool ec_neg_mask_test(const ecpt &P) {
	ecpt p, n;

	ec_set(P, p);
	ec_neg(p, n);

	u128 mask;
  
	mask = ec_gen_mask(1, 3);
	ec_neg_mask(mask, p);
	if (!ec_isequal(P, p)) {
		return false;
	}

	mask = ec_gen_mask(3, 3);
	ec_neg_mask(mask, p);
	if (!ec_isequal(n, p)) {
		return false;
	}

	mask = ec_gen_mask(4, 4);
	ec_neg_mask(mask, p);
	if (!ec_isequal(P, p)) {
		return false;
	}

	mask = ec_gen_mask(3, 10);
	ec_neg_mask(mask, p);
	if (!ec_isequal(P, p)) {
		return false;
	}

	return true;
}

static bool ec_cond_neg_test(const ecpt &P) {
	ecpt p, n;

	ec_set(P, p);
	ec_neg(p, n);

	ec_cond_neg(0, p);
	if (!ec_isequal(P, p)) {
		return false;
	}

	ec_cond_neg(1, p);
	if (!ec_isequal(n, p)) {
		return false;
	}

	ec_cond_neg(1, p);
	if (!ec_isequal(P, p)) {
		return false;
	}

	ec_cond_neg(0, p);
	if (!ec_isequal(P, p)) {
		return false;
	}

	return true;
}

static bool ec_cond_add_test(const ecpt &P) {
	ecpt a, b, c, p;

	ec_set(P, p);

	ufe tp, tb;
	ec_add(p, p, p, false, true, false, tp);

	ecpt_affine pa;
	ec_affine(p, pa);
	ec_expand(pa, p);

	ec_set(P, a);
	ec_set(P, b);

	ec_cond_add(0, b, a, b, false, true, tb);

	ecpt_affine ba;
	ec_affine(b, ba);
	ec_expand(ba, c);

	if (!ec_isequal(P, c)) {
		return false;
	}

	ec_cond_add(1, b, a, b, false, false, tb);

	ec_affine(b, ba);
	ec_expand(ba, c);

	if (!ec_isequal(p, c)) {
		return false;
	}

	return true;
}


//// Entrypoint

int main() {
	cout << "Snowshoe Unit Tester: Elliptic Curve Point Operations" << endl;

	// ec_valid, gls_morph, base point validity:
	assert(ec_base_point_test());

	// ec_add, ec_dbl <-> ec_add_ref, ec_dbl_ref:
	for (int ii = 0; ii < 8; ++ii) {
		assert(ec_dbl_add_ref_test(RS1, ii));
		assert(ec_dbl_add_ref_test(RS2, ii));
		assert(ec_dbl_add_ref_test(RS3, ii));
		assert(ec_dbl_add_ref_test(RS4, ii));
		assert(ec_dbl_add_ref_test(RS5, ii));
		assert(ec_dbl_add_ref_test(RS6, ii));
		assert(ec_dbl_add_ref_test(RS7, ii));
		assert(ec_dbl_add_ref_test(RS8, ii));
	}

	// ec_add, ec_dbl, ec_identity, ec_affine, curve order validity:
	assert(ec_curve_order_test(false, true));
	assert(ec_curve_order_test(false, false));
	assert(ec_curve_order_test(true, true));
	assert(ec_curve_order_test(true, false));

	// ec_zero:
	assert(ec_zero_test());

	// ec_expand:
	assert(ec_expand_test());

	// ec_save_xy, ec_load_xy:
	assert(ec_save_load_test(EC_G));
	assert(ec_save_load_test(EC_EG));

	// ec_add, ec_dbl, ec_identity:
	assert(ec_identity_test());

	// ec_valid:
	assert(ec_valid_test());

	// ec_neg:
	assert(ec_neg_test(EC_G));
	assert(ec_neg_test(EC_EG));

	// ec_xor_mask:
	assert(ec_xor_mask_test(EC_G));
	assert(ec_xor_mask_test(EC_EG));

	// ec_set_mask:
	assert(ec_set_mask_test());

	// ec_neg_mask:
	assert(ec_neg_mask_test(EC_G));
	assert(ec_neg_mask_test(EC_EG));

	// ec_gen_mask:
	assert(ec_gen_mask_test(EC_G));
	assert(ec_gen_mask_test(EC_EG));

	// ec_cond_neg:
	assert(ec_cond_neg_test(EC_G));
	assert(ec_cond_neg_test(EC_EG));

	// ec_cond_add:
	assert(ec_cond_add_test(EC_G));
	assert(ec_cond_add_test(EC_EG));

	cout << "All tests passed successfully." << endl;

	return 0;
}


