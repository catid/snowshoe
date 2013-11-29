// Elliptic curve point operations
#include "ecpt.cpp"
#include "misc.cpp"

#include <iostream>
#include <iomanip>
using namespace std;

/*
 * Mask a random number to produce a compatible scalar for multiplication
 */

void ec_mask_scalar(u64 k[4]) {
	// Prime order of the curve = q, word-mapped:
	// 0x0FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA6261414C0DC87D3CE9B68E3B09E01A5
	//   (      3       )(       2      )(       1      )(       0      )

	// Clear high 5 bits
	// Clears one extra bit to simplify key generation
	k[3] &= 0x07FFFFFFFFFFFFFFULL;

	// Largest value after filtering:
	// 0x07FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
	//   (      3       )(       2      )(       1      )(       0      )
}

/*
 * GLV-SAC Scalar Recoding Algorithm for m=2 [1]
 *
 * Returns low bit of 'a'.
 */

static u32 ec_recode_scalars_2(ufp &a, ufp &b, const int len) {
	u32 lsb = ((u32)a.i[0] & 1) ^ 1;
	a.w -= lsb;
	a.w >>= 1;
	a.w |= (u128)1 << (len - 1);

	const u128 an = ~a.w;
	u128 mask = 1;
	for (int ii = 1; ii < len; ++ii) {
		const u128 anmask = an & mask;
		b.w += (b.w & anmask) << 1;
		mask <<= 1;
	}

	return lsb;
}

/*
 * GLV-SAC Scalar Recoding Algorithm for m=4 [1]
 *
 * Returns low bit of 'a'.
 */

static u32 ec_recode_scalars_4(ufp &a, ufp &b, ufp &c, ufp &d, const int len) {
	u32 lsb = ((u32)a.i[0] & 1) ^ 1;
	a.w -= lsb;
	a.w >>= 1;
	a.w |= (u128)1 << (len - 1);

	const u128 an = ~a.w;
	u128 mask = 1;
	for (int ii = 1; ii < len; ++ii) {
		const u128 anmask = an & mask;
		b.w += (b.w & anmask) << 1;
		c.w += (c.w & anmask) << 1;
		d.w += (d.w & anmask) << 1;
		mask <<= 1;
	}

	return lsb;
}


//// Constant-time Point Multiplication

/*
 * Precomputed table generation
 *
 * Using GLV-SAC Precomputation with m=2 [1], assuming window size of 2 bits
 *
 * Window of 2 bits table selection:
 *
 * aa bb -> evaluated (unsigned table index), sign
 * 00 00    -3a + 0b (0)-
 * 00 01    -3a - 1b (1)-
 * 00 10    -3a - 2b (2)-
 * 00 11    -3a - 3b (3)-
 * 01 00    -1a + 0b (4)-
 * 01 01    -1a + 1b (5)-
 * 01 10    -1a - 2b (6)-
 * 01 11    -1a - 1b (7)-
 * 10 00    1a + 0b (4)+
 * 10 01    1a - 1b (5)+
 * 10 10    1a + 2b (6)+
 * 10 11    1a + 1b (7)+
 * 11 00    3a + 0b (0)+
 * 11 01    3a + 1b (1)+
 * 11 10    3a + 2b (2)+
 * 11 11    3a + 3b (3)+
 *
 * Table index is simply = (a0 ^ a1) || b1 || b0
 *
 * The differences above from [1] seem to improve the efficiency of evaulation
 * and they make the code easier to analyze.
 */

static void ec_gen_table_2(const ecpt &a, const ecpt &b, ecpt TABLE[8]) {
	ecpt bn;
	ec_neg(b, bn);

	// P[4] = a
	ec_set(a, TABLE[4]);

	// P[5] = a - b
	ufe t2b;
	ec_add(a, bn, TABLE[5], true, true, true, t2b);

	// P[7] = a + b
	ec_add(a, b, TABLE[7], true, true, true, t2b);

	// P[6] = a + 2b
	ec_add(TABLE[7], b, TABLE[6], true, true, true, t2b);

	ecpt a2;
	ec_dbl(a, a2, true, t2b);

	// P[0] = 3a
	ec_add(a2, a, TABLE[0], true, false, true, t2b);

	// P[1] = 3a + b
	ec_add(TABLE[0], b, TABLE[1], true, true, true, t2b);

	// P[2] = 3a + 2b
	ec_add(TABLE[1], b, TABLE[2], true, true, true, t2b);

	// P[3] = 3a + 3b
	ec_add(TABLE[2], b, TABLE[3], true, true, true, t2b);
}

/*
 * Table index is simply = (a0 ^ a1) || b1 || b0
 */

static CAT_INLINE void ec_table_select_2(const ecpt *table, const ufp &a, const ufp &b, const int index, ecpt &r) {
	u32 bits = (u32)(a.w >> index);
	u32 k = ((bits ^ (bits >> 1)) & 1) << 2;
	k |= (u32)(b.w >> index) & 3;

	ec_zero(r);

	for (int ii = 0; ii < 8; ++ii) {
		// Generate a mask that is -1 if ii == index, else 0
		const u128 mask = ec_gen_mask(ii, k);

		// Add in the masked table entry
		ec_xor_mask(table[ii], mask, r);
	}

	ec_cond_neg(((bits >> 1) & 1) ^ 1, r);
}

/*
 * Multiplication by variable base point
 *
 * Preconditions:
 * 	0 < k < q
 *
 * Multiplies the point by k * 4 and stores the result in R
 */

// R = k*4*P
void ec_mul(const u64 k[4], const ecpt_affine &P0, ecpt_affine &R) {
	// Decompose scalar into subscalars
	ufp a, b;
	s32 asign, bsign;
	gls_decompose(k, asign, a, bsign, b);

	// Q = endomorphism of P
	ecpt_affine Qa;
	gls_morph(P0.x, P0.y, Qa.x, Qa.y);
	ecpt Q;
	ec_expand(Qa, Q);
	ec_cond_neg(bsign, Q);

	// Set base point signs
	ecpt P;
	ec_expand(P0, P);
	ec_cond_neg(asign, P);

	// Precompute multiplication table
	ecpt table[8];
	ec_gen_table_2(P, Q, table);

	// Recode subscalars
	u32 recode_bit = ec_recode_scalars_2(a, b, 128);

	// Initialize working point
	ecpt X;
	ec_table_select_2(table, a, b, 126, X);

	ufe t2b;
	for (int ii = 124; ii >= 0; ii -= 2) {
		ecpt T;
		ec_table_select_2(table, a, b, ii, T);

		ec_dbl(X, X, false, t2b);
		ec_dbl(X, X, false, t2b);
		ec_add(X, T, X, false, false, false, t2b);
	}

	// If bit == 1, X <- X + P (inverted logic from [1])
	ec_cond_add(recode_bit, X, P, X, true, false, t2b);

	// Multiply by 4 to avoid small subgroup attack
	ec_dbl(X, X, false, t2b);
	ec_dbl(X, X, false, t2b);

	// Compute affine coordinates in R
	ec_affine(X, R);
}


//// Constant-time Generator Base Point Multiplication

/*
 * Multiplication by generator base point
 *
 * Using modified LSB-Set Comb Method from [1].
 *
 * The algorithm is tuned with ECADD = 1.64 * ECDBL in cycles.
 *
 * t = 252 bits for input scalars
 * w = window size in bits
 * v = number of tables
 * e = roundup(t / wv)
 * d = e * v
 * l = d * w
 *
 * The parameters w,v are tunable.  The number of table entries:
 *	= v * 2^(w - 1)
 *
 * Number of ECADDs and ECDBLs = e - 1, d - 1, respectively.
 *
 * Constant-time table lookups are expensive and proportional to
 * the number of overall entries.  After some experimentation,
 * 256 is too many table entries.  Halving it to 128 entries is
 * a good trade-off.
 *
 * Optimizing for operation counts, choosing values of v,w that
 * yield tables of 128 entries:
 *
 * v,w -> e,d   -> effective ECDBLs
 * 1,8 -> 32,32 -> 81.84
 * 2,7 -> 18,36 -> 74.4 <- best option
 * 4,6 -> 11,44 -> 80.52
 */

static const u64 PRECOMP_TABLE_0[] = {
0xfULL, 0x0ULL, 0x0ULL, 0x0ULL,
0x36d073dade2014abULL, 0x7869c919dd649b4cULL, 0xdd9869fe923191b0ULL, 0x6e848b46758ba443ULL,
0xc0257189412dee27ULL, 0x22d1b2a099cef701ULL, 0x467a15261c3e929dULL, 0x7fede0e4cf68d988ULL,
0x80c3dc5f34ad2f0cULL, 0x6e4c44e71fab5f84ULL, 0x9cae3727bb435cbcULL, 0x267325c8944698f8ULL,
0xa9f3f1342d5833faULL, 0xd713d9ca10dbf27ULL, 0x22c52394537fef93ULL, 0x30e11fa9329422aeULL,
0x3396304477e71d78ULL, 0x239b72e696b1d33eULL, 0x91fd62721ceb91e4ULL, 0x57f6acba2654f846ULL,
0xb347fe4dd0630434ULL, 0x1ee77493a307590dULL, 0x727d670b78421fe2ULL, 0x44c3c273df251de0ULL,
0x71dbad492800594ULL, 0x1f410b55ab343b26ULL, 0xc7a2de19aca789e4ULL, 0x32249d176df691b7ULL,
0x702fde8b105e4ce4ULL, 0x2f2baec7f8ee114eULL, 0x760e745252b16b6aULL, 0x3a8e355a529d5777ULL,
0xf0445c7cac272c35ULL, 0x66c1b4456dcb384dULL, 0x93f3e6fcb4ff83dULL, 0x1d7e4261a44beae7ULL,
0xfcc9d435a689819ULL, 0x32b55195f81677fULL, 0x21a5c3fd80210bafULL, 0x43e6cd7ebde1a73cULL,
0x2664512b4f034b84ULL, 0x224555423e27897bULL, 0xf66b3e82fa8172caULL, 0x418c43ba2b2c2cdfULL,
0xc8e0a154289d7217ULL, 0x3ccb4f21c7535eaULL, 0xad62d02a33d27dULL, 0x64f28d58eb112c20ULL,
0xc641aabb08f81e8cULL, 0x19b6ac02e86a4a74ULL, 0x3240fe9b5cfbb25fULL, 0x1c5ae3311fe73f52ULL,
0x1a2d06c48ec87eb7ULL, 0x595f030c63cb9f75ULL, 0x7e4ce069d0252eaaULL, 0x2b67e2850665113fULL,
0x53797417644d38b1ULL, 0x5c6cc0fad0961f35ULL, 0x1c47ce3bd26bfb69ULL, 0x509ac971615a6490ULL,
0x295d3d0793cf7f92ULL, 0xbabf7b7af2c6f03ULL, 0xd9377d2d2348d740ULL, 0x142f89c3f50f2c78ULL,
0x4464af962886f51fULL, 0x5a1be0fb7b7c812bULL, 0x8b1ab6e5cd7f2580ULL, 0x7ec25fb8091eb2d7ULL,
0x71fdf4c1d2ba8567ULL, 0xd2419ec9fbb0a5eULL, 0xa35fbd6da89f2d58ULL, 0x45d65de885459b51ULL,
0x6778ad3327edb696ULL, 0x643d770deacc423cULL, 0xa6063f6c5f992e06ULL, 0xde1e55a4477e352ULL,
0xd3df4250a4383b1bULL, 0x3e653721ff649bcbULL, 0x27cbfed7de7c7680ULL, 0x7508840bbc9a54a8ULL,
0x43e0480759ef1863ULL, 0x5e686fff9f8e00e9ULL, 0x62cc9dd5496d699eULL, 0x3dcbb84cc5aaf8f5ULL,
0x36e276fd00d1ca18ULL, 0x733097fdcac089b4ULL, 0x39a3d9e200c064e1ULL, 0x22c9e1ea50815782ULL,
0xa3f0cf3ef55059cULL, 0xb1213397520fde9ULL, 0x37cc47b2242bb840ULL, 0x3c56d498f324dabaULL,
0x911c26d3a8883f7eULL, 0x16bfdb623a597dc2ULL, 0xc50f51345a072322ULL, 0x3417853e70d94a89ULL,
0x3bc3895442e26a41ULL, 0x33b83a40d4bb7db1ULL, 0x517160b5a3df03b8ULL, 0x42a917d97b4053ecULL,
0x6210bdefa9e60942ULL, 0x18a3f7cd5d63a070ULL, 0xf74b507ebce9e116ULL, 0x655d223c7c63d29eULL,
0xc151db63dc28a09eULL, 0x678a968a262d4f9bULL, 0x8f631a0b2e7f7e2ULL, 0x756c53c36ffa9a22ULL,
0x5806b239a7a25c91ULL, 0x2a3598349b7e9445ULL, 0xcff2a7dafa5261ffULL, 0x2b3b19360ff8ed7ULL,
0xc5a7c2c977e23b2cULL, 0x663ce2499666bb0bULL, 0x818c546442ce6eaeULL, 0x67e824976d91a0a9ULL,
0x4ae3d06ce1eb701aULL, 0x2854f7079b7d5748ULL, 0x802442310d72f39cULL, 0x3907695186210ce2ULL,
0x4ec6a01caa4ae626ULL, 0x72974ed5a942147aULL, 0x9a0018c3383f57afULL, 0x60cd656b6ba43c86ULL,
0x5f106965ad6f60a0ULL, 0x1af09f2d0d29b7cfULL, 0x4eeeb180862ee54bULL, 0x47386eeef34f656aULL,
0x596a47a1bd9420a0ULL, 0x5611e38366cc1796ULL, 0x8df59d3e0fb621a6ULL, 0x5c2c1a8adc8cefe4ULL,
0xac1b40c5812e7e50ULL, 0x2da0c7d73d2bc48dULL, 0x8196f6c795a08fdULL, 0x4c35cdddb6dbe1d3ULL,
0x1307f06bd01dbdf8ULL, 0x2dd307afe4eb1180ULL, 0x1cc1a30f8d09f72cULL, 0x18c6880d3e28ccbeULL,
0x364615ea188011bdULL, 0x7aeaf3cc1abb3e3eULL, 0x6394782e5bc5622eULL, 0x141ca3afb157e79bULL,
0x26b300fa6f9694c9ULL, 0x164c54f0a51f3461ULL, 0x25ff971129485daeULL, 0x93318c3594df963ULL,
0xe5b0becdd878685aULL, 0xa7e933ae9ce3901ULL, 0xfbd70a101b589de0ULL, 0x6c360d48ead2a288ULL,
0x7d1a934336045c3aULL, 0x9971e128d3171bbULL, 0x897d16e25e23642bULL, 0x73c833238aeb43ULL,
0x246157b0c59feebcULL, 0x23930effceea42c2ULL, 0x1bc25fa6a2f2f87eULL, 0x43947bd545272a7ULL,
0x6c8164b134ce465aULL, 0x56289db00b3b0c66ULL, 0xd0b82225bbca8cc2ULL, 0x1062e75b29735fdULL,
0x486e795e1d05003eULL, 0x34f26d8bc540d798ULL, 0xe7d4505bd374e859ULL, 0x4cb88344c6aa3fc9ULL,
0xd006937fe8bb1e38ULL, 0x3707c25a597409deULL, 0x54d647a4e6cbd93aULL, 0x29f7ed8512e88d6eULL,
0xa2cba8a959900fd7ULL, 0x40eff4737e8f8c5aULL, 0xfd2fe210d5834093ULL, 0x525ff5c59d47781cULL,
0x6902367b210b2b7dULL, 0x346607a6772cb7a2ULL, 0x63a67723039f88d8ULL, 0x5003b4c2c8aec14cULL,
0x8c3096cd19cee08aULL, 0x56a2eee216afee7cULL, 0x3216c87f26d4ef8ULL, 0x15cf73352d1dd930ULL,
0x2091167afcec33b1ULL, 0x3232e1e814579132ULL, 0x3d894716a60d8c9dULL, 0x1dbd947433a7a0e7ULL,
0x4b962708f50de47ULL, 0x3f263060999561b5ULL, 0x445dfa684afb8065ULL, 0x7e6a41f05025964fULL,
0xcc515fdc2a996772ULL, 0x1688acc005473d04ULL, 0x40fed523adfbcaa6ULL, 0x4dd6f0927c227f04ULL,
0x6bee658a504ee8eaULL, 0x513913762537c16fULL, 0x53ea13bf8ba872d3ULL, 0x68828ca3919a38feULL,
0x6ab1c4ee45ea88c1ULL, 0x4896ace7baa4cebeULL, 0xb4b7f7a2a881e8dbULL, 0x7ce43c96aab0dda1ULL,
0x3fec69619bea6715ULL, 0x24fa602b2557cc66ULL, 0x23d372f8d0ebb104ULL, 0x457c1c459b303607ULL,
0xd3abbf07b7e02e58ULL, 0x68272102076b41ecULL, 0x2bbdf1128f5664d0ULL, 0x1d785e2044f12961ULL,
0xb09afc08021f628dULL, 0x2e7e413d1a1fcef8ULL, 0x7f7c56c95a973b7aULL, 0xd71e67e44450257ULL,
0xba7327cb15529951ULL, 0x233a365f8c8672acULL, 0xd632b507b8b60244ULL, 0xf0e3cdf4a2be9d7ULL,
0xd7f49c0568a14ddULL, 0x77d025da314730c4ULL, 0x33e0211742a1e9e7ULL, 0x5e0bd8f665266c0ULL,
0x9aebd2f4e21c94d7ULL, 0x6109a0940c5eb2cdULL, 0x9898f59fdd049337ULL, 0x7a92437a9d34e47fULL,
0x1a3f3fc106fab9aaULL, 0x1162b24aea2bb259ULL, 0xc86a152f06d72a17ULL, 0x23ff651ea14fdae0ULL,
0x653345ff325fccbdULL, 0x5281ac6e0f951094ULL, 0x828345c207be1aa0ULL, 0x61f1e49006f0931fULL,
0x3441260df3f8b144ULL, 0x548d223dbe1d9e0fULL, 0x72caf061bd8627bdULL, 0x47171b3562ee4904ULL,
0x17316a6b40d4b41fULL, 0x6f9f888cdffc867bULL, 0x98701eeaebfb1fe2ULL, 0x31dc9ce075edeb8ULL,
0x930e9a214b9702c5ULL, 0x5020f7c88a9824a0ULL, 0xb3ba00b6a8adafb1ULL, 0x3c9f65e1cb75faf6ULL,
0xe44d2f3170ef5656ULL, 0x7c299a1469a9914dULL, 0xc789162aa130c444ULL, 0x17231bc28109a8a7ULL,
0x30ef872e95d24514ULL, 0xdfed31c326eae5dULL, 0x44bb63cdfb557ba6ULL, 0x4eb5a67b005cf9bULL,
0xf22a5d5ac0caa911ULL, 0x3af1e63ab024c9adULL, 0x6e3dc0129da1c801ULL, 0x5baf6ff456cb9eaaULL,
0x4a876845e49569e7ULL, 0x5ff482757569fd03ULL, 0x97f81cc920a15d6eULL, 0x2129cdd20dbd7c16ULL,
0x2dcf87aa21fbe49aULL, 0x6380875efd5f8109ULL, 0x91cbf5dfd289baebULL, 0x7023b5752169a9f3ULL,
0x437b07ce27e83414ULL, 0x564477b46dcd4b58ULL, 0x733501b4e86df5edULL, 0x2d10b6e6ac4fdb5eULL,
0x30f72802d52080f1ULL, 0x31f4c1ae0e5de530ULL, 0xc5cc440ce3ab3f8dULL, 0x24b2985ee52e1db3ULL,
0xe3ac866d05b6dbdeULL, 0x9067409d858f0c7ULL, 0xa3cdcd75f3185a05ULL, 0x4c7e53bcb7a628eeULL,
0xb5ace6dedb4a7aa0ULL, 0x72bf61550b543a0bULL, 0x6b8ab274e3819481ULL, 0x3270a9b473a4b39aULL,
0xfe1d27507a349721ULL, 0x2982e0f0d3a8b8b4ULL, 0x1241af4ae9097d39ULL, 0x3fd5a48fb7587d9dULL,
0xf92698a6783441e4ULL, 0x5702658d0b633edbULL, 0x775b17930625c4dbULL, 0x3d50c2741ad6994cULL,
0xf62e147e3ddc8839ULL, 0x6f5ee2ee4fbabba3ULL, 0xf653ea49e73b6ddbULL, 0xd7c156b3b8e4632ULL,
0xab9968641bede913ULL, 0x65612ae285605a55ULL, 0xb55f409d704d4b0fULL, 0xedc676e75904abfULL,
0x2642c2ef9bb133bcULL, 0x3ff7e35a91fed47bULL, 0xc538cc397073664eULL, 0x4b5fa53465648b9dULL,
0x49bb424e69be2b71ULL, 0x1b40d6a429c3551eULL, 0xef221689409f26d1ULL, 0x550e086b5faa72d0ULL,
0xe288e24060619719ULL, 0x535e665982eea560ULL, 0x3d29ea146dee8c76ULL, 0x2bac0b4d685553b2ULL,
0x4689025977b93c68ULL, 0x7ffb5f13237fb487ULL, 0xf110f755b5d7e60fULL, 0x75a9fa382f80e1c3ULL,
0xdd24fecc0003fb77ULL, 0xe2ab2b459162987ULL, 0xe1f963c251317e8bULL, 0x472c71a54a0f1253ULL,
0x89cd7034b6affe91ULL, 0x67e746104adfcce5ULL, 0xc2478f6a41abaf49ULL, 0x5aece0c950f52ef0ULL,
0xde902f7206c44e1dULL, 0x4cf2372eb9554572ULL, 0x479454efe51c176dULL, 0x791a45d552b9aef8ULL,
0x95be244ebd66a248ULL, 0x2a266a5026ff8222ULL, 0xe21babaa64e4bb27ULL, 0x1b3ae0acb6708d8aULL,
0x82db122e4267d7c8ULL, 0x64b807e95e432be1ULL, 0xa72698ad20d12093ULL, 0x1956b8809ebf8e7eULL,
0x81ca7e5bfd4518bcULL, 0x52d3e0c1f5c74e06ULL, 0x4fdd44ad83241d94ULL, 0x415e5370dee0de68ULL,
0x13cee309d6c6c382ULL, 0x7416758ee0ce6007ULL, 0x7fee9ccb99be2937ULL, 0x2141ed2927414c7cULL,
0xf2c17817d46767bbULL, 0x5dbfd5c348373e00ULL, 0x90dea6fdcb52f8cbULL, 0x768f7a9a80b829d6ULL,
0xd0c41d3cc65cbdd8ULL, 0x6aa2d5f7b8d6f2d7ULL, 0x170cd847ae7b9d96ULL, 0x38152096e221a4b9ULL,
0x8f0af4d84d265b5bULL, 0x3351f6efbde387cULL, 0x945253d626cf92f2ULL, 0x54b5e94e4c61e7d6ULL,
0xd3422d930825d1eaULL, 0x497693cb9fa948e8ULL, 0xbb05949623ade10cULL, 0x394a41da3c0b7486ULL,
0xe2c024091da70f74ULL, 0x49177244e286a068ULL, 0xd7c9cf6a67342420ULL, 0x69a9559ee5c9fadaULL,
0x95dc8c05bb291750ULL, 0x76a68ea95f5d0a12ULL, 0x5cde0ad22eb36106ULL, 0x59d5e909e376ace8ULL,
0xd731dc9810733bcfULL, 0x1bea3e18c7d34d51ULL, 0x5e9a6a6f5f0810d0ULL, 0x381bb194ae27fc62ULL,
0x36f36daab1b6b4f0ULL, 0x19c7af84e22ad8b3ULL, 0xd108569802551089ULL, 0x7a91139ed69126e6ULL,
0x7bfe7e764fcacf33ULL, 0x12752cec29ae1911ULL, 0x54c62c48fa6d61d6ULL, 0x1284440b7251a832ULL,
0x751ce0f9715bd48aULL, 0x698a5f1bb64b3577ULL, 0x9c3c0d016fa0863aULL, 0x345d6e089e97335bULL,
0x87ccfa8cb18a5a6ULL, 0x774a7e21bc215eebULL, 0x90b6a8646c442017ULL, 0x31d57460cfc3941eULL,
0x2dd191a58de3b4baULL, 0x5345b68da23f81a2ULL, 0xf075fc49eacde157ULL, 0x5314b7c72c0a9a4ULL,
0xe5fe1a78d70d2dcfULL, 0x381d770736bc8c76ULL, 0xc45ec4ce181a0c1eULL, 0x719372b475a88276ULL,
0x32669d1b43955991ULL, 0x507d12301b156163ULL, 0x3ee473094a9388aeULL, 0x1489ffddf00e5accULL,
0xeacdf08b0abcd47ULL, 0x24f9e5c86770729bULL, 0xc454cec42a855819ULL, 0x1e7e64a61f5be321ULL,
0x26d918c1e14f8bdbULL, 0x501c21920449ce1cULL, 0x929fbe253f7b47bULL, 0x5340a855a1b28ff0ULL,
0x32c2e8ffe93ab96eULL, 0x28cdfcd42fb83c9aULL, 0xba7b1700fc6fe4afULL, 0x79cf20a01f81fad3ULL,
0xe55c8ce629e337f1ULL, 0x16d504f0cf1d39d3ULL, 0x8010c256bf096934ULL, 0x3e0476f55818626dULL,
0xe49a15bc222e4db7ULL, 0x7f4ab95d62e98495ULL, 0xcfba4ea5fd040ce1ULL, 0x54c78e507221cd07ULL,
0x7bd6617b34bd085aULL, 0x53cd5e2334690530ULL, 0x2f617d73b17c1f54ULL, 0x306d3238a20027b2ULL,
0xbd0c4079305199fcULL, 0xe5c92724d5703c2ULL, 0x2dfe353d1f52d4fbULL, 0x43300b1ac6159594ULL,
0x451a5fa7a52f6babULL, 0x4ab437ad7436a97fULL, 0xcb53118447289394ULL, 0x476d4525367d3569ULL,
0x3123393aa3c3f3daULL, 0x7ea9fc46e48cbc6bULL, 0xcac437e951250201ULL, 0x4a20f85babe97435ULL,
0x1848eb2ffd8991daULL, 0x739f89a48326f3efULL, 0x1b18bf4816a4175aULL, 0x4f58e466f6f071faULL,
0x7c02bb5166f7eed9ULL, 0x1aee39c9c6aa33cfULL, 0xe11b0f369b83ee69ULL, 0x2e8dd73ca132525cULL,
0x689c5fd82d5e4eULL, 0xed67cffbec9eb55ULL, 0xe2cfde41564d3cf9ULL, 0x160c66753834549fULL,
0x18f930d812226a38ULL, 0x35d3795b2b557bdbULL, 0xac4bfb8b3184108fULL, 0x157c9971868d8149ULL,
0x1bd2312a229996c9ULL, 0x744c1d176074650dULL, 0x4db58f0f76d8ca31ULL, 0x60f9edc54ec9c7b8ULL,
0x8cfdd8d6009a9b7aULL, 0x32450fbb8a180bd7ULL, 0xd2fc1508a2cb99a4ULL, 0x9a2e1c06b660f30ULL,
0xfd9d641912fca23ULL, 0x1dc31e45241b661dULL, 0x6a2d8ac2a7611674ULL, 0x1c3bb9d4a1f988edULL,
0x975181cd1d995b51ULL, 0x23cf4a66dcdafff9ULL, 0x1d362b58213e65f0ULL, 0x1e8d4ce3e669c524ULL,
0xe7edba73c655b444ULL, 0x5062192274ad6166ULL, 0xb932a0ab5afd6172ULL, 0x6cc1e7fabf0fbebaULL,
0x8e10f6541914e4a5ULL, 0x2a7ecde3da40e6fcULL, 0x97765c1f4c72cb9bULL, 0xe87b21364d1a6fdULL,
0x86ce9280fe978f8fULL, 0x70333a1d4f0a4356ULL, 0xe2d52f4ab352c207ULL, 0x71bd9ad14b825641ULL,
0xc07cd01ba84c0dbcULL, 0x5b985c2efb6220eaULL, 0xa33f05d9a466b62fULL, 0x38632d9a068a5f1aULL,
0xf34f1b39fc9df81fULL, 0x370dfab03212a848ULL, 0x69507ffce5d22a6aULL, 0xc4a9be92549827cULL,
0x8a9c0a8edc315ee3ULL, 0x711000b92783d8feULL, 0x6cc798d4abd96cafULL, 0x7c0859be24e8df18ULL,
0x75b300f8eebaec75ULL, 0x28ea828c8680d217ULL, 0xec015bcfd5baa936ULL, 0x3c1e96f916201885ULL,
0xebe0cf4dbd6c4569ULL, 0x2d1040a866de1b72ULL, 0x8e3c4d6adfbeaa9cULL, 0xdd5a2702ccb8474ULL,
0xae327cc24e1283ffULL, 0x5279f3e3b4640a06ULL, 0x2ba26dc97f8f5a5fULL, 0x7dde8425c91b480bULL,
0x157b4fc963038e16ULL, 0x42906db00ec2ae4dULL, 0x18c9d8660e5d4b94ULL, 0x4e908f918f5a45e2ULL
};
static const u64 PRECOMP_TABLE_1[] = {
0xfad2788226688f06ULL, 0x1924dedc82fbf455ULL, 0xd329997d1b3aa551ULL, 0x497c5aae63352a0cULL,
0xc8b7075a25851bbcULL, 0x5fde998a13127f2cULL, 0x4e2444a52f60e3ULL, 0x7b088f6cd796e5cdULL,
0x577bcef66fc2302dULL, 0x386b9be656afa536ULL, 0x32f5cb220f4aef1fULL, 0x7bfd92f87d54be21ULL,
0xa14e9570f4d1116eULL, 0x86bd59cfe7c08f4ULL, 0xd2e8243f6771e8beULL, 0xe3a3804a3cae6f1ULL,
0xf2f5195dbc2c3363ULL, 0x79f3637b18cfa75cULL, 0xaebb2d3ada9053b6ULL, 0x6d445340190d6947ULL,
0x8030846b78811602ULL, 0x5e29614ae17e7ad3ULL, 0x430e5e16d5071a8cULL, 0x1038bdadb1284770ULL,
0x311715d47e2dd1ffULL, 0x765e4d432a368674ULL, 0x4e9f15bb5796d8c1ULL, 0x2fead6fc4881d3fcULL,
0x7ae4008bddbc8bcULL, 0x4ebe4ad4751fa048ULL, 0xe63e7b502293abf8ULL, 0x2e8e8da1c0165c1cULL,
0xab705e354de0b139ULL, 0x46aea02f39726cc5ULL, 0xd8ed93568e0c90faULL, 0x28fdb036ff3d119aULL,
0xfaffbbb4893207c4ULL, 0x4fe7638622edcd57ULL, 0x515b36ecabb359c6ULL, 0x706a749567f5120cULL,
0x8729773630df17a6ULL, 0x262f1c0fcf126b8aULL, 0x36291313245249fdULL, 0xfd41bbd375ca912ULL,
0xa20ef6e8a9bec78cULL, 0x4c1eb58f745bd6c9ULL, 0xbca4a40173b7f72bULL, 0x380a19ef263ebfa0ULL,
0xa082b1ac77f4777eULL, 0x4f7a0a22a8470bcbULL, 0x416af9704569a77eULL, 0x3573e678efe32cc4ULL,
0x645bb1fdc33bb3bcULL, 0x2d7239d958ee13f4ULL, 0x27458df841e54646ULL, 0x7707428748bc45a7ULL,
0x85a57caf7f144e65ULL, 0x64617c55f0f3c2a7ULL, 0x7884384ec82d51b8ULL, 0x520fa64a49f43931ULL,
0xe7dcbe03283827b5ULL, 0x207d2b6a2c4ec0a0ULL, 0x17dadd10d03908bfULL, 0x7b85d1a7f18fe56cULL,
0x6212e0b015f00adaULL, 0x83cdc88ce6da191ULL, 0xf6b96d2218942406ULL, 0x1d814330e0e4fa5fULL,
0x628cf68a803f8f25ULL, 0x6488f9e309311e63ULL, 0xf3552797f35afaa8ULL, 0x7557c5f29e789eb2ULL,
0x94909c14c24961b5ULL, 0x1ecb982e76a7b180ULL, 0xcec16dd8e997486aULL, 0x6b61f2ab55f86fb9ULL,
0xb3b234f5825fdc73ULL, 0x15ceb3899c4a13f9ULL, 0x769622e4ba2d8972ULL, 0xfaab3a367311e8cULL,
0xec1f88ecbbbe470cULL, 0x34c00900873d1ad7ULL, 0xfb494e985c4aad69ULL, 0x1980a140b1448043ULL,
0x9e908cd6c5c50677ULL, 0x222223dadb2150a5ULL, 0xdf7f3e717eba6234ULL, 0x2debfb8767536ab8ULL,
0x9e2c439c7c319710ULL, 0x1a5bc0ce85bf1265ULL, 0xb62249601a8e4c2eULL, 0x427c26a1c2e4587cULL,
0xadfc2450b50dabb2ULL, 0x597094851dc7eb0dULL, 0x3e7811247c37cdd7ULL, 0x62e909e3d4335a92ULL,
0x62623e6761d1a9f2ULL, 0x29d89d527b2b4e9bULL, 0xceb2d6853abe4622ULL, 0x7f3af22e7e59b878ULL,
0x26328054619ab88ULL, 0x55232b214e398d66ULL, 0x55c27a791184f592ULL, 0x7a02bdab153aa6e5ULL,
0xbf8f6e691874f39bULL, 0x5770b51544ba3bbULL, 0x66107385520b6116ULL, 0x150637d26cc842d3ULL,
0x375f8a59dab1634fULL, 0x6945ed05ddb86d79ULL, 0xaffbfaefef4bfd35ULL, 0x563a37be1697f59cULL,
0x38f1c70daed06339ULL, 0x2a5a1f50ee8dad91ULL, 0xe3c35804057f3f9cULL, 0xcb76925147171b1ULL,
0x5c03d64c1201aedcULL, 0x3155e3e151d04bbdULL, 0x3f01065274d92758ULL, 0x391341b6a915bd7ULL,
0xb92f810acb8c6c39ULL, 0x1b10e64d04cca44dULL, 0xcee12e195959a122ULL, 0x446d77c7071b5757ULL,
0x59e68060e5c9ad5dULL, 0x60df79f4bd97e636ULL, 0xed8e7fae1ba2d483ULL, 0x58f07f69cb09dfafULL,
0x51ef76d2b65620c5ULL, 0x333d47dc30115c12ULL, 0xc9a7e52ca8b68736ULL, 0x6f3e964fb5382a86ULL,
0x17603d955b76ddafULL, 0x8ec9af5ca65fbf1ULL, 0xcc6bd8bda9ef5b2eULL, 0x73ec07167e22db75ULL,
0x64ee60e8ff31d64eULL, 0x3d523998f9e35448ULL, 0xf4159c5ca081024ULL, 0x3d8ca2ade06c8632ULL,
0xb1674189091fb477ULL, 0x6c02323351bad5adULL, 0x6aedab0585985ed9ULL, 0x208b67b4838877e1ULL,
0xff4f26ee3d1acf4ULL, 0x7259cbff1d56f698ULL, 0x776c61dda7907f71ULL, 0x2e5cd704ab5104ceULL,
0x3eba66a59e909a3ULL, 0x2917e045d580c4c4ULL, 0x4359ec718d2c41caULL, 0x49b39f5384e42e05ULL,
0x5437da45efef1b56ULL, 0x25cd152f97ad5385ULL, 0x1b0abee80bdd6b7cULL, 0xe40a22f48fe80c7ULL,
0x5502d2a3ba27ee37ULL, 0x1963f483693107dbULL, 0xd4513dbae62acadULL, 0x12b2934745acfe8fULL,
0x3808bb26d651cc8dULL, 0x24e4ddff7ff94e24ULL, 0x395e645a21fc0e5aULL, 0x2fc00bc270e7ababULL,
0x34e7c0e716d46307ULL, 0x14d255b050f0dbf1ULL, 0x508d2708412b1ba0ULL, 0x6340027b453b04c5ULL,
0x8cd0a1ad176e164aULL, 0x384e7f74942342c1ULL, 0xf2f5744dccb465d0ULL, 0x25f08533bc61a4d4ULL,
0xcb2aef0ecb1c819dULL, 0x579331bfe13f511dULL, 0xd6d378bf36815344ULL, 0x38de204f09dd6f31ULL,
0x36c8d0741dbe5e7fULL, 0x3ab41703d14d0776ULL, 0x20eeeb1582570b35ULL, 0x5fb3944b7acc2873ULL,
0xbc729ff7efb0c465ULL, 0x570b353862719dbeULL, 0xeb7d53ebb21152b4ULL, 0xeded6ef0cfcf125ULL,
0x3d47ae679c8b8470ULL, 0x4fd93bbc4bc909a2ULL, 0x16ef7f5c73583034ULL, 0x1b6923d7e08011ddULL,
0x671f03fe9bdc6e72ULL, 0x3202882d4fc1e302ULL, 0x7c454b4c96c2b651ULL, 0x5d9fe42b63c1325bULL,
0xf1fa661dd1b28140ULL, 0x67be619395ffe86ULL, 0x503553bff5eda7e1ULL, 0x394d337494db1f90ULL,
0xc276dc5de2bb2eb5ULL, 0x1b4fc3b41fb2701cULL, 0x1532632a0fcd3093ULL, 0x5d9183e4cc742f6aULL,
0x6ee17cc6dffde72aULL, 0x37c4a82f82399880ULL, 0x954dbc9233641510ULL, 0x16fb4ce792bca53cULL,
0x8941e7088b9d5c6cULL, 0x5f38a74660967984ULL, 0xe37f0d2ed0950bceULL, 0x187a1f2ef9b84a16ULL,
0xf25873678cfafc06ULL, 0x4c8a1c8fb124cd35ULL, 0x8432bd7cffc664eaULL, 0x7bc7aa364883c656ULL,
0xb6e6788a267cff05ULL, 0x45b4ebc5ac9d2a5fULL, 0x5f3ac8e8dac20cc4ULL, 0x3d260229276fefc1ULL,
0xcc4ddc09292e54d2ULL, 0x3380665f9033cb72ULL, 0x18ff9bf72155c8f3ULL, 0x20d3c684bda10482ULL,
0x190cca4e1cb0ab36ULL, 0x797544c95799271bULL, 0xcc8b014d3bd0fc82ULL, 0x5da9f99694a4c177ULL,
0xca23efb36ede6fb4ULL, 0x3bb76e1f75b82023ULL, 0x9c00bf4248a7f3b3ULL, 0x49e737cdac3391c6ULL,
0xec557930977d32c6ULL, 0x4f2af92eb642e82eULL, 0xa5809ac5b9eaa0d1ULL, 0x272981132b34c085ULL,
0xcbcd6a6ef32168bdULL, 0x9f297e94a102f01ULL, 0xabc7dc50bd56d190ULL, 0x7d5b07de3c69778fULL,
0x36ccfd6e6af7fdefULL, 0x3717cae5050048f4ULL, 0xc1bae341106e51d1ULL, 0x4362331a5503fb14ULL,
0x5848f6eb3ab36398ULL, 0x3bc30829a731c891ULL, 0x8a9a6ab8fa4ff70bULL, 0x218ea0cb4363f11aULL,
0x101e35996866b947ULL, 0x4ca3a282eadce03eULL, 0x80e642de6767c6afULL, 0x4e8b9913aa02d9c6ULL,
0xce4f2d69f03c2d87ULL, 0x20f6b7943432c3bcULL, 0x870b62ef703b76a9ULL, 0x576536f1e8bbedebULL,
0x5ee76389c4bb07e9ULL, 0x454138f29164f25bULL, 0xa1f7ebc27d752381ULL, 0x1825050c6cf66971ULL,
0xd4422d1835b7117dULL, 0x4cf8463dff47b5e3ULL, 0xacb5f6b72c6924eaULL, 0x4af14c63ae9e1ef2ULL,
0x468fc80881526fa1ULL, 0x34dfe6c12b184f6eULL, 0x623b4873865d9d8eULL, 0x7dea1a7ecce9d1cdULL,
0xf0ebb0f1bd2a476fULL, 0x67ce8c26064d7dbaULL, 0x862ee104a5220ca3ULL, 0x1e7b8cf219ebab49ULL,
0x9a377ff1f2998831ULL, 0x78409decef0b6172ULL, 0x3a5e1286fdc57f4aULL, 0xd669171802d22e4ULL,
0x61628de00dbfddf1ULL, 0x526408684467a739ULL, 0x553fc23f0f45fdd5ULL, 0x1ddfc6d4b98ae684ULL,
0x95aaa873c526c7dULL, 0x7b6aa08e3a424b30ULL, 0x921ac50bb992f831ULL, 0x65bd44e73a5b70c1ULL,
0xa45e16184cbfe242ULL, 0xcaca6af87825cb3ULL, 0x73c5207b5a3d69fdULL, 0x595b5679badad48ULL,
0x89266e8aa2ccad07ULL, 0x3df2cea734317612ULL, 0xed8a7e69099d7942ULL, 0x59e68f62319cf3d6ULL,
0x81f5135340b519fcULL, 0x51f89d8996845116ULL, 0xd9fead0e06f9cdbULL, 0x2c651aa0492fbcc2ULL,
0x3838823fba33b2eeULL, 0x3468b6dc3f012122ULL, 0x5c72ae1fc2b2e3faULL, 0x2f619c907db3f238ULL,
0xc765ff1bc66ac21fULL, 0x1cb878648799a201ULL, 0x498d00f0f9922991ULL, 0x7e8b5b7570627aacULL,
0xf84d97ed2f17b9d0ULL, 0x6ac3924a975b8d1dULL, 0xfcaed1800ef35beULL, 0x653128b241147f9cULL,
0x3cee7d2b4281cf00ULL, 0x47de2bf37490f127ULL, 0xed9e2ede42cb5928ULL, 0x5bda7ce1cad9a0b6ULL,
0x551d2c8b7acfb31bULL, 0x34ecc65f6587545fULL, 0xed9e23f9400c16e1ULL, 0x69edde1de24625acULL,
0xf081929a0c6a9b54ULL, 0x5772e5eef89b15f9ULL, 0x9c67a478dc02f333ULL, 0x25fe9ec860e4dd9cULL,
0xa48659a24288b4a7ULL, 0x21fabe29fa902501ULL, 0xd25a33c8edb2291dULL, 0x672c320623be6227ULL,
0x74618a9a119c2ddeULL, 0x6778d2ba475f1dfeULL, 0xd8b8826d62ce5022ULL, 0x7d67efe5b1f77ec7ULL,
0x35522284330b4989ULL, 0x7c65f64afa4d71edULL, 0x8214ffce07678411ULL, 0x4db895568cd4ae37ULL,
0xd0ce66678a199e50ULL, 0x508640c07445d4ebULL, 0x5893454215b7b028ULL, 0x87f1ae04074bb3eULL,
0x508c9fecf6fde67dULL, 0x1b3a4c876cb671cdULL, 0xb463efa3f7c86895ULL, 0x3eb568b98111b839ULL,
0x2e76f4ed88b2d23cULL, 0x6afc915f4d57bf14ULL, 0xb6988bdedd0b38d5ULL, 0x30c9f0fde00143e1ULL,
0x962ec0ce6ce9a6f8ULL, 0x58a8e1fbb4dd829aULL, 0xb26feb656c00e00dULL, 0x24535cd8083a2bdeULL,
0xcd222f350b18ae3dULL, 0x43194ec2c217c3beULL, 0xf2034409cb474564ULL, 0x22c93303ef02d8e4ULL,
0x6de33f5b20b4fd10ULL, 0x334254b345cba419ULL, 0xa7a782ba7c34f18fULL, 0x4fb049fc65cda83bULL,
0x20649a2d1dc0c5d3ULL, 0x3614b2c28f29d338ULL, 0x6069b5dc9e3012f2ULL, 0x22e1deab31e41833ULL,
0xc0d901e907cb180dULL, 0x1495d95c938e47c1ULL, 0x8f2fe1a6694ef611ULL, 0x326b83b598e5e782ULL,
0x556eba1fa56170ULL, 0x277c739d675d1b84ULL, 0xe1b6217937aca298ULL, 0x3d89b77e3a2a1eb4ULL,
0xc1a533a2fec647baULL, 0x15604e0ddc0f026ULL, 0x3e5592ae05b2afabULL, 0x27ecc4f05622329bULL,
0x6696c2fb609731fcULL, 0xb99ff654c97b6c8ULL, 0xa423b37da927116fULL, 0x7e5d710a298a712ULL,
0x61b49cd3ba3293fULL, 0x46759709652aff99ULL, 0x9d8d9257e7126bd2ULL, 0x381c3998e5b69f93ULL,
0x5b02a336eef80438ULL, 0x4f37064455c48655ULL, 0xa1e194a044a333baULL, 0x199b29b879b82f22ULL,
0xd568673138af7c91ULL, 0x3cdd83c1ff34e662ULL, 0x3c8f1c897f152594ULL, 0x5776215896d57292ULL,
0x5308acdb7d5f685bULL, 0x36270e6e525f58e5ULL, 0xbbe1d0fb9892b799ULL, 0x49e2de87c2dc4e7ULL,
0x497dae2a68302171ULL, 0x77e2199c889b09bfULL, 0x6b212763e88c748fULL, 0x6ac11144a6be4691ULL,
0xcdca717f01fea4eeULL, 0x7586e200343fa1bdULL, 0xd64868d39d41fb8dULL, 0x2ea49e4633ec46adULL,
0x6d42b9659f9eb87eULL, 0x2c2e3dcb4a04db36ULL, 0x222c11998d8ee6b4ULL, 0xf0aec1d78e1f539ULL,
0xf11ba40934c1aa8eULL, 0x130c50cc43ca5685ULL, 0x6eaed92bb974fe77ULL, 0x61fece5ac127b249ULL,
0x600700d780442e53ULL, 0x67ff03acf5a52b3cULL, 0x17861c3eeaac2966ULL, 0x677dc5a84e382a28ULL,
0xff5756737023ebd1ULL, 0x67ef4507a8d317cULL, 0xedcc4797d0cabf20ULL, 0x376cbe5e493facf9ULL,
0xaffe84737b31ffffULL, 0x174bbec480e36788ULL, 0xd14c942988102e3bULL, 0x100c7d0500f64c56ULL,
0x454a5dc59b499c3aULL, 0x11c2df78d3a931bdULL, 0x9feab3dc5a3a25bdULL, 0x63e8e870bcb91816ULL,
0x444d28e0cbf57623ULL, 0x18ec6887cffd950eULL, 0xf4ba22ee98078a60ULL, 0x1fb982968965d0fbULL,
0x31d53797f8041f6dULL, 0x4518d2def1ccb1a5ULL, 0x510338e692af7266ULL, 0x3387419c14422aeeULL,
0x91550cb37f0f6872ULL, 0x2173fae73f0b0aaaULL, 0xc9d65445eb9c7713ULL, 0x6c96c9fbc3167fdULL,
0xd83ff549333adf0fULL, 0x46e775c6203a7b48ULL, 0x608b5bbc62d43070ULL, 0x51ee505d489746fbULL,
0x4d6bb340d30209d4ULL, 0x22479fb4f598a90cULL, 0x7afbd4885fe52e46ULL, 0x4ff585dcd116041bULL,
0xef48edfc50107acdULL, 0x310e0c973d804044ULL, 0xe2706f3b9d1a3b1ULL, 0x57700c9227ea5380ULL,
0xc8ef2e1f82754deeULL, 0x7d0eb83f36d193eULL, 0x1399d72f41441001ULL, 0x3fca41bd27d7dab4ULL,
0xdef7378a0e4ac775ULL, 0x7c93c36a1cb2bc46ULL, 0xfb8176460652929ULL, 0x787da2f83ae2e7f4ULL,
0x9aeaa419b358604eULL, 0x2bed86175bee99eaULL, 0x7237db803d95095bULL, 0x5883ce5d43212d66ULL,
0x1f2e26e74d56ad69ULL, 0x661034047d5bf4cULL, 0x187b323517490284ULL, 0x2903942f338916aULL,
0xaf91d5ec2752e8f4ULL, 0x229591c7b3b03eb1ULL, 0x579d497bfa5ae01dULL, 0x4a84210ead4bf234ULL,
0xc6303c274146f98bULL, 0x3530728d5cba4e08ULL, 0xbe944076f6d94bd0ULL, 0x19605b6b092a8135ULL,
0xa3633dd488ea13baULL, 0x58c8ae4966dc8e63ULL, 0x20218d824066e4acULL, 0x575fce28a63b354eULL,
0xba8f7571e2520c38ULL, 0x50a85acb476d3123ULL, 0xdcd3cac0b2a9141dULL, 0xebd69b7f642bbf4ULL,
0x2c94827c37a9867fULL, 0x4bd692b25dba05d4ULL, 0x920ff9fe6696562dULL, 0x25bff65344c1edd7ULL,
0x1cee95eb08cde0bbULL, 0x731a5d56a8ed6196ULL, 0xeaf8ae53361b0232ULL, 0x7b45e2e2250c6eaULL,
0x4bfdd3a259320ff5ULL, 0x538d83b007ce401dULL, 0x9296ed89207a7b5aULL, 0x4a0f55a24b95532ULL,
0x1e620f065f6a44ceULL, 0x22aec951255050fbULL, 0xf3bde595dbc0b177ULL, 0x4f7ed6573d90190cULL,
0x86b3e046c85f95cULL, 0x37b79ae41db9951cULL, 0x12e34050c314b0fbULL, 0x2181cc2c7a6798fcULL,
0xecb8f9a8ec2b4e6fULL, 0x641f5e021f62e062ULL, 0x187bdcc5ed8a511eULL, 0x4a72b988c3b115e9ULL,
0x7e7d929656b8565dULL, 0x5d4c584c14482380ULL, 0xc13beff4bec5fcfdULL, 0x59403408a00d5dd3ULL,
0xc0e49387acb57b76ULL, 0x342b427eb0794eULL, 0x5a910c174fc1d627ULL, 0x7ae8f446eb7c4586ULL,
0xc9c85b1b23dcb561ULL, 0x12bd7c53ee30fd82ULL, 0x63e79f0ff7ebbc78ULL, 0x54773d67650bd0a0ULL
};

static const u64 PRECOMP_TABLE_2[] = {
0x3dee5bb295508114ULL, 0x12ae82ddc97f6fcfULL, 0x60f5c1e2f5beb566ULL, 0x3f99172a63932f0cULL,
0xe33eff8dbdb66890ULL, 0x139291ca41bde4bbULL, 0x34c0b221c953415bULL, 0x5a934ebf6b24fb58ULL,
0xf197f1de2d1467b1ULL, 0x3aa3c12734d1e9efULL, 0xf08498d52a27ceb5ULL, 0x3b5fe12d9ced696aULL,
0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL
};

static const ecpt_affine *GEN_TABLE_0 = (const ecpt_affine *)PRECOMP_TABLE_0;
static const ecpt_affine *GEN_TABLE_1 = (const ecpt_affine *)PRECOMP_TABLE_1;
static const ecpt *GEN_FIX = (const ecpt *)PRECOMP_TABLE_2;

static u32 ec_recode_scalar_comb(const u64 k[4], u64 b[4]) {
	const int t = 252;
	const int w = 7;
	const int v = 2;
	const int e = 252 / (w * v); // t / wv
	const int d = e * v; // ev
	const int l = d * w; // dw

	// If k0 == 0, b = q - k (and return 1), else b = k (and return 0)

	u32 lsb = ((u32)k[0] & 1) ^ 1;
	u64 mask = (s64)0 - lsb;

	u64 nk[4];
	neg_mod_q(k, nk);

	b[0] = (k[0] & ~mask) ^ (nk[0] & mask);
	b[1] = (k[1] & ~mask) ^ (nk[1] & mask);
	b[2] = (k[2] & ~mask) ^ (nk[2] & mask);
	b[3] = (k[3] & ~mask) ^ (nk[3] & mask);

	// Recode scalar:

	const u64 d_bit = (u64)1 << (d - 1);
	const u64 low_mask = d_bit - 1;

	// for bits 0..(d-1), 0 => -1, 1 => +1
	b[0] = (b[0] & ~low_mask) | d_bit | ((b[0] >> 1) & low_mask);

	for (int i = d; i < l; ++i) {
		u32 b_imd = (u32)(b[0] >> (i % d));
		u32 b_i = (u32)(b[i >> 6] >> (i & 63));
		u32 bit = (b_imd ^ 1) & b_i & 1;

		const int j = i + 1;
		u64 t[4] = {0};
		t[j >> 6] |= (u64)bit << (j & 63);

		u128 sum = (u128)b[0] + t[0];
		b[0] = (u64)sum;
		sum = ((u128)b[1] + t[1]) + (u64)(sum >> 64);
		b[1] = (u64)sum;
		sum = ((u128)b[2] + t[2]) + (u64)(sum >> 64);
		b[2] = (u64)sum;
		sum = ((u128)b[3] + t[3]) + (u64)(sum >> 64);
		b[3] = (u64)sum;
	}

	return lsb;
}

static CAT_INLINE u32 comb_bit(const u64 b[4], const int wp, const int vp, const int ep) {
	// K(w', v', e') = b_(d * w' + e * v' + e')
	u32 jj = (wp * 36) + (vp * 18) + ep;

	return (u32)(b[jj >> 6] >> (jj & 63)) & 1;
}

void ec_table_select_comb(const u64 b[4], const int ii, ecpt &p1, ecpt &p2) {
	const int t = 252;
	const int w = 7;
	const int v = 2;
	const int e = 252 / (w * v); // t / wv
	const int d = e * v; // ev
	const int l = d * w; // dw

	// DCK(v', e') = K(w-1, v', e') || K(w-2, v', e') || ... || K(1, v', e')
	// s(v', e') = K(0, v', e')

	// Select table entry 
	// p1 = s(0, ii) * tables[DCK(0, ii)][0]
	// p2 = s(1, ii) * tables[DCK(1, ii)][1]

	u32 d_0;
	d_0 = comb_bit(b, 6, 0, ii) << 5;
	d_0 |= comb_bit(b, 5, 0, ii) << 4;
	d_0 |= comb_bit(b, 4, 0, ii) << 3;
	d_0 |= comb_bit(b, 3, 0, ii) << 2;
	d_0 |= comb_bit(b, 2, 0, ii) << 1;
	d_0 |= comb_bit(b, 1, 0, ii);
	u32 s_0 = comb_bit(b, 0, 0, ii);

	ec_zero(p1);
	for (int ii = 0; ii < 64; ++ii) {
		// Generate a mask that is -1 if ii == index, else 0
		const u128 mask = ec_gen_mask(ii, d_0);

		// Add in the masked table entry
		ec_xor_mask_affine(GEN_TABLE_0[ii], mask, p1);
	}
	fe_mul(p1.x, p1.y, p1.t);
	ec_cond_neg(s_0 ^ 1, p1);

	u32 d_1;
	d_1 = comb_bit(b, 6, 1, ii) << 5;
	d_1 |= comb_bit(b, 5, 1, ii) << 4;
	d_1 |= comb_bit(b, 4, 1, ii) << 3;
	d_1 |= comb_bit(b, 3, 1, ii) << 2;
	d_1 |= comb_bit(b, 2, 1, ii) << 1;
	d_1 |= comb_bit(b, 1, 1, ii);
	u32 s_1 = comb_bit(b, 0, 1, ii);

	ec_zero(p2);
	for (int ii = 0; ii < 64; ++ii) {
		// Generate a mask that is -1 if ii == index, else 0
		const u128 mask = ec_gen_mask(ii, d_1);

		// Add in the masked table entry
		ec_xor_mask_affine(GEN_TABLE_1[ii], mask, p2);
	}
	fe_mul(p2.x, p2.y, p2.t);
	ec_cond_neg(s_1 ^ 1, p2);
}

void ec_mul_gen(const u64 k[4], ecpt_affine &R) {
	const int t = 252;
	const int w = 7;
	const int v = 2;
	const int e = 252 / (w * v); // t / wv
	const int d = e * v; // ev
	const int l = d * w; // dw

	// Recode scalar
	u64 kp[4];
	u32 recode_lsb = ec_recode_scalar_comb(k, kp);

	// Initialize working point
	ufe t2b;
	ecpt X, S, T;

	ec_table_select_comb(kp, e - 1, S, T);
	fe_set_smallk(1, S.z);
	ec_add(S, T, X, true, true, false, t2b);

	for (int ii = e - 2; ii >= 0; --ii) {
		ec_table_select_comb(kp, ii, S, T);

		ec_dbl(X, X, false, t2b);
		ec_add(X, S, X, true, false, false, t2b);
		ec_add(X, T, X, true, false, false, t2b);
	}

	// NOTE: Do conditional addition here rather than after the ec_cond_neg
	// (this is an error in the paper)
	// If carry bit is set, add 2^(w*d)
	ec_cond_add((kp[3] >> 60) & 1, X, *GEN_FIX, X, true, false, t2b);

	// If recode_lsb == 1, X = -X
	ec_cond_neg(recode_lsb, X);

	// Compute affine coordinates in R
	ec_affine(X, R);
}


//// Constant-time Simultaneous Multiplication

/*
 * Precomputed table generation
 *
 * Using GLV-SAC Precomputation with m=4 [1], assuming window size of 1 bit
 */

static void ec_gen_table_4(const ecpt &a, const ecpt &b, const ecpt &c, const ecpt &d, ecpt TABLE[8]) {
	// P[0] = a
	ec_set(a, TABLE[0]);

	// P[1] = a + b
	ufe t2b;
	ec_add(a, b, TABLE[1], true, true, true, t2b);

	// P[2] = a + c
	ec_add(a, c, TABLE[2], true, true, true, t2b);

	// P[3] = a + b + c
	ec_add(TABLE[1], c, TABLE[3], true, true, true, t2b);

	// P[4] = a + d
	ec_add(a, d, TABLE[4], true, true, true, t2b);

	// P[5] = a + b + d
	ec_add(TABLE[1], d, TABLE[5], true, true, true, t2b);

	// P[6] = a + c + d
	ec_add(TABLE[2], d, TABLE[6], true, true, true, t2b);

	// P[7] = a + b + c + d
	ec_add(TABLE[3], d, TABLE[7], true, true, true, t2b);
}

/*
 * Constant-time table selection for m=4
 */

static CAT_INLINE void ec_table_select_4(const ecpt *table, const ufp &a, const ufp &b, const ufp &c, const ufp &d, const int index, ecpt &r) {
	int k = ((u32)(b.w >> index) & 1);
	k |= ((u32)(c.w >> index) & 1) << 1;
	k |= ((u32)(d.w >> index) & 1) << 2;

	ec_zero(r);

	const int TABLE_SIZE = 8;
	for (int ii = 0; ii < TABLE_SIZE; ++ii) {
		// Generate a mask that is -1 if ii == index, else 0
		const u128 mask = ec_gen_mask(ii, k);

		// Add in the masked table entry
		ec_xor_mask(table[ii], mask, r);
	}

	ec_cond_neg(((a.w >> index) & 1) ^ 1, r);
}

/*
 * Simultaneous multiplication by two variable base points
 *
 * Preconditions:
 * 	0 < a,b < q
 *
 * Multiplies the result of aP + bQ by 4 and stores it in R
 */

// R = a*4*P + b*4*Q
void ec_simul(const u64 a[4], const ecpt_affine &P, const u64 b[4], const ecpt_affine &Q, ecpt_affine &R) {
	// Decompose scalar into subscalars
	ufp a0, a1, b0, b1;
	s32 a0sign, a1sign, b0sign, b1sign;
	gls_decompose(a, a0sign, a0, a1sign, a1);
	gls_decompose(b, b0sign, b0, b1sign, b1);

	// P1, Q1 = endomorphism points
	ecpt_affine P1a, Q1a;
	gls_morph(P.x, P.y, P1a.x, P1a.y);
	gls_morph(Q.x, Q.y, Q1a.x, Q1a.y);

	// Expand base points
	ecpt P0, Q0, P1, Q1;
	ec_expand(P1a, P1);
	ec_expand(Q1a, Q1);
	ec_expand(P, P0);
	ec_expand(Q, Q0);

	// Set base point signs
	ec_cond_neg(a0sign, P0);
	ec_cond_neg(b0sign, Q0);
	ec_cond_neg(a1sign, P1);
	ec_cond_neg(b1sign, Q1);

	// Precompute multiplication table
	ecpt table[8];
	ec_gen_table_4(P0, P1, Q0, Q1, table);

	// Recode scalar
	u32 recode_bit = ec_recode_scalars_4(a0, a1, b0, b1, 127);

	// Initialize working point
	ecpt X;
	ec_table_select_4(table, a0, a1, b0, b1, 126, X);

	ufe t2b;
	for (int ii = 125; ii >= 0; --ii) {
		ecpt T;
		ec_table_select_4(table, a0, a1, b0, b1, ii, T);

		ec_dbl(X, X, false, t2b);
		ec_add(X, T, X, false, false, false, t2b);
	}

	// If bit == 1, X <- X + P (inverted logic from [1])
	ec_cond_add(recode_bit, X, P0, X, true, false, t2b);

	// Multiply by 4 to avoid small subgroup attack
	ec_dbl(X, X, false, t2b);
	ec_dbl(X, X, false, t2b);

	// Compute affine coordinates in R
	ec_affine(X, R);
}

