// Elliptic curve point operations
#include "ecpt.inc"
#include "misc.inc"

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

static CAT_INLINE void ec_table_select_2(const ecpt *table, const ufp &a, const ufp &b, const int index, const bool constant_time, ecpt &r) {
	u32 bits = (u32)(a.w >> index);
	u32 k = ((bits ^ (bits >> 1)) & 1) << 2;
	k |= (u32)(b.w >> index) & 3;

	// If constant time requested,
	if (constant_time) {
		ec_zero(r);

		for (int ii = 0; ii < 8; ++ii) {
			// Generate a mask that is -1 if ii == index, else 0
			const u128 mask = ec_gen_mask(ii, k);

			// Add in the masked table entry
			ec_xor_mask(table[ii], mask, r);
		}
	} else {
		ec_set(table[k], r);
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
	ec_table_select_2(table, a, b, 126, true, X);

	ufe t2b;
	for (int ii = 124; ii >= 0; ii -= 2) {
		ecpt T;
		ec_table_select_2(table, a, b, ii, true, T);

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
	const int e = t / (w * v); // t / wv
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

void ec_mul_gen(const u64 k[4], bool mul_cofactor, ecpt_affine &R) {
	const int t = 252;
	const int w = 7;
	const int v = 2;
	const int e = t / (w * v); // t / wv

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

	// If multiplication by cofactor is desired,
	if (mul_cofactor) {
		// Note that this does not improve security.  It is anticipated only
		// to be useful for signing.
		ec_dbl(X, X, false, t2b);
		ec_dbl(X, X, false, t2b);
	}

	// Compute affine coordinates in R
	ec_affine(X, R);
}


//// Constant-time Simultaneous Multiplication with Generator Point

/*
 * This function is useful for EdDSA signature validation, so it is
 * interesting to optimize for this case.
 *
 * Interleaving the ECADDs for ec_mul with those from ec_mul_gen is
 * a straight-forward approach.  We want the ec_mul_gen table to
 * stay at 128 points since that is the optimal memory access time
 * trade-off.  But, there is no need to use multiple tables since
 * the ECDBLs need to be performed *anyway* for the ec_mul ops,
 * so the ECDBLs are sort-of "free."  So the optimal choice for
 * table construction is a little different from the ec_mul_gen
 * case and we need a new table for w = 8, v = 1.  Since 8 does
 * not evenly divide 252, it is not necessary to do the final
 * correction step addition which simplifies the algorithm a bit.
 *
 * For this tuning, ec_mul_gen ECADDs = 32.
 *
 * Since ECDBL/ECADD ops are linear, it is possible to interleave
 * ec_mul_gen and ec_mul even though the number of ECDBL for each
 * is different.  Introducing ECADDs for ec_mul_gen near the end
 * of the evaluation loop of ec_mul still exhibits a regular
 * pattern and will just require another 32 ECADDs.  The final
 * conditional negation from ec_mul_gen can be merged into the
 * ECADDs by inverting the sign of each added point instead to
 * avoid messing with the interleaving.
 *
 * So overall the cost should be about the same as one ec_mul
 * with just 32 extra ECADDs from table lookups, which falls
 * about mid-way between ec_mul and ec_simul for performance.
 */

static const u64 PRECOMP_TABLE_3[] = {
0xfULL, 0x0ULL, 0x0ULL, 0x0ULL,
0x36d073dade2014abULL, 0x7869c919dd649b4cULL, 0xdd9869fe923191b0ULL, 0x6e848b46758ba443ULL,
0xe692c2b52108689ULL, 0x29b1db37cf0df068ULL, 0x81b333cd2ec495dfULL, 0x2f866bf8ceb54fccULL,
0xed33e9d7b8cfb7f4ULL, 0x77be0cb0857980c4ULL, 0x73621b5857122ccdULL, 0x704a4bc8a79fd8b2ULL,
0x79e4e13c64a01e42ULL, 0x22d804171cf98ef1ULL, 0x28ae8217fb37d6ddULL, 0x2aa9137894096adeULL,
0xa12a7c195bccc055ULL, 0x60327de2d40e4d4cULL, 0x4d463cba7511e06ULL, 0x1af8959535cd0362ULL,
0x29e3aa9143bc1149ULL, 0x713965bea15c4f5bULL, 0xeb4a60a973865558ULL, 0x7762bd993319d66aULL,
0x41cd4eceaa0d0da8ULL, 0x733689b6aee3528bULL, 0x8bd8c3141f568387ULL, 0x1516b2ba94ebd332ULL,
0x63c514febb2d96dfULL, 0x218972a352400134ULL, 0x9dda12f63e3038aaULL, 0x66c63ac9c140acf2ULL,
0xb810067978c37089ULL, 0x2db36d30438e013dULL, 0x4eb3d98b82f4e828ULL, 0x101267a0e0b4597ULL,
0xe73c64dd9dce1560ULL, 0x70140cb6949f6d25ULL, 0x3d51ec7c1f1a0b77ULL, 0x2e6f36168b557585ULL,
0xddf30c5bc8767e6dULL, 0x2e35bc82a1367d92ULL, 0x2bbbb168f7bee97fULL, 0x542cb4eeb0bce77fULL,
0xa7cb9837642326d1ULL, 0x16dc1d224e2b84c8ULL, 0x76f7e18370ad322ULL, 0x1ce58585c5dbad67ULL,
0xd87a29f07eae34b0ULL, 0x417d7d06daae20afULL, 0x7933be8ba4020c46ULL, 0x511c3a2d5a28e0aULL,
0x638170ef259da988ULL, 0x5bbac242c4583b35ULL, 0xa2abdd183a3c195aULL, 0x72bc3b185027c130ULL,
0xf261b4c10931ff95ULL, 0xb52ae0632cf5a57ULL, 0xe170c99fe22a4f6bULL, 0x4aba1faf1c816acdULL,
0x5eaad107f9e358beULL, 0x77b2d4a60db4260cULL, 0xb53790ec441e8edcULL, 0x2c0bf86097e287d2ULL,
0xe881cc00ebcf9c0ULL, 0x4751a9131ff5c6e1ULL, 0x8485ba0f638daf7bULL, 0x7e600f67a2566517ULL,
0x6c3498f39756606aULL, 0x2b1c7d0fac019c98ULL, 0x43eebe827bb163a6ULL, 0x57e94fb4f3a3eb5aULL,
0xe58349c9ba104ecULL, 0x2199b74f4e6f315bULL, 0x9e6363e941229c88ULL, 0x2eaef58807c3eb2cULL,
0xf2934b1549f03c4aULL, 0x694c49d8062ada4fULL, 0xa73a37f92b014459ULL, 0x5dd01b27f8002943ULL,
0xc65defded541ec91ULL, 0x1c22ff7825a64c97ULL, 0x52ef08df835147abULL, 0x4c6024f28a176b5fULL,
0xb91732c5f37a8b22ULL, 0x5baba2f3ec3f4d7eULL, 0xbe23511b5f265fdfULL, 0x37bbf6921f553b16ULL,
0x49f5740fac11b32bULL, 0x61dbddb3fe57f683ULL, 0xa814e93b287f7cc2ULL, 0x258c063807c85377ULL,
0x30f55d8e3d94bd7eULL, 0x3b87d935125f0fa3ULL, 0x930aebe79cc81956ULL, 0x621df22a74fac617ULL,
0x3eccaef5ba7dbab8ULL, 0x494517814b8073b5ULL, 0xb0447e592b1e5202ULL, 0x6bcb37791046ce10ULL,
0x2be8e9fd50a42224ULL, 0x6996f9db6cd065b1ULL, 0xc386ac61bca15e3fULL, 0x48b9efce0202e328ULL,
0x80c10f28e75c2871ULL, 0x76994c74b2639ce3ULL, 0xde07c0e82f2c17a7ULL, 0x19149653b7311706ULL,
0x1d30252e26d4134cULL, 0x1bd243b921953723ULL, 0xa099ff71e9e34897ULL, 0x70983a37615f90d0ULL,
0xee87840f22c8b68fULL, 0x96c7e80297bebf3ULL, 0x2f6855ed63f557ccULL, 0x4da4d2f761153048ULL,
0x399d816c2e082be5ULL, 0x4ec14f9b72eb3a4eULL, 0xeb951c5a672313e7ULL, 0x41bd029fd2b56d6eULL,
0x94a44ad8febbfd4cULL, 0x58437e7dc9ead019ULL, 0x88e7827a66cd9513ULL, 0x496cbe2c640fe578ULL,
0x50aa74cc8ae93118ULL, 0x73d4bf8aeb440efcULL, 0xbea7c962be371d0ULL, 0x2884b2e59dac0899ULL,
0xe8d5352bf53dedb6ULL, 0x24ae13d6fd2c22d9ULL, 0xa89aa05c55de69f3ULL, 0x8fe82c072dd4e67ULL,
0x3191010895c6be96ULL, 0x311a47355bac5c7eULL, 0x3472d3bb023b3c2bULL, 0x2599527b1275748dULL,
0x601cf96330cd6a25ULL, 0x70981e92610b61d0ULL, 0x6a47abcdbfb3b5caULL, 0x755b7d5eb445ee1bULL,
0x3c8d4e871770d6eeULL, 0x6d46d0eea417355ULL, 0x53e823a0f7ff1484ULL, 0x470642ea232463a6ULL,
0x72a175ab41831a7cULL, 0x31395f627a8b27c1ULL, 0x6c7d1a2c09702db1ULL, 0x79fab6f1474aaf47ULL,
0x3b1a58ff0922aff3ULL, 0x4f0a25074400cef4ULL, 0x8a8ee482d04055a5ULL, 0x49cda7f84c57166ULL,
0xfa5a9155dffeb757ULL, 0x363067500a4f84eeULL, 0x5447c2a78172c167ULL, 0x7c102c86f29164a6ULL,
0x8ceda0c1786a3e15ULL, 0x1891b247d096db4ULL, 0xcc5ce4be544a1960ULL, 0x7176086423552376ULL,
0x19f8055f631df57cULL, 0x20e073de066f9894ULL, 0x76c71d20632f2562ULL, 0x1aa41f4d0664cab9ULL,
0x2be06a00d326bb1fULL, 0x3e47da7957832c24ULL, 0x918e3f5a7b36ed79ULL, 0x242106dd9ca01031ULL,
0x18a0697d42afd536ULL, 0x11f33fdd622714c2ULL, 0xb1d7be810c06bed7ULL, 0x7b9cc40626457058ULL,
0x3da028baeb632baeULL, 0x3de075d3105d8285ULL, 0x9aecb8d882b40092ULL, 0x73175b6a8357ad0eULL,
0xcedb306d65f0be19ULL, 0x128efb20b310967aULL, 0x3d629fd265c50970ULL, 0x78e82c99bcc37ebcULL,
0x5bd73bb3f57cf051ULL, 0x32cba8a9eb098243ULL, 0x7bb072abb066a2d5ULL, 0xf221220013e8e14ULL,
0x2fcdce462f930795ULL, 0x402c9d364043be20ULL, 0x45ecf791782d2acaULL, 0x202904d26441eb30ULL,
0x32e3f7a1be59b95dULL, 0x4de0cd47b63960b2ULL, 0x1324d540401ee759ULL, 0x195f5e3e2ff65a8ULL,
0xa9533dde52f503d2ULL, 0x6c0fd23d14ce520cULL, 0xdaecab86a0a008adULL, 0x1b0a22e40819cfebULL,
0x6ee5e9988fbb95b5ULL, 0x6ef95f18b20f19d7ULL, 0x1f2c22e8a3f08953ULL, 0x6c90ca047d1811c5ULL,
0x8537610fc2fcb9a0ULL, 0x6c1fdd548c4c93d1ULL, 0xda06b44b53bdc79aULL, 0x645929f91d3cc07eULL,
0x82c99f77a1381b4ULL, 0x6405c786adf01d82ULL, 0xcb5969eacd3558c3ULL, 0x2132f36ecf5d64efULL,
0x46cc872e464fda06ULL, 0x236a92bba26eff9aULL, 0x3c293cb0606dc04eULL, 0x41716fec61af1c29ULL,
0x2dc4670ba4e564e4ULL, 0x2153ab5c05bf5cbULL, 0x4eb260a192a68ad6ULL, 0x4549608d42f230dfULL,
0x8e61a754985402d0ULL, 0x69cf45006e98860aULL, 0x601f108e9632a4bcULL, 0xaa73604993a0731ULL,
0x7ef954d84705d87cULL, 0x6fd19ddaa16f1375ULL, 0xa398968b090f7fb6ULL, 0x4a3e1908269f4735ULL,
0xe916703f6fdd6ccaULL, 0x60bac745a363f431ULL, 0x9c5a003699a1153ULL, 0x7c64904469d3ee4cULL,
0x96c9f2f6bd228058ULL, 0xd8dc5deb258664dULL, 0x362fa8a4fda8091bULL, 0x6022b605e799330dULL,
0x1558ea9c651231d0ULL, 0x2adc35c359863ae7ULL, 0xd7784bf9de06d2c9ULL, 0x2c45edc058c6c07cULL,
0x2a272f140cf6a039ULL, 0x6cf5567783b932bdULL, 0xb366286876763a09ULL, 0x1ed3fb5745b102a7ULL,
0x39bfa0406ec1d258ULL, 0x5340f69d79366580ULL, 0xb16e1482f5239340ULL, 0x6f7a614a5861efb1ULL,
0x8a873fdb7081d704ULL, 0x25011a5393ef5052ULL, 0x2b8ed62a09412261ULL, 0xa511bd753027e4dULL,
0x6ffb58e77b1a075bULL, 0x1a0997c83a30f3e2ULL, 0xfac380fbdf73ec1bULL, 0x5347439c7e936e10ULL,
0x69c17ab675db020dULL, 0x5af4805f99a172daULL, 0x6e6ea3b002befd41ULL, 0x2122d08c68606296ULL,
0x934501006dd12750ULL, 0x1998299f076df9a1ULL, 0x92a3c5b9b07319bdULL, 0x6e1da08c14284de4ULL,
0x46896b4b00abb531ULL, 0x350d6776c51932a9ULL, 0xabb0e649510ec69cULL, 0x27d59f94a1013afcULL,
0x9d9efedec634ba37ULL, 0x78c2d703a4b55c48ULL, 0x54f79179cccfb0b5ULL, 0x1f49f591c2c63785ULL,
0xc83f76814420a8a9ULL, 0x34613ff8436bffe5ULL, 0x3ee8f858a065ac74ULL, 0x3474de16e38f3ac3ULL,
0xad640e4bc9ade2aULL, 0x443a14c5568b66e5ULL, 0xe016f7c12189001aULL, 0x4d71fa65db967e48ULL,
0x2840420705bfa45eULL, 0xa0cc8cf1645138fULL, 0x51f631c9928753e4ULL, 0x5dad8874f340c7acULL,
0x7115aee0861b37dcULL, 0x65783a8443495b0cULL, 0x3e03cd00994c9129ULL, 0x202b0bacba8ac5acULL,
0x40dce2cfa0d1ca2cULL, 0x1bd9b233ae2fb858ULL, 0x309156b972266e96ULL, 0x66a0e17e938789d5ULL,
0x87cf9baf62880e94ULL, 0x114e5f0f1fe1b829ULL, 0x609c915831c3c07aULL, 0x1b21919950fd249aULL,
0xda758fbf680f83e3ULL, 0xb3a08e40a24b8c2ULL, 0xdd777684f4383b02ULL, 0x47bac7465e9ceaf7ULL,
0xd0db3ac726a1fc64ULL, 0x2b9aeed793a7c10eULL, 0xb22613cb6d9df6c9ULL, 0x1886dc74df384dULL,
0x3e91c34d7aa57aeaULL, 0x69d9bdaaf993aedbULL, 0x8422a4f24d9bb9daULL, 0x2410c1750863bedULL,
0xe30c350602f96d88ULL, 0x74ae58c3c571d036ULL, 0x197b09b142846429ULL, 0x463dd9816c666996ULL,
0xd3d4ce26fe9337c9ULL, 0x78748a39dbdd968cULL, 0x68f8765c9312987dULL, 0x170635739f8e35f7ULL,
0x69202261cff9aff4ULL, 0x189d2b7acb168047ULL, 0xa484de22faba45dbULL, 0xd75bf94230b2b2fULL,
0x112af14864ea5a27ULL, 0x779d50a67c9c72ffULL, 0x26313c47d178df80ULL, 0x4d659bb726593146ULL,
0xbf1c32db7cfa7a2cULL, 0x79397fe8ccdd49fdULL, 0xc4ba634f14ee2c1cULL, 0x4db880a7c6752efbULL,
0xc2a7df12c7b9f8f9ULL, 0x5278bccf4dbfaa52ULL, 0x4062efe069b7efe8ULL, 0x3f4387f9e049b11bULL,
0xb1e2392f05b688b1ULL, 0x1cd8b2428921af0fULL, 0xf8b51ef15c8cd5e6ULL, 0x143a4275e683d1baULL,
0xe98ea99ee13fcf6dULL, 0x6d567f2a1887e8baULL, 0x4143e2e109d13320ULL, 0x626c20c162702c91ULL,
0x68967c0cdae4aca1ULL, 0x5bd98faed4ecc0bdULL, 0x97d47540514d24e6ULL, 0x7e6b9d24fafe80f8ULL,
0xf5bff5a891ede735ULL, 0x600dbbf318966aa8ULL, 0x62a13cb7c782c7eeULL, 0x396c66e698820d80ULL,
0x9d897d208ab10349ULL, 0x5989e13a92dd696bULL, 0xf1ce1c7713bee7f9ULL, 0x509b0ad491c31784ULL,
0xfa2d57c00f140aefULL, 0x6c71549ad4f91cf9ULL, 0xff9cd3b801351c82ULL, 0x2a46f6a7e42b835ULL,
0x4047df78b8b7dc5eULL, 0x6d163e3ae85a22a1ULL, 0x90fd2f56f0330665ULL, 0x1ac33a648bffd5e8ULL,
0xa572edcc658e01e0ULL, 0x60c3ef82559cc2d6ULL, 0x4d07761af11c66b0ULL, 0x70113bfb53d7e291ULL,
0xe4a79aee160f24dfULL, 0x1de9a052171c5fdULL, 0x49ad4d9775bba82fULL, 0x5fa25a756df35a26ULL,
0x7aaf1d4392b93f28ULL, 0x766f697bd1e055d1ULL, 0x66279bd01ef815ecULL, 0x72213eb562d07af3ULL,
0xda01c15067aa3166ULL, 0x7753ad7c1b5254d3ULL, 0x4632802573117ceaULL, 0x6c32c3d056988ff7ULL,
0x6c0b502a27f9317aULL, 0x26a5d817f2ac2859ULL, 0x88db19c4c2575c17ULL, 0x6160ad345ce1911dULL,
0x4b616f28ef1acdebULL, 0x346566f4b1f23520ULL, 0xd6c12b5bd9f0d5a5ULL, 0x285236c467e1b69eULL,
0xee5022f45d8884c2ULL, 0x648187a66cee4120ULL, 0xd21b231d8f11cfe9ULL, 0x60d80627a007e5dfULL,
0xa5fed209cf08a3ffULL, 0x734afb90aecebd0fULL, 0x72d10ba167f13809ULL, 0x323904255b4ecbdaULL,
0x1b1ab20ba2aa5309ULL, 0x16e5cb7b70bd6b15ULL, 0x9d49dbbf898068daULL, 0x610e7f28d5c8773dULL,
0xbd2ac30cb0e05179ULL, 0x4903b9dbbe2f7bddULL, 0xb82d06dc7d78430aULL, 0x531d540c811dcde5ULL,
0x9b63650c014103ULL, 0x2825fdb792f4b3c3ULL, 0xb9111d805ffd66f7ULL, 0x180f5016eff97ab6ULL,
0xf597e6e466abc784ULL, 0x2e5e8507e1bd4fa7ULL, 0x19a5c489e0182e26ULL, 0x63f2c26ce5a72e40ULL,
0xa5dc9673ff9a79c0ULL, 0x7f8255bf092ef306ULL, 0x7e155d2866cd6a89ULL, 0x32f4488cee7cc92bULL,
0x45abdc745e099457ULL, 0x48330ed29d5fb570ULL, 0x2127414334b981d3ULL, 0x73520059f7930049ULL,
0x19772b40e5d25214ULL, 0x6df2fc1900010d6eULL, 0xd3175b2943ea940fULL, 0x5159c01c16002f1ULL,
0xb153e0fde9c68466ULL, 0x37316c19540738b4ULL, 0x17c7a66da0e63fa9ULL, 0x4454853e29968806ULL,
0x23ca0c21fcf80990ULL, 0x5a2fefc1b1a6e898ULL, 0xd8a24d5cdc997be6ULL, 0x46f64a09ee78ecb1ULL,
0x36fe3f977d44f733ULL, 0x124a9da2ee5690e8ULL, 0x65328938364d2f6bULL, 0x160fbe5355822e7aULL,
0xb6f5a869b1c17b44ULL, 0x67e4de26dfe4f4e3ULL, 0xed31b78913183956ULL, 0x2c62ccbea4c7da32ULL,
0x1fd1cdeecd101730ULL, 0x2150a71f8e3a7501ULL, 0x293650ff21de4af8ULL, 0x6cc03c682d16ce87ULL,
0x81459452abcc3bf1ULL, 0x587b97f214fb22d1ULL, 0x260d195ede9a0a54ULL, 0x3364f6b106b7dd6bULL,
0xa88ad6435d617f44ULL, 0x2355f9fb9632479aULL, 0x1845dc57f85de50dULL, 0x50fbfd2d219a27ffULL,
0x52096036ab699940ULL, 0x65c9c827c29b4287ULL, 0x190953f69b9f5eaeULL, 0x12e656eb41a2955fULL,
0x58f92f8a657e2856ULL, 0x16bd0d4f1e541580ULL, 0x735123c880e2c380ULL, 0x219d0128109313dbULL,
0xfb4cf233b3e5666cULL, 0xa2f8e51fe33ed44ULL, 0xbd24dd02788e30aaULL, 0x5eda3d73d67edcfeULL,
0x9a6b72bbd5306c25ULL, 0x168be5ec834be619ULL, 0x522dd529d605c494ULL, 0x25500f1b353fc0ULL,
0x24adfba3cb9f1b3cULL, 0x2be2ba5b13e9888dULL, 0xefc0e10e4824adc5ULL, 0x70878fa33bcb59cdULL,
0xbcba9087016c1bc3ULL, 0x5b41c03d1e31f40fULL, 0xecf78c5c7639ef25ULL, 0x7c3c2dd2d495fc92ULL,
0xa651602fdb2903a1ULL, 0x61522d4b51ea963bULL, 0x2475a95eaf7f2847ULL, 0x6e0e858004c8a02ULL,
0xb2f475ff25e03693ULL, 0x7bb563e250534b27ULL, 0x26c7385856bf5eccULL, 0x134fd9925eca85b8ULL,
0xbf3692ff808759c1ULL, 0x1c155eec26e55d58ULL, 0x55cf25b67d184e6aULL, 0x3abd66502df1c2f2ULL,
0xcb7a4385c079c347ULL, 0x638ab02664b9ab95ULL, 0x85c5710c4228b303ULL, 0x290c4ff9d4bd3322ULL,
0x18b6d4b54ae09599ULL, 0x2c6b7df3c468c854ULL, 0xc259b6294eef466cULL, 0x612595765012eac7ULL,
0x2c9e88fd66e705eaULL, 0x59a8c6785421a523ULL, 0x33a3cf8452667f6dULL, 0x62bbc43954533cd5ULL,
0x97b5a6fc2a121de9ULL, 0x3536ce27d7b4db1eULL, 0x3508112afe64137ULL, 0x20df22ff1b4ffc79ULL,
0x7fe66d225f42baa1ULL, 0x448de6617aa69d87ULL, 0x4522e8a2edeb7848ULL, 0xef365c3437a8f9dULL,
0x7c1d9d5989844c0ULL, 0x13aec7211f384708ULL, 0x9bbc1f63aa0090c7ULL, 0x393223699758ccc2ULL,
0x81e4b763692f49b9ULL, 0x242897e57317126ULL, 0xa134b85714e43dafULL, 0x48b027b70199c4eeULL,
0xd0c546f80f494033ULL, 0x3c6bbf1949bdf71dULL, 0x99d331ac48ea6906ULL, 0x622c2b27914b8119ULL,
0xbb81a8a691120bf8ULL, 0x2c5c4163a7789161ULL, 0x89205cd4bbc6295ULL, 0x6126bfc0e8de65d8ULL,
0x5dacfb777572136fULL, 0x590f6251272b94ddULL, 0x73c8b73fb9ba96d0ULL, 0x66e952c0b81d25d7ULL,
0x551ad368e6d4f8ebULL, 0x41e20a74db18b275ULL, 0xf9f835169c12c678ULL, 0x4d02a8a7191f8a4eULL,
0xdd955e19782c7df7ULL, 0x63789edf7e8bc884ULL, 0xa0693736e09a368bULL, 0x2ca2c60d6ca962baULL,
0x1e6a8a9a4566c6e0ULL, 0x3c2695ea2b806e85ULL, 0x908b0371f92a31c7ULL, 0x33b25ff50d4c2051ULL,
0xb3ca398bffbe16afULL, 0x53c594da786302fbULL, 0x9715569a4d5763dcULL, 0x67673912ad9302fcULL,
0xcde3c285de6b305aULL, 0x4e823af06fc13538ULL, 0x728764f31665e9eaULL, 0x1adc884a00a68140ULL,
0x609b1028e2d255e7ULL, 0x78a570b42cdb656eULL, 0xb7c3e93e44d84b2bULL, 0x4b2b08e5ca5dd3cdULL,
0x4e229a522da57343ULL, 0x753de3c3b5e7d47aULL, 0xe45151ce80a64704ULL, 0x5cf92edd8d9f9b4aULL,
0x32e50326e3ecfc5aULL, 0x7606968908a00a7dULL, 0x4c8f067e97bed258ULL, 0x1a15b9f88545b8c4ULL,
0x3ea55674c3d69795ULL, 0x44993dfe83fca6e9ULL, 0x7d83ca39ada462d4ULL, 0x592dfb4a5e031de9ULL,
0xc86052c36bce3be2ULL, 0x76041c90ee2e0fcaULL, 0x7961d3cc1b77c5a8ULL, 0x15061c7febe100bcULL,
0xba6fe327f8e86d7cULL, 0x33850d3b47bd51a9ULL, 0xc00edabd4f822a3eULL, 0x456386ebfc0de0f6ULL,
0x8d9d2d0bc1daa720ULL, 0x59c73b7b0ca7afb9ULL, 0x7c730af762f03bfcULL, 0x46b6db79b8d82162ULL,
0x8b4b99091659762aULL, 0x144edc732fdcacc4ULL, 0x105cd77032658432ULL, 0x1026448d611f3d75ULL,
0x71258a4cc2331c2cULL, 0x4ef12183fb2c3cf3ULL, 0x273502ac8229624eULL, 0x2244ad3f7401cde8ULL,
0xb3edd5d9ecc83cd3ULL, 0x65d27166602f4430ULL, 0x32def6f96bfc94e7ULL, 0x481853e5447d6ed2ULL,
0x6c6b2a4d182befe1ULL, 0xae090e3a78cc6baULL, 0x7ab89549798c3055ULL, 0x4424e7fee68df44aULL,
0xa22a1dde61c4bcdaULL, 0x7f61e981dfb1e6fdULL, 0xea65364b1e26cebbULL, 0x14400f0827e32d8cULL,
0xadea183f733cca66ULL, 0x43ffacfeb7308026ULL, 0x556063a0ff9b6715ULL, 0x46360b2ccfbcbc16ULL,
0xab0c96189fe2622bULL, 0x9e539f6f0c3be24ULL, 0x71d64227529a6182ULL, 0x3f04d3a0be2c0011ULL,
0x474f64d5bc44cc1eULL, 0x48fc31c4e5c42917ULL, 0xcab4b6c5ada07f5aULL, 0x6765f3c38531329aULL,
0x1c44ce2fd9a508fcULL, 0x6a95c0fa18499970ULL, 0xc124f00a71811a70ULL, 0x37699f892e059956ULL,
0x11cf868f6b7f8342ULL, 0x558603795e3e8df2ULL, 0x6ed434ad61bf00d3ULL, 0x843221ba0ff3bb1ULL,
0x1fe49eb46a9e2418ULL, 0x11800f94ec4edd80ULL, 0x5b9bc04b00e5540dULL, 0x62456c4f033eea9ULL,
0x849652cf634c0c6aULL, 0x4af94ea80b16aeeaULL, 0x6cb82878e6facea2ULL, 0x77762f31499fdffbULL,
0x9151500eb0492b41ULL, 0x276a8d98a360fa8cULL, 0x3812384550fc1776ULL, 0x5d1f0b100e025e1eULL,
0x7fd455fe5e53427dULL, 0x6089f56c2a7de8beULL, 0x423e5d80dd74ce14ULL, 0x2797b91f1eae8633ULL,
0x6047e19df907e159ULL, 0x6854b34877c7cf3dULL, 0xa0f9ab9797525f6cULL, 0x1afaa0e17988cc8aULL,
0x3b0f698cec073832ULL, 0x63503857639b62f5ULL, 0x49c526f1872c7a58ULL, 0x132061288bedb713ULL,
0xe62725136479b2c1ULL, 0x1da6d569da9e70f4ULL, 0xb937116a7e460218ULL, 0x27712fb79f5182a1ULL,
0x1667a5d62e3b35f2ULL, 0x6944398b36d0e12fULL, 0xb00154c5de63205fULL, 0x2a2f74343febe92dULL,
0x973304ad6c01f423ULL, 0x5f5b8a1dc473b605ULL, 0xa356f5e8e64e7b0dULL, 0xbe8e6c0bdb2ce4ULL,
0x50f604d8d7c6e4d1ULL, 0x4cdcfc1f286ad550ULL, 0xda71782f6a178817ULL, 0x5e7bececb620cd3fULL,
0x23d46bd892ddaa95ULL, 0x3e18c33b88e3a5bULL, 0xd261e6ad57b5601cULL, 0x26674e7ba7b351a1ULL,
0x4ab0274412d283a0ULL, 0x71793f18b0407c96ULL, 0x36640a4cf9cc2020ULL, 0x33b032797bf64acULL,
0xddb390e1c70add33ULL, 0x4ea540a2e98bbfd0ULL, 0x8084ec31666ffe74ULL, 0x6c0b56e4767a2f4aULL,
0xebaf9d8854db2992ULL, 0x656dd17ad53830b5ULL, 0xb47d6faecf4de613ULL, 0x52067591ce164ea1ULL,
0xe9aa92b9daa7ce41ULL, 0x6913d4a8f78129e7ULL, 0xa8156d704a9a7aaULL, 0x345b7a148952142eULL,
0x51ff94e0ad4a4529ULL, 0x5cc29992004f3dd6ULL, 0xa37a2748454bae83ULL, 0x76dc68973b70a735ULL,
0x1732e867c05f4c3ULL, 0x10b02134132a54c7ULL, 0x58c998cc52df5803ULL, 0x19f2ad5efd2651adULL,
0xeac96c6272e08de5ULL, 0x78e199a16dae5264ULL, 0x6629b6f83bcc2bdeULL, 0x40b8e8785ea0b265ULL,
0xc637bdd0a19491b6ULL, 0x358130ec5a75beb6ULL, 0x285cf8d460b35f7fULL, 0x54c6d42d31d0e58cULL,
0x1f0dbcc9fa802ec7ULL, 0x39b32da4c36db104ULL, 0x91d40b161039a071ULL, 0x21da2e5011795f12ULL,
0xb6ed28be9c23f454ULL, 0x6c4bf1d81ece1de6ULL, 0xa1226e2387524be4ULL, 0x72b384a1b53413f9ULL,
0x892bde386b219d1bULL, 0x366ec0418b0a8c8bULL, 0xf3ecd71ef7c5b14dULL, 0x1881d4434fa8d568ULL,
0x13bbe301bf631c86ULL, 0x6f39fd8d0006c9f9ULL, 0xcd1d85102c69d7dfULL, 0x402b03a64c5c921eULL,
0x2ff37c5bf81676cULL, 0x4ac2ef8d1efba68ULL, 0xa8f5fc015227b1ebULL, 0x7a748c013c9d522dULL,
0xa8ed2a74b4f6e1ULL, 0x52c3eb3694981da5ULL, 0x263ca40cc6150d41ULL, 0x69d3deec48021ce9ULL,
0x8a923839122109fdULL, 0x278bcab3ddb7e2fbULL, 0xe686212f5306f798ULL, 0x6dda731055efa3c9ULL,
0x5285bbde66d3b37fULL, 0x65db550a960a3099ULL, 0x401d7ab5b96c639bULL, 0x7f72f954d21f3ea5ULL,
0x6183565ec48cf3f3ULL, 0x77628c047979444dULL, 0x615989afe0031dc9ULL, 0x6424ae2688d82aa0ULL,
0x3660641944ad8ed9ULL, 0x16dcb994e8c7cae5ULL, 0x38f10f9dbd924a6cULL, 0x7116064182a230e1ULL,
0x1b57bd3af144dbbaULL, 0x6ffe779348b58a8dULL, 0xc0bb30891bcace13ULL, 0x22b9c2f948feab67ULL,
0x8dd0b11b3a78179ULL, 0x44f5863174bd6543ULL, 0xcd36139618a7f8c2ULL, 0x48de27484afd3026ULL,
0x672f7964406a52dULL, 0x7d74ec938450a17cULL, 0x79be5066ae021a45ULL, 0x1b8335514f106762ULL,
0x36347273d066a6c4ULL, 0x729a07aa666b9bebULL, 0xcd5c9679b04a78c0ULL, 0x1871e463922f7acfULL,
0x84aae36bbdbe5382ULL, 0x2492c3380302086ULL, 0x7a2f7807924ea10cULL, 0x5f5ed39eb6684087ULL,
0xe4d6b596114994eaULL, 0x7b5d7b880cd514c1ULL, 0x759c7c738fc06e43ULL, 0x65c0cc5ba1eccc2bULL,
0xb3595896d47494beULL, 0x635d7dc1664e624fULL, 0xa307af0a338193b3ULL, 0x75f2a87e5f945e7aULL,
0x699e2ed19546484ULL, 0x3ddd5502c523b38eULL, 0x9c238a3d4699213ULL, 0x15c2d4f3288f1ac2ULL,
0x794b9fde5539fdbdULL, 0xa4c8edc5122fd8ULL, 0x584337f0f710607eULL, 0x4fbdf8fe52bbf0c6ULL,
0xbeed9782dc8484ceULL, 0x3c472eec38398e88ULL, 0xfe143f8c6216583cULL, 0xf457a23a19e8513ULL,
0x41ab2c4dc2f58383ULL, 0x4e30b2af359c2c7ULL, 0x701e766dcafd105aULL, 0x3b0801a0c3bde485ULL,
0x85c8da6e017b931fULL, 0x5232f2ae37f5de0cULL, 0x67d90613db73a5bfULL, 0x6b847a4b2317c98cULL,
0x82feac69a7cd7015ULL, 0x424b79c5585cfe26ULL, 0x7c1693460b2d0634ULL, 0x70121a6c9a7816f8ULL,
0x4b1a9e9394677651ULL, 0x5bbc81731836bce1ULL, 0xbfd55a98f807269aULL, 0x5b0dcf6c50e5616ULL,
0x3689355edaa6eefcULL, 0x433575dc39d6f3eeULL, 0x53c376fb71e6d0faULL, 0x7b75ab4dd7460312ULL,
0x5b4ae4b231bf04b8ULL, 0x5a9677f87b599ba8ULL, 0xe0326bc4fa2ca744ULL, 0x11ad8c06d5340001ULL,
0x57ce96720ace4297ULL, 0x4c1c7bad23245ce3ULL, 0x50d85ed4f867cf41ULL, 0x57d3671eab5ad319ULL,
0x7223dc981dd6b0b2ULL, 0x5d5892bdd9550776ULL, 0x74ac25392e70a059ULL, 0x2358c235433a1997ULL,
0x977add80da922cc8ULL, 0x613cedca577712a4ULL, 0x9dc42052d7a63aebULL, 0x46e7022975158bb9ULL,
0x51fd01eb7fc05d32ULL, 0x319921802c6ce600ULL, 0x8578eb315d1fecd3ULL, 0x4c0611a9e939520cULL,
0x3c19db5d8fea9e50ULL, 0x26a00f12f03a8070ULL, 0x4a0fa0f0d75dbb23ULL, 0x236528ed6438dcb8ULL,
0xf64c4a1fa82ae5feULL, 0x24bf9649017e3859ULL, 0xd2efefb6e3885fe0ULL, 0x430aaf1a36537cf9ULL,
0xefd92a7772283727ULL, 0x12f495f5c9c691dULL, 0x5b772f67c7ff429fULL, 0x33d5c5e60a3f64eeULL,
0x284d22038b8acb6fULL, 0x4f69f38145ba86deULL, 0x306eb189c67d6d75ULL, 0x4f0f79eb2d924523ULL,
0x98505c85c16c4018ULL, 0x684bc14ae520d034ULL, 0xd9bc8dbc651e3e48ULL, 0x62bbd623a2532c35ULL,
0x9e5d7a332ac5a2e7ULL, 0x134ac10a66909f2bULL, 0x294ded79913c9c5cULL, 0x5a4bfcf34aae4228ULL,
0x652d44970cf7dc69ULL, 0x35f08a571316b46eULL, 0x3679c35cb7fd0d19ULL, 0x1c1a844974083cULL,
0xd3636c17b3670754ULL, 0x2365a8e51418582eULL, 0xa49012aa91ef513cULL, 0x699a784c2104cb09ULL,
0x4f3fb4ddc036cb3eULL, 0x23d67a308148fa0dULL, 0xd9aa34fc6325ed26ULL, 0x3da149343f59bacdULL,
0x4dc557d1ceade366ULL, 0x6c9bb47f8dc8c6e2ULL, 0x2966e21030419d59ULL, 0x25ad1874cf53330cULL,
0xf8193baab903d60ULL, 0x401239915da48fbULL, 0x27584c2ee8059851ULL, 0x6acdc10404b6cdb0ULL,
0x9144fb3d23df8e3bULL, 0x4964a5a54b5cbf87ULL, 0xe8b7ae5657c57309ULL, 0x798dfb0a66eed55ULL,
0x21e36fef9aa453eeULL, 0x1d9e4ffcc4818d89ULL, 0x502b1d5a349f40dbULL, 0xdaecb8641228475ULL,
0xafad7aae90b62197ULL, 0x1ed18919e6d99efdULL, 0x1796a98ff289a29dULL, 0x28fa365c0654734dULL,
0x768f1d887cc53bbcULL, 0x6039a4e94561d1a9ULL, 0xc33b84f782ff1eb0ULL, 0xdbbaad1a410de8ULL,
0x26d098dd9a82e9e2ULL, 0x65ad0678ec7b3f4aULL, 0xfc1f64cb51483f1dULL, 0x5f2583bdec75110aULL,
0xbe71cd89f6f08295ULL, 0x1c1920166cb08e9ULL, 0xeb47a27907ee65d5ULL, 0x2f897a5c90e5b195ULL,
0x139a89fc039e042bULL, 0x22ab1e0ad85f2880ULL, 0x148f5c77b2867979ULL, 0x5876f0dd99a78b21ULL,
0xc1c67a449ef66098ULL, 0x5745671ca2f8ad5cULL, 0x10aee4c9a164bfbcULL, 0x4ed05a9445307402ULL,
0x4ddc882a167dfd2aULL, 0x55462eec69343ebbULL, 0x67f019200d101b1aULL, 0x38ab946a070d294fULL,
0xd6b80326b09e9d68ULL, 0xd30138eb1371936ULL, 0xa2d69161c3f5404cULL, 0x17e3a61335c82996ULL,
0x5ec96947147c5b40ULL, 0x6f7596f84153bf3bULL, 0xd106f0a4d8513d9eULL, 0x7062da7f4d727f6aULL,
0x16940356986e5c24ULL, 0x2a1b40cec77569bdULL, 0xffb69837b5452fb2ULL, 0x369e611333110a91ULL,
0xa2909146ce928867ULL, 0x6c16ee795a75f70dULL, 0xf821924710665bacULL, 0x2103ecb59c2e3230ULL,
0x46422886893d69b8ULL, 0x30d0a84850472264ULL, 0xd54601c8a2249625ULL, 0x35a80aa88f045fb9ULL,
0x54146941bb47d51aULL, 0x1743514dd74a38aeULL, 0x7fe16754f27d2b65ULL, 0x5da646f96dd67edULL,
0x76cbd638a948ddceULL, 0x6d90356cefd71933ULL, 0xb90f9d9d4bbd6c1aULL, 0x32a0d63fa9b9d1d9ULL,
0x45cc2bc2a672987aULL, 0x2b6ce4a8ea94fe8bULL, 0x6c040a5623a45536ULL, 0x4daeedfba3935f53ULL,
0xc011a713124d456ULL, 0x6dcdfa3248fff9bdULL, 0x27724eedf4498ad0ULL, 0x186204c65e6fb75eULL,
0xfc3cbcd671904180ULL, 0x4a0b9e4a30d4a61fULL, 0x78b18826171fa37dULL, 0x4526c2a804f96d1fULL,
0xca318c731a9120adULL, 0x60c00a0d6924cd58ULL, 0xd1ccf5a0715a6e78ULL, 0x2c2933debb738cebULL,
0x48c3d9a1a1c75481ULL, 0x655e13993d052aaaULL, 0xb1cd31cdf0d8637aULL, 0x86930936d8a3528ULL,
0x1da4274ed80ffb84ULL, 0x507c19ef54a22771ULL, 0xec4785805dd70009ULL, 0x8c54674896801c6ULL,
0x7976d5fb96bba9abULL, 0x6fc6b18cb344b902ULL, 0x1f061e28740d854bULL, 0x4ff326bf97e2179eULL,
0x5e948d170a67da6ULL, 0x2fd78e41cc1700abULL, 0x2103d23e7631f84ULL, 0x63ad8c716289e1a2ULL,
0x204e67963fa7c353ULL, 0x379cdd8bc74c0ed8ULL, 0xb10b05c4f273f986ULL, 0x328a498c5a445a4aULL,
0xcf6c0773cacea123ULL, 0x51f5b579e609968ULL, 0x618a0564c515421dULL, 0x66c0a1186378cc88ULL,
0xc9ea5316922d7a3aULL, 0x138c210fbf19494bULL, 0xf534c206a03b002bULL, 0x36c1173d05b1b257ULL,
0x6dc760e563b366cfULL, 0x7dcd442d281452e2ULL, 0x4d36631606e0264cULL, 0x476dc66db8c8f2d1ULL,
0x89eaf4801ff7dd54ULL, 0x1e9fdeca6657c319ULL, 0xa1a466ecb2c3a8f4ULL, 0x4d94d7e7fd22bbc6ULL,
0x6bf5b58257ed387bULL, 0x4c836f8e006da76cULL, 0x228098496beca652ULL, 0x1f6f3d91d7e33eb3ULL,
0x828edc57329a0b1fULL, 0x78332cb8e4cc90a1ULL, 0xe471963f58a8482eULL, 0x26109ff020788830ULL,
0xcc47e645d3e6dd60ULL, 0x4328bac901799627ULL, 0xbd4de219e202b90fULL, 0x2842ef57ff1e7b3fULL,
0x4367d81c26e753f0ULL, 0x42a23a2f38af711cULL, 0x2ef5ec8d2b20c2a1ULL, 0x46dad53966b3027fULL,
0x67742f32ef02195eULL, 0x2ad4a601e1cad8ddULL, 0xf4e138571d573e33ULL, 0x34dded3f3e112262ULL,
0x29dc4ae6adbb6b5aULL, 0x4dffdcd5c5d3d7d3ULL, 0xd606653fbb258010ULL, 0x565921e9d6dd8b46ULL,
0x8abde2d1c877aae7ULL, 0x731fd6ee427f3440ULL, 0x403a81fb168b1a38ULL, 0x5ac4a26c7d9a840ULL,
0x94be06c4746f4d26ULL, 0x7e557cb4b4be9065ULL, 0xe1e5fe953f658cadULL, 0x6b03a9a02bc91ccdULL,
0xd80e81d6354d26cULL, 0x77a9364f604a3855ULL, 0xa2d15ba6481b65e2ULL, 0x6b48ca5461de5450ULL,
0xdf32838fba98620ULL, 0x12dbbcc92d04bbd5ULL, 0x462ec81a8bdc56b4ULL, 0x7fb9b2d0aec54030ULL,
0xdb99db0336ce3656ULL, 0x4dbf964df23aedcbULL, 0xf182e5fabd433280ULL, 0x2b22fd68d32c726bULL,
0x974141cced5f91dfULL, 0x15104380fe0103baULL, 0xace0a79c61c8b0bfULL, 0x38e42addd79f7573ULL,
0x4a927d03f7eeda77ULL, 0x5390589168f06c23ULL, 0x6b5dc677531c6abeULL, 0x396d5ac5e6cce75dULL,
0x6b3c947deb9a5dc0ULL, 0xeba709dc9a12a89ULL, 0xb0f25aac1fabf1b7ULL, 0x4fdf324318914cULL
};

static const ecpt_affine *SIMUL_GEN_TABLE = (const ecpt_affine *)PRECOMP_TABLE_3;

static u32 ec_recode_scalar_comb1(const u64 k[4], u64 b[4]) {
	//const int t = 252;
	const int w = 8;
	const int v = 1;
	const int e = 32; // t / wv
	const int d = e * v;
	const int l = d * w;

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

static CAT_INLINE u32 comb_bit1(const u64 b[4], const int wp, const int ep) {
	// K(w', v', e') = b_(d * w' + e * v' + e')
	u32 jj = (wp << 5) + ep;

	return (u32)(b[jj >> 6] >> (jj & 63)) & 1;
}

// NOTE: Not constant time because it does not need to be for ec_simul_gen
void ec_table_select_comb1(const u32 recode_lsb, const u64 b[4], const int ii, ecpt &p1) {
	// DCK(v', e') = K(w-1, v', e') || K(w-2, v', e') || ... || K(1, v', e')
	// s(v', e') = K(0, v', e')

	// Select table entry 
	// p1 = s(0, ii) * tables[DCK(0, ii)][0]
	// p2 = s(1, ii) * tables[DCK(1, ii)][1]

	u32 d_0;
	d_0 = comb_bit1(b, 7, ii) << 6;
	d_0 |= comb_bit1(b, 6, ii) << 5;
	d_0 |= comb_bit1(b, 5, ii) << 4;
	d_0 |= comb_bit1(b, 4, ii) << 3;
	d_0 |= comb_bit1(b, 3, ii) << 2;
	d_0 |= comb_bit1(b, 2, ii) << 1;
	d_0 |= comb_bit1(b, 1, ii);
	u32 s_0 = comb_bit1(b, 0, ii);

	ec_expand(SIMUL_GEN_TABLE[d_0], p1);

	// Flip sign here rather than at the end to interleave more easily
	ec_cond_neg(s_0 ^ recode_lsb ^ 1, p1);
}

/*
 * Simultaneous multiplication by two base points,
 * where one is variable and the other is the generator point
 *
 * Preconditions:
 * 	0 < a,b < q
 *
 * Multiplies the result of aG + bQ by 4 and stores it in R
 */

// R = a*4*G + b*4*Q
void ec_simul_gen(const u64 a[4], const u64 b[4], const ecpt_affine &Q, ecpt_affine &R) {
	// Decompose scalar into subscalars
	ufp b1, b2;
	s32 b1sign, b2sign;
	gls_decompose(b, b1sign, b1, b2sign, b2);

	// Q2 = endomorphism of Q
	ecpt_affine Qa;
	gls_morph(Q.x, Q.y, Qa.x, Qa.y);
	ecpt Q2;
	ec_expand(Qa, Q2);
	ec_cond_neg(b2sign, Q2);

	// Set base point signs
	ecpt Q1;
	ec_expand(Q, Q1);
	ec_cond_neg(b1sign, Q1);

	// Precompute multiplication table
	ecpt qtable[8];
	ec_gen_table_2(Q1, Q2, qtable);

	// Recode subscalars
	u64 a1[4];
	const u32 comb_lsb = ec_recode_scalar_comb1(a, a1);
	u32 recode_bit = ec_recode_scalars_2(b1, b2, 128);

	// Initialize working point
	ecpt X;
	ec_table_select_2(qtable, b1, b2, 126, false, X);

	ufe t2b;
	ecpt T;
	for (int ii = 124; ii >= 32; ii -= 2) {
		ec_table_select_2(qtable, b1, b2, ii, false, T);

		ec_dbl(X, X, false, t2b);
		ec_dbl(X, X, false, t2b);
		ec_add(X, T, X, false, false, false, t2b);
	}

	// For the last 32 doubles, interleave ec_mul_gen adds

	for (int ii = 30; ii >= 0; ii -= 2) {
		ec_dbl(X, X, false, t2b);

		ec_table_select_comb1(comb_lsb, a1, ii+1, T);
		ec_add(X, T, X, true, false, false, t2b);

		ec_dbl(X, X, false, t2b);
		ec_table_select_comb1(comb_lsb, a1, ii, T);

		ec_add(X, T, X, true, false, false, t2b);

		ec_table_select_2(qtable, b1, b2, ii, false, T);
		ec_add(X, T, X, false, false, false, t2b);
	}

	// If bit == 1, X <- X + Q1 (inverted logic from [1])
	if (recode_bit != 0) {
		ec_add(X, Q1, X, true, false, false, t2b);
	}

	// Multiply by 4 to avoid small subgroup attack
	ec_dbl(X, X, false, t2b);
	ec_dbl(X, X, false, t2b);

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

