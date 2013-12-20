## Snowshoe
### Portable, Secure, Fast Elliptic Curve Math Library in C

[This project](http://github.com/catid/snowshoe) aims to provide a simple C API for various types of optimized elliptic curve point multiplication:

+ Fixed-base (Public key generation, Signature generation) `mul_gen`
+ Variable-base (Diffie-Hellman key agreement) `mul`
+ Variable double-base Simultaneous (EC-DH-FS key agreement) `simul`

Each multiplication routine is fast, constant-time, simple, easy to analyze,
portable, well-documented, and uses no dynamic memory allocation.

It is designed for a "128-bit" security level to be used with 256-bit keys.

On side-channel attack resilience: Passive timing attacks are defeated
by regular memory access and code execution patterns.  Power analysis attacks
are mitigated by these methods, but no blinding is in place to avoid leaking
information in power traces.

Additionally to speed up signature verification a variable single-base
simultaneous function `simul_gen` is provided that is not constant-time.

It is intended to be a reliable and robust library that provides the fastest
low-complexity, open-source implementation of these math routines available,
which should be considered wherever you would use e.g. Curve25519/Ed25519.

SUPERCOP Level 0 copyright/patent protection: There are no known present or future claims by a copyright/patent holder that the distribution of this software infringes the copyright/patent. In particular, the author of the software is not making such claims and does not intend to make such claims.


## Benchmarks

On Sandy Bridge and newer processors, the RDTSC instruction does not return
actual clock cycle counts.  Instead it returns a clock at a fixed frequency,
which is determined during testing.  The actual processor cycles are typically
at a higher Turbo Boost frequency under load, which varies unpredictably.
The median of 10,000 measurements is taken, and repeat measurements indicate
that they are accurate to ~100 cycles.  I was able to demonstrate that the
Sandy Bridge cycle counts on my laptop and desktop matched with TB off,
which validates the benchmarking methodology.

+ To disable TB on Mac, use [DisableTurboBoost](https://github.com/nanoant/DisableTurboBoost.kext) (included under ./tests)
+ To disable TB on Windows, edit the power settings to peg the processor at 100% frequency.  You may verify this with CPUz.


##### libsnowshoe.a on Macbook Air (1.7 GHz Core i7-4650U Haswell, June 2013):

+ RDTSC instruction runs at 2.29988 GHz so a correction factor of 1.7/2.29988 = 0.7391690 is needed.

`make snowshoetest` results (TB off) after applying the correction factor:

+ EC-DH client: `123748` corrected cycles, `73.2569` avg usec
+ EC-DH server: `123751` corrected cycles, `73.2452` avg usec
+ EdDSA sign: `61304` corrected cycles, `-134.132` avg usec
+ EdDSA verify: `137672` corrected cycles, `82.6574` avg usec

These corrected cycle counts can be directly compared with the latest SUPERCOP
Haswell results.  Not sure why time went backwards for signing. :-/

`make ecmultest` results (TB on, demonstrating usual walltime):

+ EC-DH client: `86284` median cycles, `37.9291` avg usec
+ EC-DH server: `86286` median cycles, `37.9505` avg usec
+ EdDSA sign: `41894` median cycles, `18.5784` avg usec
+ EdDSA verify: `95956` median cycles, `42.4539` avg usec


##### libsnowshoe.a on Macbook Air (1.7 GHz Core i5-2557M Sandy Bridge, July 2011):

RDTSC instruction runs at 1.70004 GHz so no correction factor is needed.

`make ecmultest` results (TB off):

+ ec_mul: `136788` median cycles, `81.2626` avg usec
+ ec_mul_gen: `62000` median cycles, `36.7403` avg usec
+ ec_simul: `199796` median cycles, `118.398` avg usec
+ ec_simul_gen: `159172` median cycles, `94.3428` avg usec

`make snowshoetest` results (TB off):

+ EC-DH client: `137528` median cycles, `81.1455` avg usec
+ EC-DH server: `137316` median cycles, `81.0305` avg usec
+ EC-DH-FS client gen: `73756` median cycles, `43.8528` avg usec
+ EC-DH-FS server proc: `139072` median cycles, `82.5851` avg usec
+ EC-DH-FS client proc: `203432` median cycles, `120.562` avg usec
+ EdDSA sign: `63488` median cycles, `37.6539` avg usec
+ EdDSA verify: `161892` median cycles, `95.7937` avg usec

`make ecmultest` results (TB on, demonstrating usual walltime):

+ ec_mul: `86216` median cycles, `51.1041` avg usec
+ ec_mul_gen: `38982` median cycles, `23.1696` avg usec
+ ec_simul: `124536` median cycles, `73.8478` avg usec
+ ec_simul_gen: `99698` median cycles, `59.3` avg usec

`make snowshoetest` results (TB on, demonstrating usual walltime):

+ EC-DH client: `86212` median cycles, `51.1257` avg usec
+ EC-DH server: `86126` median cycles, `51.0684` avg usec
+ EC-DH-FS client gen: `46312` median cycles, `27.4853` avg usec
+ EC-DH-FS server proc: `86808` median cycles, `51.4173` avg usec
+ EC-DH-FS client proc: `126934` median cycles, `75.1333` avg usec
+ EdDSA sign: `40408` median cycles, `23.9925` avg usec
+ EdDSA verify: `100882` median cycles, `59.7367` avg usec


##### libsnowshoe.a on iMac (2.7 GHz Core i5-2500S Sandy Bridge, June 2011):

RDTSC instruction runs at 2.69393 GHz so no correction factor is needed.

`make ecmultest` results (TB off):

+ ec_mul: `137448` median cycles, `51.1046` avg usec
+ ec_mul_gen: `62732` median cycles, `23.3941` avg usec
+ ec_simul: `199656` median cycles, `74.1709` avg usec
+ ec_simul_gen: `158080` median cycles, `58.8016` avg usec

These results serve to justify the benchmarking methodology.  Notice that the
cycle counts exactly match the laptop version.

`make snowshoetest` results (TB off):

+ EC-DH client: `136608` median cycles, `50.8648` avg usec
+ EC-DH server: `136428` median cycles, `50.7681` avg usec
+ EC-DH-FS client gen: `74584` median cycles, `27.7514` avg usec
+ EC-DH-FS server proc: `138048` median cycles, `51.3688` avg usec
+ EC-DH-FS client proc: `201932` median cycles, `75.1432` avg usec
+ EdDSA sign: `63672` median cycles, `23.7501` avg usec
+ EdDSA verify: `161520` median cycles, `60.083` avg usec

`make ecmultest` results (TB on, demonstrating usual walltime):

+ ec_mul: `100484` median cycles, `38.2353` avg usec
+ ec_mul_gen: `46533` median cycles, `17.9244` avg usec
+ ec_simul: `145766` median cycles, `55.5407` avg usec
+ ec_simul_gen: `115574` median cycles, `44.0471` avg usec

`make snowshoetest` results (TB on, demonstrating usual walltime):

+ EC-DH client: `99710` median cycles, `38.1784` avg usec
+ EC-DH server: `99582` median cycles, `38.1185` avg usec
+ EC-DH-FS client gen: `54438` median cycles, `20.964` avg usec
+ EC-DH-FS server proc: `100756` median cycles, `38.7462` avg usec
+ EC-DH-FS client proc: `147472` median cycles, `56.6324` avg usec
+ EdDSA sign: `46838` median cycles, `17.9698` avg usec
+ EdDSA verify: `118128` median cycles, `45.1675` avg usec


##### libsnowshoe.lib on Windows 7 laptop (2.67 GHz Core i7 620M Westmere, Jan 2010):

RDTSC instruction runs at 2.66008 GHz, CPU pegged at 3.3254 GHz implies a correction factor of 3.3254/2.66008 = 1.2501128

MSVC2010 build results (CPU at 100%) before applying the correction factor:

+ EC-DH client: `152939` median cycles, `61.2776` avg usec
+ EC-DH server: `152744` median cycles, `61.2637` avg usec
+ EC-DH-FS client gen: `81463` median cycles, `32.2964` avg usec
+ EC-DH-FS server proc: `154260` median cycles, `60.9112` avg usec
+ EC-DH-FS client proc: `233256` median cycles, `88.8928` avg usec
+ EdDSA sign: `64234` median cycles, `25.5279` avg usec
+ EdDSA verify: `183615` median cycles, `72.5732` avg usec

MSVC2010 build results (CPU at 100%) after applying the correction factor:

+ EC-DH client: `191191` median cycles, `61.2776` avg usec
+ EC-DH server: `190947` median cycles, `61.2637` avg usec
+ EC-DH-FS client gen: `101838` median cycles, `32.2964` avg usec
+ EC-DH-FS server proc: `192842` median cycles, `60.9112` avg usec
+ EC-DH-FS client proc: `291596` median cycles, `88.8928` avg usec
+ EdDSA sign: `80300` median cycles, `25.5279` avg usec
+ EdDSA verify: `229539` median cycles, `72.5732` avg usec


#### Usage

This git repo uses submodules so be sure to run `git submodule update --init` to download all the code.

You can either compile-in the software or link to it.  I recommend statically linking
the code, since that enables full optimization and speeds up your compilation.

To build the project you only need to compile `src/snowshoe.cpp`, which includes
all of the other source files.  Or link to a prebuilt static library under `bin/`

To use the project you only need to include [include/snowshoe.h](https://github.com/catid/snowshoe/blob/master/include/snowshoe.h), which declares the C exports from the source files.

An example project that uses Snowshoe for signatures and handshakes is [Tabby](https://github.com/catid/tabby).


#### Building: Mac

To build the static library, install command-line Xcode tools and simply run the make script:

~~~
make release
~~~

This produces `libsnowshoe.a` with optimizations.


#### Building: Windows

Download LLVM from [http://llvm.org/builds/](http://llvm.org/builds/) for Windows to C:\LLVM\.
Download Mingw64 from [http://mingw-w64.sourceforge.net/](http://mingw-w64.sourceforge.net/) for Windows 64-bit to C:\mingw64\.

~~~
copy Makefile.mingw64 Makefile
c:\mingw64\bin\mingw32-make.exe release
~~~

This produces `libsnowshoe.lib` with vector extensions, which can be linked to an MSVC2010 build.


#### Example Usage: EC-DH

[Elliptic Curve Diffie-Hellman](http://en.wikipedia.org/wiki/Elliptic_curve_Diffie%E2%80%93Hellman) is the baseline for key agreement over the Internet.  Implementing it with this library is straight-forward:

Verify binary API compatibility on startup:

~~~
	if (snowshoe_init()) {
		throw "Buildtime failure: Wrong snowshoe static library";
	}
~~~

Allocate memory for the keys:

~~~
	char sk_c[32], sk_s[32];
~~~

Fill `sk_c` and `sk_s` with random bytes here.  Snowshoe does not provide a random number generator.

Now generate the server public/private key pair:

~~~
	char pp_s[64];
	snowshoe_secret_gen(sk_s);
	if (snowshoe_mul_gen(sk_s, pp_s)) {
		throw "Secret key was generated wrong (developer error)";
	}
~~~

`snowshoe_secret_gen` will mask off some bits of the random input string to make it suitable for use as a private key.

`snowshoe_mul_gen` takes a prepared private key and multiplies it by the base point (described below) to produce public point (X,Y) coordinates encoded as 64 bytes.

Generate client public/private key pair:

~~~
	char pp_c[64];
	snowshoe_secret_gen(sk_c);
	if (snowshoe_mul_gen(sk_c, pp_c)) {
		throw "Secret key was generated wrong (developer error)";
	}
~~~

Client side: Multiply client secret key by server public point

~~~
	char sp_c[64];
	if (snowshoe_mul(sk_c, pp_s, sp_c)) {
		// Reject server input pp_s that is invalid
		return false;
	}
~~~

Server side: Multiply server secret key by client public point

~~~
	char sp_s[64];
	if (snowshoe_mul(sk_s, pp_c, sp_s)) {
		// Reject client input pp_c that is invalid
		return false;
	}
~~~

Server and client both arrive at `sp_c == sp_s`, which is the secret key for the session.

The error checking used above should be replaced with more suitable reactions to failures in production code.  It is important to check if the functions return non-zero for failure, since this indicates that the other party has provided bad input in an attempt to attack the cryptosystem.


#### Example Usage: EC-DH-FS

An improvement on EC-DH is the concept of "forward secrecy", where if the server's long-term
secret key is revealed it is not any easier to decrypt past logged communication.

This example sketches how to implement a protocol that can achieve this goal.  It requires
that the client generate a new public key periodically to delete the past sessions, and
the server should also regenerate an ephemeral key periodically.

After the key agreement completes, the server should provide a proof that it knows the
secret key in its response to the client.  This message can also carry the first encrypted
packet.

~~~
	char h[32], d[32], a[32];
	char sk_c[32], sk_s[32], sk_e[32];
	char pp_c[64], pp_s[64], pp_e[64];
	char sp_c[64], sp_s[64];

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

	generate_k(sk_c);
	snowshoe_secret_gen(sk_c);
	if (snowshoe_mul_gen(sk_c, pp_c)) {
		return false;
	}

	// h = H(pp_s, pp_e, pp_c, client_nonce, server_nonce)
	generate_k(h);

	// Online: Server handles client request

	// d = sk_e + h * sk_s (mod q)
	snowshoe_mul_mod_q(h, sk_s, sk_e, d);
	if (snowshoe_mul(d, pp_c, sp_s)) {
		return false;
	}

	// Online: Client handles server response

	// a = h * sk_c (mod q)
	snowshoe_mul_mod_q(h, sk_c, 0, a);
	if (snowshoe_simul(sk_c, pp_e, a, pp_s, sp_c)) {
		return false;
	}

	for (int ii = 0; ii < 64; ++ii) {
		if (sp_c[ii] != sp_s[ii]) {
			return false;
		}
	}
~~~

The error checking used above should be replaced with more suitable reactions to failures in production code.  It is important to check if the functions return non-zero for failure, since this indicates that the other party has provided bad input in an attempt to attack the cryptosystem.


#### Example Usage: EdDSA

Here is a sketch of how to implement [EdDSA](http://ed25519.cr.yp.to/ed25519-20110926.pdf):

Note that the "key massaging" would be different for my group order.
Instead of Ed25519 key masking, use `snowshoe_secret_gen`.

Key generation:

+ Generate a random number k < 2^256
+ hi,lo = H(k)
+ a = snowshoe_secret_gen(lo)
+ A = a*G

Sign message M:

+ r = H(hi,M) (mod q)
+ t = H(R,A,M) (mod q)
+ R = r*G
+ s = r + t*a (mod q)
+ Produce: R, s

Verify:

+ u = H(R,A,M) (mod q)
+ nA = -A
+ R =?= s*G + u*nA

~~~
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
	if (snowshoe_mul_gen(a, pp_A)) {
		return false;
	}

	// Sign:

	snowshoe_mod_q(h_hi_m, r);
	snowshoe_mod_q(h_r_a_m, t);
	if (snowshoe_mul_gen(r, pp_R)) {
		return false;
	}
	snowshoe_mul_mod_q(a, t, r, s); // s = a * t + r (mod q)

	// Verify:

	snowshoe_mod_q(h_r_a_m, u);
	snowshoe_neg(pp_A, pp_A);
	if (snowshoe_simul_gen(s, u, pp_A, pp_Rtest)) {
		return false;
	}

	// Test if 4 * pp_R == pp_Rtest
	if (snowshoe_equals4(pp_Rtest, pp_R)) {
		return false;
	}
~~~

Note that this library provides a number of helpful math functions for doing math modulo q.


#### Comparisons

This section demonstrates performance comparisons with other fast timing-attack resilient
elliptic curve math libraries at the roughly "128-bit" security level.  Note that only
Curve25519 and Ed25519 are fully in the public domain at this time.

##### Curve25519 [http://cr.yp.to/ecdh/curve25519-20060209.pdf](http://cr.yp.to/ecdh/curve25519-20060209.pdf):

- Availability : Equivalently free, open-source, and portable

On 2.9 GHz i5-4570S Haswell, TB off, using SUPERCOP:

- ecmul_gen : (slower) `171kcy`
- ecmul : (slower) `162kcy`
- ecsimul_gen : (not implemented)
- ecsimul : (not implemented)

##### Ed25519 [http://ed25519.cr.yp.to/ed25519-20110926.pdf](http://ed25519.cr.yp.to/ed25519-20110926.pdf):

- Availability : Equivalently free, open-source, and portable

On 2.9 GHz i5-4570S Haswell, TB off, using SUPERCOP:

- ecmul_gen : (similar) `67kcy`
- ecmul : (not implemented)
- ecsimul_gen : (slower) `207kcy`
- ecsimul : (not implemented)

##### monfp127e2 [http://eprint.iacr.org/2013/692.pdf](http://eprint.iacr.org/2013/692.pdf):

- Availability : Free, open-source, but not portable (uncommented assembly only)

On i7-3520M Ivy Bridge, TB off, using SUPERCOP:

- ecmul_gen : (not implemented)
- ecmul : (slower) `145kcy`
- ecsimul_gen : (not implemented)
- ecsimul : (not implemented)

##### kumfp127g [http://eprint.iacr.org/2012/670.pdf](http://eprint.iacr.org/2012/670.pdf):

- Availability : Free, open-source, but not portable (uncommented assembly only)
- This code is also extremely complex and looks tricky to audit.

On i7-3520M Ivy Bridge, TB off, using SUPERCOP:

- ecmul_gen : (slower) `108kcy`
- ecmul : (faster) `110kcy`
- ecsimul_gen : (not implemented)
- ecsimul : (not implemented)

##### Hamburg's implementation [http://mikehamburg.com/papers/fff/fff.pdf](http://mikehamburg.com/papers/fff/fff.pdf):

- Availability : Not available online

On i7 2720QM Sandy Bridge with TB on, but scaled cycle counts:

- ecmul_gen : (faster) `60kcy`
- ecmul : (similar) `153kcy`
- ecsimul_gen : (similar) `169kcy`
- ecsimul : (not implemented)

##### Longa's implementation [http://eprint.iacr.org/2013/158](http://eprint.iacr.org/2013/158):

- Availability : Not available online?

On 3.4 GHz i7-2600 Sandy Bridge with TB off:

- ecmul_gen : (faster) `48kcy`
- ecmul : (faster) `96kcy`
- ecsimul_gen : (faster) `116kcy`
- ecsimul : (not implemented)

On 3.4 GHz i7-3770 Ivy Bridge with TB off:

- ecmul_gen : (faster) `46kcy`
- ecmul : (faster) `92kcy`
- ecsimul_gen : (faster) `111kcy`
- ecsimul : (not implemented)

##### Crypto++ Library 5.6.2

On iMac (2.7 GHz Core i5-2500S Sandy Bridge, June 2011):

Using `cryptopp562/cryptest.exe`, the best 256-bit ECDSA implementation takes
`2,040,000 cycles` for signing, and `7,750,000 cycles` for verification.
Similarly all other operations are roughly 20x slower than Snowshoe.

##### LibTomCrypt Library 1.17

On iMac (2.7 GHz Core i5-2500S Sandy Bridge, June 2011):

After building and installing LibTomFastMath 0.12 and fixing some syntax errors
for Mac and disabling LTC_KSEED,
`CFLAGS="-DTFM_DESC -DUSE_TFM" EXTRALIBS=-ltfm make test timing` gave me a
binary, which passed all tests.

These results are roughly 20x slower than Snowshoe:

+ ECC-256 make_key    took         1738168 cycles
+ ECC-256 encrypt_key took         3464327 cycles
+ ECC-256 decrypt_key took         1785863 cycles
+ ECC-256 sign_hash took           1823073 cycles
+ ECC-256 verify_hash took         2436467 cycles

##### Kasper's OpenSSL implementation [http://static.googleusercontent.com/media/research.google.com/en/us/pubs/archive/37376.pdf](http://static.googleusercontent.com/media/research.google.com/en/us/pubs/archive/37376.pdf):

The NIST P-224 curve point multiplication Kasper presents in the paper ran in
`457,813` cycles for P-224 ("112 bit" security level), which is less secure
than this curve by about 14 "bits."  I do not have access to a Core 2 Duo to
do a direct comparison, but SUPERCOP has Curve25519 benchmarked at `306,936`
cycles on the same architecture, and for Sandy Bridge Curve25519 performs a
scalar multiplication in `205,696` cycles.  So I would expect Snowshoe to run
`ec_mul` in `204,190` cycles on the same machine.  Since Kasper does not
mention disabling Turbo Boost, Snowshoe is anywhere from 2.24x to 3x faster.

Kasper's results also make it possible to compare Snowshoe to OpenSSL.
In OpenSSL 1.0.1e, the NIST P-256 curve math in `crypto/ec/ecp_nistp256.c` is
similar in code complexity to Snowshoe.  However, OpenSSL does not take
advantage of efficient endormorphisms, fast twisted Edwards group laws, nor
fields with fast reductions or Karatsuba methods for field multiplication, so
it is much slower.  Kasper demonstrates that the improvements in the paper are
roughly 3x faster than normal OpenSSL, pegging Snowshoe roughly 6x faster than
OpenSSL 1.0.1e for the same security level.

And indeed building OpenSSL 1.0.1e with `make tests`, after configuring for
64-bit builds and enabling the timing tests macros in `crypto/ec/ectest.c` show
that it takes roughly 100 usec to generate a key (4x slower) and that it takes
roughly 400 usec for EC-DH (8x slower).

##### Nightcracker's Ed25519 implementation [https://github.com/nightcracker/ed25519](https://github.com/nightcracker/ed25519):

These are the Ed25519 results from the tester:

~~~
testing seed generation performance: 6483.760800us per seed
testing key generation performance: 107.181400us per keypair
testing sign performance: 103.336800us per signature
testing verify performance: 263.802000us per signature
testing keypair scalar addition performance: 93.632100us per keypair
testing public key scalar addition performance: 106.290500us per key
testing key exchange performance: 229.176800us per shared secret
~~~

Compared to Snowshoe on the same computer, this version of Ed25519 is 3.8x
slower for key generation, 4.33x slower for signing, 4.34x slower for
signature verification, and 4.36x slower for EC-DH.


## Details

~~~
"Simplicity is the ultimate sophistication."
                        -Leonardo da Vinci
~~~


#### Curve specification:

+ Field math: Fp^2 with p = 2^127-1
+ Curve shape: E' : a * u * x^2 + y^2 = 1 + d * u * x^2 * y^2
+ u = 2 + i
+ a = -1
+ d = 109
+ Group size: #E' = 4 * q,
+ q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA6261414C0DC87D3CE9B68E3B09E01A5


#### Compatible with SafeCurve guidelines [9]:

+ Base point and its endomorphism are of order q
+ Resists Pollard rho method with a large prime group order q
+ Not anomalous and embedding degree is large
+ CM field discriminant is large
+ Fully rigid curve design: All parameters are chosen for simplicity
+ The curve supports constant-time multiplication with Montgomery Ladders
+ Prevents small-subgroup/invalid-curve attacks by validating input points
+ The doubling and unified addition formulae used are complete
+ Supports Elligator construction, though not implemented here


#### Implementation features:

+ Uses constant-time arithmetic, group laws, and point multiplication
+ Critical math routines have proofs of correctness
+ All math routines are unit-tested for expected edge cases


#### Performance features:

+ Most efficient extension field arithmetic: Fp^2 with p = 2^127-1 [3]
+ Most efficient point group laws: Extended Twisted Edwards [3]
+ Efficient 2-dimensional GLS endomorphism [12]
+ Most efficient constant-time scalar multiplication: Windowed GLV-SAC [1]


## Choosing a Finite Field

## Choosing an Optimal Extension Field

This was easy.  The only real option is Fp^2.  Larger exponents provide
diminishing returns in terms of real security, and math over Fp^2 is the same
as grade-school complex math (a + b*i).

Multiplication over Fp^2 takes about 66 cycles on a Sandy Bridge processor,
whereas multiplication on well-chosen Fp of the same size takes 55 cycles. [17]
This performance loss is more than compensated for by faster field inversions
and allowing for efficient endomorphisms.


#### Why p = 2^127-1?

I care mainly about server performance, with 64-bit Linux VPS in mind, and
the Intel x86-64 instruction set has a fast 64x64->128 multiply instruction
that is exploited to simplify the code and speed up math on this field.

There are other fields with fast reductions, but the Fp used here (2^127-1)
seems to be the best option in terms of practical performance after reading
a lot of whitepapers.

Some useful details about Fp:

~~~
p = 1 (mod 3), p = 2 (mod 5), p = 7 (mod 8).

-1 = p - 1 is a non-square in Fp.  3 is also non-square in Fp.
~~~

#### Another Option: 2^31 - 1

Good for ARM systems, but extension fields are not ideal:

+ Fp^7 -> 217-bit keys, too far from targetted security level
+ Weil Descent attacks apply to Fp^8/9

#### Another Option: 2^61 - 1

Requires a 4-GLV decomposition to use properly, and I think 2 dimensional
decomposition is complicated enough for now.  Also it would take a huge
security hit to use.  Maybe not a huge deal.


## Choosing Endomorphism Dimensions

Trying to go for a 4-dimensional GLV-GLS combination seems like it adds more
complexity than I am comfortable with, so I decided to go for a 2-dimensional
decomposition to be more conservative despite good examples like [1].

There are effectively two choices for 2-dimensional endomorphisms at the moment:

+ Q-Curves: Practically implemented late 2013 [16]
+ GLS Curves: Known and analyzed since 2008 [12]


## Choosing X-Only or (X, Y) Coordinates

At this time, the best group laws are for Montgomery curves with X-Only
coordinates and extended twisted Edwards curves for (X, Y) coordinates.
Wisdom has dictated in the past that X-Only is better for performance.

However for 2-dimensional approaches this is no longer the case:

According to the group laws presented in [14] and the results
from [16], a `2`-dimensional Montgomery scalar multiplication performs
`1` add and `1` double-add per bit.  Yet it is shown in [1] that twisted
Edwards group laws allow for a protected table-based approach that requires
`1` double and `0.5` adds per bit with a practical window size of 2 bits.

With our finite field, squares require 0.67x the time of a multiply, so the
operation count is roughly `13.36` multiplies per bit for Montgomery curves,
and `(5.68 + 5.68 + 9) / 2 = 10.18` multiplies per bit for Edwards curves.
Note that the estimate of `9` multiplies above is pessimistic.

This demonstrates that (X, Y) coordinates are always preferred for speed for
all types of applications.  Furthermore (X, Y) coordinates are much more
flexible, as they are useful for signature verification.

The size difference between 32 byte public keys (X-Only) and 64 byte public
keys (X,Y) is minimal in practice.  RSA public keys are much larger and
are still used for real applications, for example.


## Choosing a Curve Shape

There was a choice between curves that have efficient endomorphisms:

GLV curves are rare and secure varieties cannot be found over our field [10].

Q-curve endomorphisms [10] are useful in the context of scalar multiplication
but thankfully not in the context of Pollard rho [4], so using these
endomorphisms does not reduce the effective security of the cryptosystem,
whereas GLS curve endomorphisms necessarily reduce security a little.

These Q-curves require that the twisted Edwards form has a large value for
the 'a' coefficient, whereas we would prefer a = -1, since it is multiplied
in both the doubling and the addition point operations.

GLS curves have a similar disadvantage in that curve operations are done
over the twist of an Fp curve, and so the twisted Edwards constant a != -1,
it is effectively `-(2+i)`, which requires 4 Fp additions per point operation.

Therefore I am using GLS curves despite the slight drop in security as is
done in [3].

The possibility of using a change of variables and specially chosen Q-curve
parameters to use the efficient dedicated addition formula from [5] for a = -1
with a correction post-step (LaineyCurves.md) was also explored.  However,
introducing additional complexity and dangerous incomplete addition laws to
save one multiplication per addition seems like a desperate choice to me for
little practical gain, especially since the number of additions is already
halved with a 2-bit window during evaluation.

Since multiplication by `u` appears in both the doubling and addition formulae,
it was chosen to be as efficient as possible to multiply: `u = (2 + i)`.
This is the same as presented in [3].  An alternative would be to use
`u = (1+2^64*i)` but this requires practically the same operation count while
adding more complexity.  Nicer options such as `(1+i)` are unfortunately square.
Sticking with the existing literature is best here.

And I chose a = -1.  a = 1 is also possible without making much difference,
except that I would have to diverge from existing literature without gaining
any real advantage.

The only remaining free parameter is d.  Using a MAGMA script
(magma_curvegen.txt) the following alternative curves were found:

~~~
u = 2 + i
a = -1

E: a * u * x^2 + y^2 = d * u * x^2 * y^2 + 1

d=109 : #E = 4*q, q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA6261414C0DC87D3CE9B68E3B09E01A5
d=139 : #E = 4*q, q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBE279B04A75463D09403332A27015D91
d=191 : #E = 4*q, q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF826EDB49112B894254575EA3A0C8BDC5
d=1345: #E = 4*q, q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF80490915366733181B4DC41442AAF491
d=1438: #E = 4*q, q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC96E66D7F2A4B799044761AE30653065
d=1799: #E = 4*q, q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF8FF32A5C1ACEC774E308CDB3636F2311
d=2076: #E = 4*q, q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF81EBFEA8A9E1FB42ED4A6EBB16B24A91
d=2172: #E = 4*q, q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF819920B3F8F71CD85DD3F4242C1B0E11
d=2303: #E = 4*q, q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF9B3E69111FF31FA521F8B59CC48B4101
d=2377: #E = 4*q, q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF94B9FB29B4A87B1DAEFA7A69FC19FD11
d=2433: #E = 4*q, q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF8F4C87E0F8EB73ABCB41D9C4CF92FC41
...
~~~

Again I decided to go with the smallest d = 109 to avoid diverging from the
literature without any real benefits.


## In-depth Curve Properties and Validation

This section explores all of the interesting properties of the curve, starting
with a formal definition of the curve.

If K is a finite field, then for all x there is a y coordinate for which
(x, y) belongs to either curve E or its twist E^d and sometimes both.

~~~
Curve E(Fp) : a * x^2 + y^2 = d * x^2 * y^2 + 1  (twisted Edwards form)
p = 2^127 - 1
a = -1
d = 109

Isomorphic to: y^2 = x^3 + A*x + B
A = 42535295865117307932921825928971026216
B = 85070591730234615865843651857942050915

r = #E(Fp) = 170141183460469231745929266840803953228
Factorization(r) = [ <2, 2>, <1789, 1>, <3041, 1>,
					<2427899077100477, 1>, <3220268376816859, 1> ]
t = TraceOfFrobenius(E) = -14241963124919847500 = p + 1 - r (verified)
~~~

~~~
Twist E'(Fp) is isomorphic to: y^2 = x^3 + A'*x + B'
A' = 112267898042735844839809317154290646844
B' = 120885228940773737117756558037765573846

r' = #E'(Fp) = 170141183460469231717445340590964258228
Factorization(r') = [ <2, 2>, <3, 1>, <11, 1>, <181, 1>, <443, 1>,
					<80447, 1>, <199822028697017221643157029, 1> ]
t' = TraceOfFrobenius(E') = 14241963124919847500 = p + 1 - r' (verified)
~~~


##### Formal Curve Definition

Now the cryptographically interesting twist of E(Fp):

~~~
Curve E(Fp^2) : a * u * x^2 + y^2 = d * u * x^2 * y^2 + 1  (twisted Edwards form)
u = 2 + i
i = sqrt(-D), D = 1 (reduces to usual complex i)

Isomorphic to: y^2 = x^3 + A*x + B
A = 170141183460469231731687303715884104864*i + 127605887595351923798765477786913078648
B = 85070591730234615865843651857942031430*i + 170141183460469231731687303715884101830

r = #E(Fp^2) = 28948022309329048855892746252171976962839764946219840790663900086538002237076
Factorization(r) = h * q (cofactor times large prime q)
h = 4
q = 7237005577332262213973186563042994240709941236554960197665975021634500559269
t = Tr(E(Fp)) = 14241963124919847500
Verified: r = (p - 1) ^ 2 + t^2
Verified: r * P = [0, 1] = point-at-infinity
~~~

~~~
Twist E'(Fp^2) is isomorphic to: y^2 = x^3 + A'*x + B'
A' = 74770288137151926641346873498907920755*i + 33812952767805664868511289460037871005
B' = 105685070695019342075176919984401662667*i + 99525886244248837335664845733467848745
r' = #E'(Fp^2) = 28948022309329048855892746252171976963114662652758564302138142702555026159984
(factorization omitted, but it is composite)
Verified: r' = (p + 1) ^ 2 - Tr(E'(Fp))^2
~~~


##### GLS Endormorphism

This curve has a useful endomorphism introduced in [12]:

The endomorphism is defined in Weierstrass coordinates as:

~~~
endomorphism(x, y) = (u * x^p / u^p, sqrt(u^3) * y^p / sqrt(u^(3p)))
u = 2 + i
~~~

To simplify, recall: a^p = conj(a) over Fp^2

~~~
endomorphism(x, y) = (wx * conj(x), wy * conj(y))
wx = u / u^p
   = 102084710076281539039012382229530463437*i + 34028236692093846346337460743176821146
wy = sqrt(u^3 / u^(3*p))
   = 81853545289593420067382843140645586077*i + 145299643018878500690299925732043646621

endomorphism(P) = lambda * P
~~~

Now to derive the value of lambda:

~~~
By [12], lambda = sqrt(-1) mod q = Modsqrt(q-1, q), q | #E(Fp^2)

NOTE: This differs from [12, section 2.1] in notation.  Here #E(Fp^2) is used in place
of #E'(Fp^2).  And `q` is used in place of `r`.

-> lambda = 6675262090232833354261459078081456826396694204445414604517147996175437985167
= 0xEC2108006820E1AB0A9480CCBB42BE2A827C49CDE94F5CCCBF95D17BD8CF58F
~~~

Using [6] I have verified that (p - 1) * InverseMod(t, q) mod q = lambda.

Using the attached script magma_endomorphism.txt I have verified that
`lambda * P = endomorphism(P)`.

Hurrah, the math works!

A formula for the endomorphism in twisted Edwards coordinates is provided in [3]:

~~~
endomorphism(x, y) = (wx' * conj(x), conj(y)) = lambda * (x, y)

wx' = sqrt(u^p/u)
= 68985359527028636873539608271459718931*i + 119563271493748934302613455993671912329
= 0x33E618D29DA66430D2B569B107BC1713 * i + 0x59F30C694ED33218695AB4D883DE0B89
~~~


##### Scalar Decomposition

Using the GLS method [12] Snowshoe will perform point multiplication by splitting the scalar to be
multiplied into two half-sized scalars, and then computing the sum of the half-sized products.

Lattice decomposition with precomputed basis vectors is the usual approach to decompose the scalar.
Suitable short basis vectors are provided in [15] and they were scaled to work modulo `q` rather than
`4*q`.

Given:

~~~
// Fp modulus=
p := 2^127-1;

// group order=
q := 7237005577332262213973186563042994240709941236554960197665975021634500559269;
r := 4 * q;

// group order trace=
t := 14241963124919847500;

// lambda=
l := 6675262090232833354261459078081456826396694204445414604517147996175437985167;

// b1 vector <x1, y1> = <(p-1)/2, -t/2>
// b2 vector <x2, y2> = <-t/2, (1-p)/2>
// Principal constants:
// A := (p - 1) div 2;
// B := t div 2;

A := 2^126 - 1; // Replace multiplication by A with shift and subtract
B := 0x62D2CF00A287A526; // Relatively small 64-bit constant
~~~

After substitution and arranging expressions to eliminate any negative results until the final
expressions for k1, k2:

~~~
a1 := A * k;
a2 := B * k;

//qround := q div 2 + 1;
qround := 0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD3130A0A606E43E9E74DB471D84F00D3;

z1 := (a1 + qround) div q;
z2 := (a2 + qround) div q;

// The final results may be negative after the subtractions:
k1 := k - (z1 * A + z2 * B);
k2 := z1 * B - z2 * A;

// k = (k1 + k2 * l) mod q
~~~

Note that multiplications by A can be achieved efficiently by a shift and subtract,
and that the constant B is a small 64-bit number - these are relatively cheap.

This decomposition guarantees that k1, k2 are 126 bits or less.  As shown in [15],
`||k_i|| <= ||(p-1)/2||, i={1,2}`, so their magnitude can be represented in 126 bits.
The division by 2 is due to scaling the basis vectors from [15] as suggested by that
paper to work over `q` modulus rather than `4*q` modulus.

The values k1, k2 can each be negative.  So a fast constant-time conditional point negation
must be performed on the base points before scalar multiplication.  An alternative approach
was taken in [16] to offset the lattice decomposition to avoid negative values, but it is
actually faster and easier to do conditional negation.

Since the decomposition works modulo `q` instead of modulo `4*q`, the base point(s) must be
premultiplied by 4 before scalar multiplication proceeds to avoid small-subgroup attack.
This is equivalent to clearing the low 2 bits of the scalar, except that it is more efficent
to reduce the sub-scalar length by 1 bit and premultiply by 2 doubles in the case of variable-
base scalar multiplication than it is to use the unscaled bases and perform an extra ec_add
during evaluation.

See magma_decompose.txt for more details.


##### Base Point

The base point should be in the q torsion group such that q * P is the identity
element (0, 1).

Since the curve has cofactor 4, there are actually 4 different types of base
points on the curve (in Weierstrass coordinates):

~~~
+ q * P = (18*i + 36, 0)
+ q * P = (68041245184837551457177097034226058488*i + 51049969137815840137255103340829023597, 0)
+ q * P = (102099938275631680274510206681658047221*i + 119091214322653391594432200375055082094, 0)
+ q * P = (0, 1)
~~~

And in twisted Edwards coordinates:

~~~
+ q * P = (0, 0)
+ q * P = (0, 1)
~~~

Note that for all of these, `4 * q * P = (0, 1)`.  And it is expected that the endomorphism of the
base point will be in the q torsion group by construction.

For rigidity [9], I chose the first twisted Edwards curve point of q torsion with a small, real X coordinate:

~~~
X = 0 * i + 15
Y = 0x6E848B46758BA443DD9869FE923191B0 * i + 0x7869C919DD649B4C36D073DADE2014AB
T = X * Y
Z = 1
~~~

Values of X from [0..14] are either not on the curve, or `q * (X, Y) = (0, 0)` instead of `(0, 1)`.

The endomorphism of the generator point in twisted Edwards coordinates:

~~~
X' = 0xA7B74573CBFDEDC58A1315F74055A23 * i + 0x453DBA2B9E5FEF6E2C5098AFBA02AD11
Y' = 0x117B74B98A745BBC226796016DCE6E4F * i + 0x7869C919DD649B4C36D073DADE2014AB
T' = X' * Y'
Z' = 1
~~~

Using the provided ecpt unit tester these points are validated to be on the curve and
of q torsion: `q * P = (0, 1)`.  See magma_generator.txt for full details.


##### Exceptional Input Scalars

Snowshoe's `ec_mul_gen` and `ec_simul_gen` generator scalars can take on the 0 value,
which will result in `k * G = (0, 1)`.  However the `ec_mul`, `ec_simul_gen` and
`ec_simul` input point scalars result in invalid output when the scalar is 0.

To address this problem, the input scalars for all of the Snowshoe multiplication
functions are checked to ensure they are in the expected range 0 < k < q.
After writing the [Tabby](https://github.com/catid/tabby) library, there does not
appear to be any reason to allow for 0 input.


## References

I am attempting to faithfully reproduce the results of some brilliant people:

##### [1] ["Keep Calm and Stay with One" (Hernandez Longa Sanchez 2013)](http://eprint.iacr.org/2013/158)
Introduces GLV-SAC exponent recoding

##### [2] ["Division by Invariant Integers using Multiplication" (Granlund Montgomery 1991)](http://pdf.aminer.org/000/542/596/division_by_invariant_integers_using_multiplication.pdf)
Division on fixed field in constant time

##### [3] ["Analysis of Efficient Techniques for Fast Elliptic Curve Cryptography on x86-64 based Processors" (Longa Gebotys 2010)](http://eprint.iacr.org/2010/335)
Ted1271 curve used as a basis for this one

##### [4] ["Elliptic and Hyperelliptic Curves: a Practical Security Analysis" (Bos Costello Miele 2013)](http://eprint.iacr.org/2013/644.pdf)
Discusses effect of morphisms on practical attacks in genus 1 and 2

##### [5] ["Twisted Edwards Curves Revisited" (Hisil Wong Carter Dawson 2008)](http://www.iacr.org/archive/asiacrypt2008/53500329/53500329.pdf)
Introduces Extended Twisted Edwards group laws

##### [6] [MAGMA Online Calculator](http://magma.maths.usyd.edu.au/calc/)
Calculating constants required, verification, etc

##### [7] ["Fault Attack on Elliptic Curve with Montgomery Ladder Implementation" (Fouque Lercier Real Valette 2008)](http://www.di.ens.fr/~fouque/pub/fdtc08.pdf)
Discussion on twist security

##### [8] ["Curve25519: new Diffie-Hellman speed records" (Bernstein 2006)](http://cr.yp.to/ecdh/curve25519-20060209.pdf)
Example of a conservative elliptic curve cryptosystem

##### [9] ["Safe Curves" (DJB 2013)](http://safecurves.cr.yp.to/)
Covers criterion for secure elliptic curves

##### [10] ["Families of fast elliptic curves from Q-curves" (Smith 2013)](http://eprint.iacr.org/2013/312.pdf)
New endomorphisms

##### [11] ["Signed Binary Representations Revisited" (Okeya et al 2004)](http://eprint.iacr.org/2004/195.pdf)
Introduces wMOF and explains wNAF - Neither produce regular patterns of zeroes so are not applicable

##### [12] ["Endomorphisms for Faster Elliptic Curve Cryptography on a Large Class of Curves" (Galbraith Lin Scott 2008)](http://eprint.iacr.org/2008/194)
Introduces GLS method

##### [13] ["Endomorphisms for Faster Elliptic Curve Cryptography on a Large Class of Curves" (Galbraith Lin Scott 2009)](http://www.iacr.org/archive/eurocrypt2009/54790519/54790519.pdf)
More information on GLS method
 
##### [14] ["EFD: Genus-1 curves over large-characteristic fields" (Lange et al)](http://www.hyperelliptic.org/EFD/g1p/index.html)
Explicit Forms Database with comparisons between cost of basic curve operations

##### [15] ["Easy scalar decompositions for efficient scalar multiplication on elliptic curves and genus 2 Jacobians" (Smith 2013)](http://hal.inria.fr/docs/00/87/49/25/PDF/easy.pdf)
Provides short bases for scalar decomposition for curves with efficient endomorphisms

##### [16] ["Faster Compact Diffie-Hellman: Endomorphisms on the x-line" (Costello Hisil Smith 2013)](http://eprint.iacr.org/2013/692.pdf)
Demonstrates several new techniques including Q-curve endomorphisms

##### [17] ["Fast and compact elliptic-curve cryptography" (Hamburg 2012)](http://eprint.iacr.org/2012/309.pdf)
Explores SAB-set comb multiplication and special prime moduli


## Credits

I appreciate suggestions and criticism from Mike Hamburg on an early version of the software.

