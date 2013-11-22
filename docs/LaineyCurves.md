# Lainey Curves

Introducing a sub-class of Q-curves that offer efficient arithmetic and also
efficient 2-dimensional endomorphisms that do not accelerate the Pollard rho
attack (according to [4]) but do offer fast scalar multiplication.

Definition:

~~~
Choose a base field Fp, ie. p = 2^127-1.

Let D be a non-square in Fp, so that i := sqrt(D), ie. D = -1.

By [10, sect.5] there is a class of curves of the Weierstrass form:

	y^2 = x^3 + Wa * x + Wb  over Fp^2

Wa = -6 * (5 - 3 * s * i)
Wb = 8 * (7 - 9 * s * i)
~~~

Note that `s` is a free parameter and is varied while searching for curves of
cryptographically interesting group order.  Also note that `s` is a complex
number.

This curve has an efficient 2-dimensional endomorphism that can be even more
efficiently computed in twisted Edwards coordinates:

~~~
By [10, sect.8], this curve can be birationally equivalent to a
twisted Edwards curve of the form:

	E : a * x^2 + y^2 = 1 + d * x^2 * y^2

C = 9 * (1 + s * i)

B = SquareRoot(2 * C)
A = 12 / B

a = (A + 2) / B
d = (A - 2) / B
~~~

Based on the efficient group laws from [5] we would like `a = -1`.

Since the twisted Edwards form is much nicer for performing point operations,
it is interesting to derive the value of `s` from a desired value for `a`:

~~~
Given a:

B = (2 + SquareRoot(4 + 48 * a)) / (2 * a)

C = B * B / 2

s = (C / 9 - 1) / i
~~~

Since `a` and `d` are linked to the choice of s, it is extremely unlikely that
`a = -1` will be useable.  And in practice it is not.

Further, using formulas other than `a = -1` introduces an unwanted extra Fp^2
multiplication into the Twisted Edwards group laws in addition to a
multiplication by `a`.

However, for special values of `a` there is a helpful change of variables:

Dedicated Doubling group law from [5]:

~~~
X3 = 2 * X1 * Y1 * (2 * Z1^2 - Y1^2 - a * X1^2)
Y3 = (Y1^2 + a * X1^2) * (Y1^2 - a * X1^2)
T3 = 2 * X1 * Y1 * (Y1^2 - a * X1^2)
Z3 = (Y1^2 + a * X1^2) * (2 * Z1^2 - Y1^2 - a * X1^2)
~~~

With the change of variables X1' = X1 * sqrt(u), these formulae compute instead:

~~~
X3' = SquareRoot(u) * 2 * X1 * Y1 * (2 * Z1^2 - Y1^2 - a * u * X1^2)
Y3' = (Y1^2 + a * u * X1^2) * (Y1^2 - a * u * X1^2)
T3' = SquareRoot(u) * 2 * X1 * Y1 * (Y1^2 - a * u * X1^2)
Z3' = (Y1^2 + a * u * X1^2) * (2 * Z1^2 - Y1^2 - a * u * X1^2)
~~~

Note we can select a = -1 for more efficient evaluation.

Compared to the original coordinates:

~~~
X3' = SquareRoot(u) * X3  (so the change of variables is maintained consistently)
Y3' = Y3
T3' = SquareRoot(u) * T3
Z3' = Z3
~~~

And we have almost recovered the original result except that `T3` is multiplied
by `SquareRoot(u)`.

If we choose `u = (v^-1)^2`, for some small v, then we can fix `T3'`
after our favorite implementation of the group laws finishes by multiplying `T3`
by `v`, which effectively removes the unwanted `SquareRoot(u)` term.

Interestingly if `T3` is not needed, we can skip the correction step entirely.

Similarly for the Dedicated Addition group law from [5]:

~~~
X3 = (X1 * Y2 - Y1 * X2) * (T1 * Z2 + Z1 * T2)
Y3 = (Y1 * Y2 + a * X1 * X2) * (T1 * Z2 - Z1 * T2)
T3 = (T1 * Z2 + Z1 * T2) * (T1 * Z2 - Z1 * T2)
Z3 = (Y1 * Y2 + a * X1 * X2) * (X1 * Y2 - Y1 * X2)
~~~

With the change of variables `X1' = X1 * SquareRoot(u)`, these formulae compute:

~~~
X3 = SquareRoot(u) * (X1 * Y2 - Y1 * X2) * (T1 * Z2 + Z1 * T2)
Y3 = (Y1 * Y2 + a * u * X1 * X2) * (T1 * Z2 - Z1 * T2)
T3 = (T1 * Z2 + Z1 * T2) * (T1 * Z2 - Z1 * T2)
Z3 = SquareRoot(u) * (Y1 * Y2 + a * u * X1 * X2) * (X1 * Y2 - Y1 * X2)
~~~

Compared to the original coordinates:

~~~
X3' = SquareRoot(u) * X3  (so the change of variables is maintained consistently)
Y3' = Y3
T3' = T3
Z3' = SquareRoot(u) * Z3
~~~

This time `Z3` will require correction by multiplying by `v` as before.

The Unified Addition law from [5] can be adapted similarly.

Lainey curves may be found with a MAGMA script such as magma_laineygen.txt:

~~~
p := 2^127 - 1;
K<i> := GF(p^2);

ds := 12;

print ds;

pwr := 2 ^ ds;

h := pwr*i;

hi := (h^-1);
u := hi^2;
a := -u;
B := (2 + SquareRoot(4 + 48 * a)) / (2 * a);
C := B * B / 2;
s := (C / 9 - 1) / i;

// Validate:
//C1 := 9 * (1 + s * i);
//B1 := SquareRoot(2 * C);
//a1 := (12 / B1 + 2) / B1;
// Note a1 != a if guess +/- in B step wrong.  This is fine because we end up
// on the same curve either way.

wa := -6 * (5 - 3 * s * i);
wb := 8 * (7 - 9 * s * i);

E := EllipticCurve([K | wa, wb]);

time EP := Order(E);

if EP mod 8 eq 0 then
   TEST := IsPrime(EP div 8);
else
   if EP mod 4 eq 0 then
      TEST := IsPrime(EP div 4);
   else
      if EP mod 2 eq 0 then
           TEST := IsPrime(EP div 2);
      else
           TEST := IsPrime(EP);
      end if;
   end if;
end if;

print h;
if TEST eq true then
   print "FOUND A LAINEY CURVE!";
else
   print "bad curve";
end if;
~~~

As a side note, I initially attempted to work with GLS curves because they
seemed to offer nicer values for `a`.  The problem with GLS curves is that the
value for `a` must be a non-square in Fp, which means this trick cannot be
applied to GLS curves.  Once I realized this, I attempted to apply the trick
to Q-curves with endomorphisms and was delighted to find that it works, and
since my girlfriend has some nice curves, it seemed suitable to name this
new sub-class of Q-curves after her.

The first Lainey curve:

~~~
hi =
84904762116725386595141009873777590016*i
u =
43696183973333720444381044205952108540
a =
126444999487135511287306259509931997187
B =
113969212755026911009472686468072940927*i + 2105352
C =
133072242013244136953902226708658819287*i + 4432519676016
s =
170141183460469231731687303223381919504*i +
147117836248503195452857039191094173264

C1 =
133072242013244136953902226708658819287*i + 4432519676016

B1 =
113969212755026911009472686468072940927*i + 2105352

a1 =
126444999487135511287306259509931997187

d1 =
37604001849568726182388673756441218534*i +
69731271666979100710077158271303961945

wa =
96003300566019042176117149701433532847*i + 8865039351984
wb =
126410348117331526490593312341918185793*i +
170141183460469231731687268255726697727

EP =
28948022309329048855892746252171976962752096054440303231503648918509526300376

[ <2, 3>, <3618502788666131106986593281521497120344012006805037903937956114813690787547, 1>; ]

r = 8 * q

q = 0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFCAD47CC7B804E4DA6654BD32BAEB36DB
~~~

For these curves with CM discriminant = 2, the cofactor must be 8 instead of 4.  This is not a huge deal, but the bigger problem is that I have been unable to get the Q-curve endomorphism to work for values of s in Q.  Smith's paper indicates it *should* work, but I must be doing something wrong.  So I have had to abandon this project and stick with slower and less secure GLS curves.  At some point I will revisit this and maybe implement a 4-GLV method using the Lainey curve trick, which would be a new software speed record by a small margin AFAICS.

