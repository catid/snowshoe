# Snowshoe project layout

Snowshoe is essentially one large C file, but it has been broken out into pieces so that each one can be separately unit-tested.

~~~
.
├── snowshoe.cpp
├── snowshoe.hpp
├── ecmul.inc
├── misc.inc
├── ecpt.inc
├── endo.inc
├── fe.inc
└── fp.inc
~~~

Each file includes the previous one:

+ `fp.inc` : Fp finite field arithmetic
+ `fe.inc` : Fp^2 optimal extension field, includes `fp.inc`
+ `endo.inc` : Endomorphism implementation, includes `fe.inc`
+ `ecpt.inc` : Elliptic curve point operations, includes `endo.inc`
+ `ecmul.inc` : Elliptic curve scalar multiplication, includes `ecpt.inc` and `misc.inc`
+ `snowshoe.cpp` : Defines library interface
+ `snowshoe.h` : Declares library interface

This way the unit testers can include e.g. `fp.inc` and use a minimal subset of the code to test those routines.

