# Snowshoe project layout

Snowshoe is essentially one large C file, but it has been broken out into pieces so that each one can be separately unit-tested.

~~~
.
├── Snowshoe.cpp
├── Snowshoe.hpp
├── ecmul.cpp
├── ecpt.cpp
├── endo.cpp
├── fe.cpp
└── fp.cpp
~~~

Each file includes the previous one:

+ `fp.cpp` : Fp finite field arithmetic
+ `fe.cpp` : Fp^2 optimal extension field, includes `fp.cpp`
+ `endo.cpp` : Endomorphism implementation, includes `fe.cpp`
+ `ecpt.cpp` : Elliptic curve point operations, includes `endo.cpp`
+ `ecmul.cpp` : Elliptic curve scalar multiplication, includes `ecpt.cpp`
+ `Snowshoe.cpp` : Defines library interface
+ `Snowshoe.hpp` : Declares library interface

This way the unit testers can include e.g. `fp.cpp` and use a minimal subset of the code to test those routines.

