# Change your compiler settings here

# Clang seems to produce faster code
#CCPP = g++
#CC = gcc
#OPTFLAGS = -O3 -fomit-frame-pointer -funroll-loops
CCPP = clang++
CC = clang
OPTFLAGS = -O4
DBGFLAGS = -g -O0 -DDEBUG
CFLAGS = -Wall -fstrict-aliasing -I./libcat -I./include
LIBNAME = libsnowshoe.a
LIBS =


# Object files

shared_test_o = Clock.o EndianNeutral.o

snowshoe_o = snowshoe.o

fp_test_o = fp_test.o $(shared_test_o)
fe_test_o = fe_test.o $(shared_test_o)
endo_test_o = endo_test.o $(shared_test_o)
ecpt_test_o = ecpt_test.o $(shared_test_o)
ecmul_test_o = ecmul_test.o $(shared_test_o)
snowshoe_test_o = snowshoe_test.o $(shared_test_o)


# Release target (default)

release : CFLAGS += $(OPTFLAGS)
release : library


# Debug target

debug : CFLAGS += $(DBGFLAGS)
debug : LIBNAME = libsnowshoe_debug.a
debug : library


# Library.ARM target

library.arm : CCPP = /Volumes/casedisk/prefix/bin/arm-unknown-eabi-g++
library.arm : CPLUS_INCLUDE_PATH = /Volumes/casedisk/prefix/arm-unknown-eabi/include
library.arm : CC = /Volumes/casedisk/prefix/bin/arm-unknown-eabi-gcc
library.arm : C_INCLUDE_PATH = /Volumes/casedisk/prefix/arm-unknown-eabi/include
library.arm : library


# Library target

library : CFLAGS += $(OPTFLAGS)
library : $(snowshoe_o)
	ar rcs $(LIBNAME) $(snowshoe_o)


# tester executables

fptest : CFLAGS += -DUNIT_TEST $(OPTFLAGS)
fptest : $(fp_test_o)
	$(CCPP) $(LIBS) -o fptest $(fp_test_o)

fetest : CFLAGS += -DUNIT_TEST $(OPTFLAGS)
fetest : $(fe_test_o)
	$(CCPP) $(LIBS) -o fetest $(fe_test_o)

endotest : CFLAGS += -DUNIT_TEST $(OPTFLAGS)
endotest : $(endo_test_o)
	$(CCPP) $(LIBS) -o endotest $(endo_test_o)

ecpttest : CFLAGS += -DUNIT_TEST $(OPTFLAGS)
ecpttest : $(ecpt_test_o)
	$(CCPP) $(LIBS) -o ecpttest $(ecpt_test_o)

ecmultest : CFLAGS += -DUNIT_TEST $(OPTFLAGS)
ecmultest : clean $(ecmul_test_o)
	$(CCPP) $(LIBS) -o ecmultest $(ecmul_test_o)
	./ecmultest

snowshoetest : CFLAGS += -DUNIT_TEST $(OPTFLAGS)
snowshoetest : clean $(snowshoe_test_o) library
	$(CCPP) $(LIBS) -L. -lsnowshoe -o snowshoetest $(snowshoe_test_o)
	./snowshoetest


# Shared objects

Clock.o : libcat/Clock.cpp
	$(CCPP) $(CFLAGS) -c libcat/Clock.cpp

EndianNeutral.o : libcat/EndianNeutral.cpp
	$(CCPP) $(CFLAGS) -c libcat/EndianNeutral.cpp


# Library objects

snowshoe.o : src/snowshoe.cpp
	$(CCPP) $(CFLAGS) -c src/snowshoe.cpp


# Executable objects

fp_test.o : tests/fp_test.cpp
	$(CCPP) $(CFLAGS) -c tests/fp_test.cpp

fe_test.o : tests/fe_test.cpp
	$(CCPP) $(CFLAGS) -c tests/fe_test.cpp

endo_test.o : tests/endo_test.cpp
	$(CCPP) $(CFLAGS) -c tests/endo_test.cpp

ecpt_test.o : tests/ecpt_test.cpp
	$(CCPP) $(CFLAGS) -c tests/ecpt_test.cpp

ecmul_test.o : tests/ecmul_test.cpp
	$(CCPP) $(CFLAGS) -c tests/ecmul_test.cpp

snowshoe_test.o : tests/snowshoe_test.cpp
	$(CCPP) $(CFLAGS) -c tests/snowshoe_test.cpp


# Cleanup

.PHONY : clean

clean :
	-rm fptest fetest endotest ecpttest ecmultest snowshoetest libsnowshoe.a $(shared_test_o) $(fp_test_o) $(fe_test_o) $(endo_test_o) $(ecpt_test_o) $(ecmul_test_o) $(snowshoe_test_o) $(snowshoe_o)

