# Makefile for kmercount

# Copyright (C) 2012-2023 Torbjorn Rognes and Frederic Mahe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
# Department of Informatics, University of Oslo,
# PO Box 1080 Blindern, NO-0316 Oslo, Norway

COMMON = -g
EXTRAOBJ =
LINKOPT =
LIBS = -lpthread
WARNINGS = -Wall -Wextra
CLANG = clang

# Run "make RELEASE=1" to compile for release
ifdef RELEASE
	COMMON += -DNDEBUG
endif

# Run "make PROFILE=1" to compile for profiling
ifdef PROFILE
	COMMON += -pg
endif

# Run "make COVERAGE=1" to compile for coverage tests
ifdef COVERAGE
	COMMON += -fprofile-arcs -ftest-coverage -O0
	LIBS += -lgcov
	LINKFLAGS += --coverage
else
	COMMON += -O3
endif

# Machine specific
ifeq ($(CXX), aarch64-linux-gnu-g++)
	MACHINE = aarch64
else ifeq ($(CXX), powerpc64le-linux-gnu-g++)
	MACHINE = ppc64le
else
	MACHINE = $(shell uname -m)
endif

ifeq ($(MACHINE), x86_64)
	COMMON += -march=x86-64 -mtune=generic -std=c++11
else ifeq ($(MACHINE), aarch64)
	COMMON += -march=armv8-a+simd -mtune=generic \
	          -flax-vector-conversions -std=c++11
else ifeq ($(MACHINE), arm64)
	COMMON += -std=c++11
else ifeq ($(MACHINE), ppc64le)
	COMMON += -mcpu=power8 -maltivec -std=gnu++11
endif

# OS specific
ifeq ($(CXX), x86_64-w64-mingw32-g++)
	LIBS += -lpsapi
	LINKOPT += -static
	BIN = kmercount.exe
else
	WARNINGS += -pedantic
	BIN = kmercount
endif

LINKFLAGS = $(COMMON) $(LINKOPT)

CXXFLAGS = $(COMMON) $(WARNINGS)

PROG = kmercount

OBJS = arch.o bloomflex.o db.o main.o util.o fatal.o kmercount.o

DEPS = Makefile \
	arch.h bloomflex.h db.h pseudo_rng.h main.h util.h fatal.h

all : $(PROG)

kmercount : $(OBJS) $(DEPS)
	$(CXX) $(LINKFLAGS) -o $@ $(OBJS) $(LIBS)

clean :
	rm -f kmercount *.o *~ gmon.out *.gcno *.gcda *.gcov

.o : .cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<
