/*
    kmercount

    Copyright (C) 2012-2023 Torbjorn Rognes and Frederic Mahe

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

#include <cinttypes>

#ifndef PRIu64
#ifdef _WIN32
#define PRIu64 "I64u"
#else
constexpr char PRIu64[] = "lu";
#endif
#endif

#ifndef PRId64
#ifdef _WIN32
#define PRId64 "I64d"
#else
constexpr char PRId64[] = "ld";
#endif
#endif

#include <algorithm>
#include <array>
#include <cassert>
#include <climits>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>  // std::exit
#include <cstring>
#include <fcntl.h>
#include <getopt.h>
#include <iostream>
#include <random>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>  // replace with <fstream> to improve portability
#include <vector>

#ifdef __APPLE__
#include <sys/resource.h>
#include <sys/sysctl.h>
#elif defined _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#include <sys/sysinfo.h>
#endif

#ifdef __aarch64__
#include <arm_neon.h>
#elif defined __x86_64__

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#ifdef __SSSE3__
#include <tmmintrin.h>  // could be removed?
#endif

#ifdef __SSE4_1__
#include <smmintrin.h>  // could be removed?
#endif

#ifdef __POPCNT__
#include <popcntintrin.h>
#endif

#elif defined __PPC__

#ifdef __LITTLE_ENDIAN__
#include <altivec.h>
#else
#error Big endian ppc64 CPUs not supported
#endif

#else

#error Unknown architecture
#endif

static_assert(INT_MAX > INT16_MAX, "Your compiler uses very short integers.");

#include "arch.h"
#include "bloomflex.h"
#include "db.h"
#include "fatal.h"
#include "pseudo_rng.h"
#include "util.h"


/* constants */

const std::string program_version {"0.0.1"};
constexpr char sepchar {' '};
constexpr char dash_filename {'-'};


/* structures and data types */

struct queryinfo
{
  uint64_t qno;
  int64_t len;
  char * seq;
};

using queryinfo_t = struct queryinfo;
extern queryinfo_t query;

/* common data */

struct Parameters {
  bool opt_help {false};
  bool opt_version {false};
  std::string input_filename {dash_filename};
  std::string opt_output_file {dash_filename};
};

extern std::string opt_log;  // used by multithreaded functions
extern int64_t opt_threads;

extern std::FILE * outfile;
extern std::FILE * logfile;


/* inline functions */

inline auto nt_extract(char * seq, uint64_t pos) -> unsigned char
{
  // Extract compressed nucleotide in sequence seq at position pos
  static constexpr unsigned int max_nt_per_uint64 {32};  // 32 nt fit in 64 bits
  static constexpr unsigned int drop_remainder {5};  // (len+31) % 32 (drop remainder)
  static constexpr unsigned int max_range {3};
  // outputs four possible values: 0, 1, 2 or 3
  return (((reinterpret_cast<uint64_t*>(seq))[pos >> drop_remainder]) >> \
          ((pos & (max_nt_per_uint64 - 1)) << 1)) & max_range;
}

inline auto nt_bytelength(unsigned int len) -> unsigned int
{
  // Compute number of bytes used for compressed sequence of length len
  // (minimum result is 8 bytes)
  static constexpr unsigned int max_nt_per_uint64 {32};  // 32 nt fit in 64 bits
  static constexpr unsigned int drop_remainder {5};  // (len + 31) % 32 (drop remainder)
  static constexpr unsigned int bytes_per_uint64 {8};  // times 8 to get the number of bytes
  return ((len + max_nt_per_uint64 - 1) >> drop_remainder) * bytes_per_uint64;
}

/* functions in kmercount.cc */

void kmercount(struct Parameters const & p);
