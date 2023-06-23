/*
    Copyright (C) 2023 Torbjorn Rognes

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

/* Kmercounter using Bloom filter and rapid hash function */

/* Works with kmers of length up to k=32 */
/* Uses a 64 bit hash */

/*
  k=31
  using bits 0-61, 2 bits per nucleotide
  bit 0-1 encodes the first (1) nucleotide in the sequence
  bit 60-61 encodes the last (31) nucleotide
  A=00, C=01 G=10 T=11

  Example sequence:
  AAGAAATGAGAAGTAATCAGAAAACCACTTAAGG

  Kmer 1:                         Encoding:
  AAGAAATGAGAAGTAATCAGAAAACCACTTA 0f4500870e08b020

  Kmer 2:
  AGAAATGAGAAGTAATCAGAAAACCACTTAA 03d14021c3822c08

  Kmer 3:
  GAAATGAGAAGTAATCAGAAAACCACTTAAG 20f4500870e08b02

  Kmer 4:
  AAATGAGAAGTAATCAGAAAACCACTTAAGG 283d14021c3822c0

*/

#include "main.h"

static unsigned int k = 31; // Default 31

static uint64_t unique = 0;

struct hashentry
{
  uint64_t kmer;
  uint64_t count;
};

static const unsigned int shift_factor = 2;

static const uint64_t hashvalues[4] =
  {
    /* These pseudo-random constants should perhaps be choosen wisely? */
    0xba64e57c490e2ef4,
    0x4938a808abe1edcf,
    0x715849e4da68576a,
    0x02db58f212586265
  };

static const uint64_t hashvalues_rot60[4] =
  {
    0x4ba64e57c490e2ef,
    0xf4938a808abe1edc,
    0xa715849e4da68576,
    0x502db58f21258626
  };

inline uint64_t rotate_left_64(uint64_t x, unsigned int r)
{
  /* rotate value in x by r places to the left */
  /* 0 <= r <= 64 */
  /* operates on 64 bit values */

  return (x << r) | (x >> (64 - r));
}

inline uint64_t reverse_nucleotides(unsigned int k, uint64_t kmer)
{
  uint64_t res = 0;
  uint64_t x = kmer;
  for (unsigned int i = 0; i < k; i++)
    {
      uint64_t t = x & 3;
      x >>= 2;
      res <<= 2;
      res |= t;
    }
  return res;
}

uint64_t hash_full(unsigned int k, uint64_t kmer)
{
  /* compute 64 bit rolling hash of given k-mer from scratch */

  uint64_t hash = 0;

  for (unsigned int i = 0; i < k; i++)
    {
      /* rotate hash */
      hash = rotate_left_64(hash, shift_factor);

      /* insert value coming in */
      uint64_t in = (kmer >> (i << 1)) & 3;
      hash ^= hashvalues[in];
    }

  return hash;
}

inline uint64_t hash_update(uint64_t k, uint64_t h, uint64_t out, uint64_t in)
{
  /* update 64 bit rolling hash with a new nucleotide */
  uint64_t hash = h;

  /* remove value going out (rotation may be precomputed for fixed k)*/
  hash ^= rotate_left_64(hashvalues[out], shift_factor * (k - 1ULL));

  /* rotate hash */
  hash = rotate_left_64(hash, shift_factor);

  /* insert value coming in */
  hash ^= hashvalues[in];

  return hash;
}

inline uint64_t hash_update_31(uint64_t h, uint64_t out, uint64_t in)
{
  /* update 64 bit rolling hash with a new nucleotide */
  uint64_t hash = h;

  /* remove value going out (rotation may be precomputed for fixed k)*/
  hash ^= hashvalues_rot60[out];

  /* rotate hash */
  hash = rotate_left_64(hash, shift_factor);

  /* insert value coming in */
  hash ^= hashvalues[in];

  return hash;
}

void hash_insert(uint64_t hash,
		 uint64_t kmer,
		 hashentry *  seqhashtable,
		 uint64_t seqhashsize)
{
  uint64_t seqhashindex = hash % seqhashsize;

  while (1)
    {
      uint64_t kmerfound = seqhashtable[seqhashindex].kmer;

      if (kmerfound == (uint64_t) - 1)
	{
	  /* free slot, not seen before, insert new, zero count */
	  seqhashtable[seqhashindex].kmer = kmer;
	  seqhashtable[seqhashindex].count = 0;
	  unique++;
	  return;
	}

      if (kmerfound == kmer)
	{
	  /* slot in use, with match */
	  return;
	}

      /* in use, but no match, try next */
      seqhashindex = (seqhashindex + 1) % seqhashsize;
    }
}

inline void hash_count(uint64_t hash,
		       uint64_t kmer,
		       hashentry *  seqhashtable,
		       uint64_t seqhashsize)
{
  uint64_t seqhashindex = hash % seqhashsize;

  while (1)
    {
      uint64_t kmerfound = seqhashtable[seqhashindex].kmer;

      if (kmerfound == (uint64_t) - 1)
	{
	  /* no match, ignore this kmer */
	  return;
	}
      else if (kmerfound == kmer)
	{
	  /* match, count it */
	  seqhashtable[seqhashindex].count++;
	  return;
	}

      /* in use, not matching, try next bucket */
      seqhashindex = (seqhashindex + 1) % seqhashsize;
    }
}

void kmer_check(unsigned int seqlen, char * seq, bloomflex_s * bloom, hashentry * seqhashtable, uint64_t seqhashsize)
{
  uint64_t kmer = 0;
  uint64_t h = 0;

  if (seqlen < k)
    return;

  /* first kmer */
  uint64_t * p = (uint64_t *) seq;
  uint64_t mem = *p++;
  kmer = mem;
  kmer &= 0x3fffffffffffffff;
  mem >>= 62;

  h = hash_full(k, kmer);
  if (bloomflex_get(bloom, h))
    hash_count(h, kmer, seqhashtable, seqhashsize);

  if (k == 31)
    {
      // tailored for k=31
      for(unsigned int i = 31; i < seqlen; i++)
	{
	  if ((i & 31) == 0)
	    mem = *p++;

	  uint64_t out = kmer & 3;
	  uint64_t in = mem & 3;
	  kmer >>= 2;
	  kmer |= in << 60;
	  mem >>= 2;

	  h = hash_update_31(h, out, in);
	  if (bloomflex_get(bloom, h))
	    hash_count(h, kmer, seqhashtable, seqhashsize);
	}
    }
  else
    {
      for(unsigned int i = k; i < seqlen; i++)
	{
	  if ((i & 31) == 0)
	    mem = *p++;

	  uint64_t out = kmer & 3;
	  uint64_t in = mem & 3;
	  kmer >>= 2;
	  kmer |= in << (2*(k-1));
	  mem >>= 2;

	  h = hash_update(k, h, out, in);
	  if (bloomflex_get(bloom, h))
	    hash_count(h, kmer, seqhashtable, seqhashsize);
	}
    }
}

void fprintseq(FILE * fp, uint64_t kmer)
{
  char sym_nt[5] = "ACGT";
  char buffer[33];
  for (unsigned int i = 0; i < k; i++)
    buffer[i] = sym_nt[(kmer >> 2*i) & 3];
  buffer[k] = 0;
  fprintf(fp, "%s", buffer);
}


void kmer_insert(unsigned int seqlen, char * seq, bloomflex_s * bloom, hashentry * seqhashtable, uint64_t seqhashsize)
{
  /* find and hash all kmers in a sequence */
  /* 1 <= k <= 32 */

  uint64_t kmer = 0;
  uint64_t h = 0;

  if (seqlen == k)
    {
      kmer = *((uint64_t*) seq);
      h = hash_full(k, kmer);
      bloomflex_set(bloom, h);
      hash_insert(h, kmer, seqhashtable, seqhashsize);
    }
  else
    {
      fprintf(logfile, "\nFatal error: Sequence length (%u) is different from given k (%u).\n", seqlen, k);
      exit(1);
    }
}

int compare_kmers(const void * a, const void * b)
{
  const struct hashentry * x = (struct hashentry *)(a);
  const struct hashentry * y = (struct hashentry *)(b);

  // Compare kmer counts, sort by decending order
  if (x->count > y->count)
    return -1;
  else if (x->count < y->count)
    return +1;
  else
    {
#if 1
      if (x->kmer < y->kmer)
	return -1;
      else if (x->kmer > y->kmer)
	return +1;
      else
	return 0;
#else
      /* Too slow */
      /* Equal count, sort by sequence */
      uint64_t e = reverse_nucleotides(k, x->kmer);
      uint64_t f = reverse_nucleotides(k, y->kmer);
      if (e < f)
	return -1;
      else if (e > f)
	return +1;
      else
	return 0;
#endif
    }
}

void print_results(hashentry * seqhashtable, uint64_t seqhashsize)
{
  fprintf(logfile, "\n");

  progress_init("Sorting results:  ", 1);
  qsort(seqhashtable,
	seqhashsize,
	sizeof(struct hashentry),
	compare_kmers);
  progress_done();

  /* Print kmers and counts to output file */
  unsigned int x = 0;
  uint64_t y = 0;
  progress_init("Writing results:  ", seqhashsize);
  for (uint64_t i = 0; i < seqhashsize; i++)
    {
      struct hashentry * e = seqhashtable + i;
      if (e->kmer != (uint64_t)-1)
	{
	  if (e->count > 0)
	  {
	    fprintseq(outfile, e->kmer);
	    fprintf(outfile,"\t%" PRIu64 "\n", e->count);
	    x++;
	    y += e->count;
	  }
	  else
	    break;
	}
      else
	break;
      progress_update(i);
    }
  progress_done();

  fprintf(logfile, "Matching kmers:    %u\n", x);
  fprintf(logfile, "Total matches:     %" PRIu64 "\n", y);
}


void kmercount(const char * kmer_filename,
	       const char * seq_filename,
	       int opt_k)
{
  k = opt_k;

  /* Read FASTA with kmers */
  fprintf(logfile, "Reading kmer file\n");
  struct db_s * kmer_db = db_read(kmer_filename);
  unsigned int kmer_count = db_getsequencecount(kmer_db);

  /* set up Bloom filter, 1 byte per kmer, 4 of 8 bits set */
  bloomflex_s * bloom = bloomflex_init(kmer_count, 4);

  /* set up hashtable */
  const uint64_t seqhashsize = 2 * kmer_count;
  struct hashentry * seqhashtable = nullptr;
  seqhashtable = new hashentry [seqhashsize] { };
  for (uint64_t j = 0; j < seqhashsize; j++)
    {
      seqhashtable[j].count = 0;
      seqhashtable[j].kmer = (uint64_t) -1;
    }

  /* compute hash for all kmers and store them in bloom & hash table */
  progress_init("Indexing kmers:   ", kmer_count);
  for(unsigned int i = 0; i < kmer_count; i++)
    {
      char * seq;
      unsigned int seqlen;
      db_getsequenceandlength(kmer_db, i, & seq, & seqlen);
      kmer_insert(seqlen, seq, bloom, seqhashtable, seqhashsize);
      progress_update(i);
    }
  progress_done();

  fprintf(logfile,   "Unique kmers:      %" PRIu64 "\n", unique);

  db_free(kmer_db);

  fprintf(logfile, "\n");

  /* Read FASTA sequence file */
  fprintf(logfile, "Reading sequence file\n");
  struct db_s * seq_db = db_read(seq_filename);
  int seq_count = db_getsequencecount(seq_db);
  uint64_t seq_nucleotides = db_getnucleotides(seq_db);

  /* Compute hash for all kmers in db and count */
  progress_init("Counting matches: ", seq_nucleotides);
  uint64_t nt_processed = 0;
  for(int i = 0; i < seq_count; i++)
    {
      char * seq;
      unsigned int seqlen;
      db_getsequenceandlength(seq_db, i, & seq, & seqlen);
      kmer_check(seqlen, seq, bloom, seqhashtable, seqhashsize);
      nt_processed += seqlen;
      progress_update(nt_processed);
    }
  progress_done();

  print_results(seqhashtable, seqhashsize);

  delete [] seqhashtable;
  seqhashtable = nullptr;
  bloomflex_exit(bloom);
  db_free(seq_db);
}
