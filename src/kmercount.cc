/* Works with kmers of length up to k=32 */
/* Uses a 64 bit hash */

#include "main.h"

#define FIXEDK31
const unsigned int k = 31;

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

uint64_t hash_full(unsigned int k, uint64_t kmer)
{
  /* compute 64 bit rolling hash of given k-mer from scratch */

  uint64_t hash = 0;

  for (unsigned int i = 0; i < k; i++)
    {
      /* rotate hash */
      hash = rotate_left_64(hash, shift_factor);
      /* insert value coming in */
      uint64_t in = (kmer >> ((k - 1 - i) << 1)) & 3;
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

static uint64_t duplicates = 0;
static uint64_t unique = 0;
FILE * fp_index = nullptr;

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
	  /* free slot, not seen before, insert new */
	  seqhashtable[seqhashindex].kmer = kmer;
	  seqhashtable[seqhashindex].count = 0;
	  unique++;
	  return;
	}

      if (kmerfound == kmer)
	{
	  /* slot in use, with match */
	  duplicates++;
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
  /* zero two upper bits */
  uint64_t * p = (uint64_t *) seq;
  uint64_t mem = *p++;
  kmer = mem;
  kmer &= 0x3fffffffffffffff;
  mem >>= 62;

  h = hash_full(k, kmer);
  if (bloomflex_get(bloom, h))
    hash_count(h, kmer, seqhashtable, seqhashsize);

  for(unsigned int i = k; i < seqlen; i++)
    {
      if ((i & 31) == 0)
	mem = *p++;

      uint64_t out = kmer & 3;
      uint64_t in = mem & 3;
      kmer >>= 2;
      kmer |= in << (2*(k-1));
      mem >>= 2;

#ifdef FIXEDK31
      h = hash_update_31(h, out, in);
#else
      h = hash_update(k, h, out, in);
#endif

      if (bloomflex_get(bloom, h))
	hash_count(h, kmer, seqhashtable, seqhashsize);
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
      for (unsigned int i = 0; i < seqlen; i++)
	{
	  uint64_t in = nt_extract(seq, i);
	  uint64_t out = kmer & 3;
	  kmer >>= 2ULL;
	  kmer |= (in << (2ULL * (k - 1ULL)));

	  if (i >= k - 1)
	    {
	      if (i == k - 1)
		h = hash_full(k, kmer);
	      else
#ifdef FIXEDK31
		h = hash_update_31(h, out, in);
#else
		h = hash_update(k, h, out, in);
#endif

	      bloomflex_set(bloom, h);
	      hash_insert(h, kmer, seqhashtable, seqhashsize);
	    }
	}
    }
}

void kmercount(const char * kmer_filename, const char * seq_filename)
{
  (void) kmer_filename;

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
      seqhashtable[j].kmer = -1;
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

  fprintf(logfile,   "Unique kmers:      %llu\n", unique);
  if (duplicates > 0)
    fprintf(logfile, "Duplicates:        %llu\n", duplicates);

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

  fprintf(logfile, "\n");


  /* Print kmers and counts to output file */
  unsigned int x = 0;
  uint64_t y = 0;
  progress_init("Writing results", seqhashsize);
  for (uint64_t i = 0; i < seqhashsize; i++)
    {
      struct hashentry * e = seqhashtable + i;
      if (e->kmer != (uint64_t)-1)
	if (e->count > 0)
	  {
#if 0
	    fprintf(outfile, ">kmer%u\n", x+1);
	    fprintseq(outfile, e->kmer);
	    fprintf(outfile, "\n");
#else
	    fprintseq(outfile, e->kmer);
	    fprintf(outfile,"\t%llu\n", e->count);
#endif
	    x++;
	    y += e->count;
	  }
      progress_update(i);
    }
  progress_done();

  fprintf(logfile, "Matching kmers:    %u\n", x);
  fprintf(logfile, "Total matches:     %llu\n", y);

  delete [] seqhashtable;
  seqhashtable = nullptr;
  bloomflex_exit(bloom);
  db_free(seq_db);
}
