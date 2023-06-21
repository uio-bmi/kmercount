/* Works with kmers of length up to k=32 */
/* Uses a 64 bit hash */

#include "main.h"


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

bool is_indexing = false;
static uint64_t duplicates = 0;
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

      if (kmerfound == kmer)
	{
	  /* slot in use, with match */
	  duplicates++;
	  seqhashtable[seqhashindex].count++;
	  return;
	}

      if (kmerfound == (uint64_t) - 1)
	{
	  /* free slot, not seen before, insert new */
	  seqhashtable[seqhashindex].kmer = kmer;

	  if (is_indexing)
	    {
	      fwrite(& kmer, 8, 1, fp_index);
	      seqhashtable[seqhashindex].count = 1;
	    }
	  else
	    {
	      seqhashtable[seqhashindex].count = 0;
	    }
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
  /* find and hash all kmers in a sequence */
  /* 1 <= k <= 32 */

  const unsigned int k = 31;
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

      //h = hash_update(k, h, out, in);
      h = hash_update_31(h, out, in);
      if (bloomflex_get(bloom, h))
	hash_count(h, kmer, seqhashtable, seqhashsize);
    }
}

void kmer_insert(unsigned int seqlen, char * seq, bloomflex_s * bloom, hashentry * seqhashtable, uint64_t seqhashsize)
{
  (void) bloom;

  /* find and hash all kmers (minimizers) in a sequence */
  /* 1 <= k <= 32 */

  const unsigned int k = 31;
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
		//h = hash_update(k, h, out, in);
		h = hash_update_31(h, out, in);

	      bloomflex_set(bloom, h);
	      hash_insert(h, kmer, seqhashtable, seqhashsize);
	    }
	}
    }
}

void prepare_index(struct Parameters const & p)
{
  /* Read FASTA with kmers, produce index file */
  db_read(p.input_filename.c_str(), p);

  fp_index = fopen("index.bin", "w");
  if (!fp_index)
    fatal("Cant open index.bin");

  unsigned int seqcount = db_getsequencecount();

  bloomflex_s * bloom = bloomflex_init(seqcount, 8);

  /* set up hashtable */

  const uint64_t seqhashsize = 2 * seqcount;
  struct hashentry * seqhashtable = nullptr;
  seqhashtable = new hashentry [seqhashsize] { };
  for (uint64_t j = 0; j < seqhashsize; j++)
    {
      seqhashtable[j].count = 0;
      seqhashtable[j].kmer = -1;
    }


  /* compute hash for all amplicons and store them in a hash table */

  progress_init("Analyse sequences:", seqcount);

  for(unsigned int i = 0; i < seqcount; i++)
    {
      char * seq;
      unsigned int seqlen;
      db_getsequenceandlength(i, & seq, & seqlen);
      kmer_insert(seqlen, seq, bloom, seqhashtable, seqhashsize);
      progress_update(i);
    }

  progress_done();

  printf("Duplicates: %llu\n", duplicates);

  delete [] seqhashtable;
  seqhashtable = nullptr;

  bloomflex_exit(bloom);

  fclose(fp_index);
}

void perform_kmer_count(struct Parameters const & p)
{
  (void)p;

  /* Read index file */

  struct stat buf;
  int ret = stat("index.bin", &buf);
  if (ret != 0)
    fatal("Can not stat index.bin");

  fp_index = fopen("index.bin", "r");
  if (!fp_index)
    fatal("Can not open index.bin");

  int kmercount = buf.st_size / 8;

  printf("Reading %d kmers...", kmercount);
  uint64_t * kmers = (uint64_t *) xmalloc(kmercount * 8);
  fread(kmers, 8, kmercount, fp_index);
  printf(" Done\n");

  /* Bloom filter, 1 byte per kmer */
  bloomflex_s * bloom = bloomflex_init(kmercount, 4);

  /* set up hashtable, 2 slots per kmer */

  const uint64_t seqhashsize = 2 * kmercount;
  struct hashentry * seqhashtable = nullptr;
  seqhashtable = new hashentry [seqhashsize] { };
  for (uint64_t j = 0; j < seqhashsize; j++)
    {
      seqhashtable[j].count = 0;
      seqhashtable[j].kmer = -1;
    }


  /* compute hash for all amplicons and store them in a hash table */

  progress_init("Inserting kmers/hash:", kmercount);
  for(int i = 0; i < kmercount; i++)
    {
      kmer_insert(31, (char*) (kmers + i), bloom, seqhashtable, seqhashsize);
      progress_update(i);
    }
  progress_done();


  /* Read FASTA */
  printf("Reading FASTA file...");
  db_read(p.input_filename.c_str(), p);
  int db_seqcount = db_getsequencecount();
  printf("Done\n");


  /* Compute hash for all kmers in db and count */
  progress_init("Analyse sequences:", db_seqcount);
  for(int i = 0; i < db_seqcount; i++)
    {
      char * seq;
      unsigned int seqlen;
      db_getsequenceandlength(i, & seq, & seqlen);
      kmer_check(seqlen, seq, bloom, seqhashtable, seqhashsize);
      progress_update(i);
    }
  progress_done();


  /* Print out kmer, hash, and count, only 50 first */
  int x = 0;
  for (uint64_t i = 0; i < seqhashsize; i++)
    {
      struct hashentry * e = seqhashtable + i;
      if ((e->kmer != (uint64_t)-1) && (e->count > 0))
	{
	  printf("kmer %016llx hash %016llx count %llu\n",
		 e->kmer, hash_full(31, e->kmer), e->count);
	  x++;
	}
      if (x >= 1)
	break;
    }

  delete [] seqhashtable;
  seqhashtable = nullptr;
  bloomflex_exit(bloom);
  fclose(fp_index);
}

void kmercount(struct Parameters const & p)
{
  if (0)
    {
      is_indexing = true;
      prepare_index(p);
    }
  else
    {
      is_indexing = false;
      perform_kmer_count(p);
    }
}
