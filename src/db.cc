/*
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

#include "main.h"


constexpr unsigned int memchunk {1 << 20};  // 1 megabyte
constexpr unsigned int linealloc {2048};

static const signed char map_nt[256] =
  {
    // AaNn = 0, Cc = 1, Gg = 2, TtUu = 3, rest = -1
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1,  0, -1,
    -1, -1, -1, -1,  3,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1,  0, -1,
    -1, -1, -1, -1,  3,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
  };

struct seqinfo_s
{
  char * seq;
  unsigned int seqlen;
};

struct db_s
{
  unsigned int sequences;
  uint64_t nucleotides;
  unsigned int longest;
  char * datap;
  struct seqinfo_s * seqindex {nullptr};
};

unsigned int db_getsequencecount(struct db_s * d)
{
  return d->sequences;
}

struct db_s * db_read(const char * filename)
{
  struct db_s * d = (struct db_s *) xmalloc(sizeof(struct db_s));

  d->sequences = 0;
  d->nucleotides = 0;
  d->longest = 0;
  d->datap = nullptr;
  d->seqindex = nullptr;

  /* allocate space */

  uint64_t dataalloc {memchunk};
  d->datap = static_cast<char *>(xmalloc(dataalloc));
  uint64_t datalen {0};

  /* open input file or stream */

  assert(filename != nullptr);  // filename is set to '-' (stdin) by default

  std::FILE * input_fp { fopen_input(filename) };
  if (input_fp == nullptr)
    {
      fatal(error_prefix, "Unable to open input data file (", filename, ").\n");
    }

  /* get file size */

  struct stat fs;

  if (fstat(fileno(input_fp), & fs) != 0)
    {
      fatal(error_prefix, "Unable to fstat on input file (", filename, ").\n");
    }
  bool is_regular = S_ISREG(fs.st_mode);
  uint64_t filesize = is_regular ? fs.st_size : 0;
  uint64_t filepos = 0;

  if (! is_regular)
    {
      fprintf(logfile, "Waiting for data... (hit Ctrl-C and run again with '-h' for help)\n");
    }

  size_t linecap = linealloc;
  char * line {new char[linecap]};
  ssize_t linelen = xgetline(& line, & linecap, input_fp);
  if (linelen < 0)
    {
      line[0] = 0;
      linelen = 0;
    }
  filepos += linelen;

  unsigned int lineno {1};

  progress_init("Reading sequences:", filesize);

  while(line[0] != 0)
    {
      /* read header */
      /* the header ends at a space, cr, lf or null character */

      if (line[0] != '>') {
        fatal(error_prefix, "Illegal header line in fasta file.");
      }


      /* store the line number */

      while (datalen + sizeof(unsigned int) > dataalloc)
        {
          dataalloc += memchunk;
          d->datap = static_cast<char *>(xrealloc(d->datap, dataalloc));
        }
      memcpy(d->datap + datalen, & lineno, sizeof(unsigned int));
      datalen += sizeof(unsigned int);


      /* get next line */

      linelen = xgetline(& line, & linecap, input_fp);
      if (linelen < 0)
        {
          line[0] = 0;
          linelen = 0;
        }
      filepos += linelen;

      lineno++;


      /* store a dummy sequence length */

      unsigned int length {0};

      while (datalen + sizeof(unsigned int) > dataalloc)
        {
          dataalloc += memchunk;
          d->datap = static_cast<char *>(xrealloc(d->datap, dataalloc));
        }
      uint64_t datalen_seqlen = datalen;
      memcpy(d->datap + datalen, & length, sizeof(unsigned int));
      datalen += sizeof(unsigned int);


      /* read and store sequence */

      uint64_t nt_buffer {0};
      unsigned int nt_bufferlen {0};
      const unsigned int nt_buffersize {4 * sizeof(nt_buffer)};
      static constexpr int new_line {10};
      static constexpr int carriage_return {13};
      static constexpr int start_chars_range {32};  // visible ascii chars: 32-126
      static constexpr int end_chars_range {126};

      while ((line[0] != 0) && (line[0] != '>'))
        {
          unsigned char c {0};
          char * pl = line;
          while((c = static_cast<unsigned char>(*pl++)) != 0U)
            {
              signed char m {0};
              if ((m = map_nt[static_cast<unsigned int>(c)]) >= 0)
                {
                  nt_buffer |= ((static_cast<uint64_t>(m))) << (2 * nt_bufferlen);
                  length++;
                  nt_bufferlen++;

                  if (nt_bufferlen == nt_buffersize)
                    {
                      while (datalen + sizeof(nt_buffer) > dataalloc)
                        {
                          dataalloc += memchunk;
                          d->datap = static_cast<char *>(xrealloc(d->datap, dataalloc));
                        }

                      memcpy(d->datap + datalen, & nt_buffer, sizeof(nt_buffer));
                      datalen += sizeof(nt_buffer);

                      nt_bufferlen = 0;
                      nt_buffer = 0;
                    }
                }
              else if ((c != new_line) && (c != carriage_return))
                {
                  if ((c >= start_chars_range) && (c <= end_chars_range)) {
                    fatal(error_prefix, "Illegal character '", static_cast<char>(c),
                          "' in sequence on line ", lineno, ".");
                  }
                  else {
                    fatal(error_prefix, "Illegal character (ascii no ", static_cast<char>(c),
                          ") in sequence on line ", lineno, ".");
                  }
                }
            }

          linelen = xgetline(& line, & linecap, input_fp);
          if (linelen < 0)
            {
              line[0] = 0;
              linelen = 0;
            }
          filepos += linelen;

          lineno++;
        }

      /* fill in real length */

      memcpy(d->datap + datalen_seqlen, & length, sizeof(unsigned int));

      if (length == 0)
        {
          fatal(error_prefix, "Empty sequence found on line ", lineno - 1, ".");
        }

      d->nucleotides += length;

      if (length > d->longest) {
        d->longest = length;
      }


      /* save remaining padded 64-bit value with nt's, if any */

      if (nt_bufferlen > 0)
        {
          while (datalen + sizeof(nt_buffer) > dataalloc)
            {
              dataalloc += memchunk;
              d->datap = static_cast<char *>(xrealloc(d->datap, dataalloc));
            }

          memcpy(d->datap + datalen, & nt_buffer, sizeof(nt_buffer));
          datalen += sizeof(nt_buffer);

          nt_buffer = 0;
          nt_bufferlen = 0;  // that value is never read again, all tests pass without it
        }

      d->sequences++;

      if (is_regular) {
        progress_update(filepos);
      }
    }
  progress_done();

  fclose(input_fp);


  /* create indices */

  d->seqindex = (struct seqinfo_s *) xmalloc(d->sequences * sizeof(struct seqinfo_s));
  struct seqinfo_s * seqindex_p = d->seqindex;

  char * pl = d->datap;
  progress_init("Indexing database:", d->sequences);
  for(auto i = 0ULL; i < d->sequences; i++)
    {
      /* get line number */
      pl += sizeof(unsigned int);

      /* and sequence */
      unsigned int seqlen = *(reinterpret_cast<unsigned int*>(pl));
      seqindex_p->seqlen = seqlen;
      pl += sizeof(unsigned int);
      seqindex_p->seq = pl;
      pl += nt_bytelength(seqlen);

      seqindex_p++;
      progress_update(i);
    }
  progress_done();

  delete [] line;
  line = nullptr;
  linecap = 0;

  fprintf(logfile, "Database info:     %" PRIu64 " nt", d->nucleotides);
  fprintf(logfile, " in %u sequences,", d->sequences);
  fprintf(logfile, " longest %u nt\n", d->longest);

  return d;
}


void db_getsequenceandlength(struct db_s * d,
			     uint64_t seqno,
                             char ** address,
                             unsigned int * length)
{
  *address = d->seqindex[seqno].seq;
  *length = d->seqindex[seqno].seqlen;
}

uint64_t db_getnucleotides(struct db_s * d)
{
  return d->nucleotides;
}


void db_free(struct db_s * d)
{
  if (d->datap)
    xfree(d->datap);
  d->datap = nullptr;

  if (d->seqindex)
    xfree(d->seqindex);
  d->seqindex = nullptr;
}
