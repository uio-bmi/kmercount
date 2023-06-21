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
constexpr long unsigned int n_chars {INT8_MAX + 1};  // 128 ascii chars
constexpr unsigned int max_sequence_length {999999999};
// for longer sequences, 'zobrist_tab_byte_base' is bigger than 8 x
// 2^32 (512 x max_sequence_length) and cannot be addressed with
// uint32 pointers, which leads to a segmentation fault
constexpr unsigned int max_header_length {16777216 - 1};  // 2^24 minus 1

auto make_nt_map () -> std::array<signed char, n_chars> {
    // set the 128 ascii chars to '-1' except Aa, Cc, Gg, Tt and Uu
  std::array<signed char, n_chars> ascii_map {{0}};
    ascii_map.fill(-1);
    ascii_map['A'] = 1;
    ascii_map['a'] = 1;
    ascii_map['C'] = 2;
    ascii_map['c'] = 2;
    ascii_map['G'] = 3;
    ascii_map['g'] = 3;
    ascii_map['T'] = 4;
    ascii_map['t'] = 4;
    ascii_map['U'] = 4;
    ascii_map['u'] = 4;
    ascii_map['N'] = 1;
    ascii_map['n'] = 1;
    return ascii_map;
    }

const auto map_nt = make_nt_map();

const std::array<char, 32> sym_nt =
  {'-', 'A', 'C', 'G', 'T', ' ', ' ', ' ',
   ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
   ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
   ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};


static unsigned int sequences {0};
static uint64_t nucleotides {0};
static uint64_t headerchars {0};
static unsigned int longest {0};
static unsigned int longestheader {0};
static char * datap {nullptr};

struct seqinfo_s
{
  char * header;
  char * seq;
  uint64_t abundance;
  uint64_t hdrhash;
  uint64_t seqhash;
  int headerlen;
  unsigned int seqlen;
  unsigned int clusterid;
  int abundance_start;
  int abundance_end;
  int dummy; /* alignment padding only */
};

using seqinfo_t = struct seqinfo_s;
extern seqinfo_t * seqindex;

seqinfo_t * seqindex {nullptr};


auto db_getnucleotidecount() -> uint64_t
{
  return nucleotides;
}


auto db_getsequencecount() -> unsigned int
{
  return sequences;
}


auto db_getlongestsequence() -> unsigned int
{
  return longest;
}



void db_read(const char * filename, struct Parameters const & p)
{
  (void) p;

  /* allocate space */

  uint64_t dataalloc {memchunk};
  datap = static_cast<char *>(xmalloc(dataalloc));
  uint64_t datalen {0};

  longest = 0;
  longestheader = 0;
  sequences = 0;
  nucleotides = 0;
  headerchars = 0;

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

      auto headerlen = static_cast<unsigned int>
        (strcspn(line + 1, " \r\n"));

      headerchars += headerlen;

      if (headerlen > longestheader) {
        longestheader = headerlen;
      }

      if (longestheader > max_header_length) {
        fatal(error_prefix, "Headers longer than 16,777,215 symbols are not supported.");
      }

      /* store the line number */

      while (datalen + sizeof(unsigned int) > dataalloc)
        {
          dataalloc += memchunk;
          datap = static_cast<char *>(xrealloc(datap, dataalloc));
        }
      memcpy(datap + datalen, & lineno, sizeof(unsigned int));
      datalen += sizeof(unsigned int);


      /* store the header */

      while (datalen + headerlen + 1 > dataalloc)
        {
          dataalloc += memchunk;
          datap = static_cast<char *>(xrealloc(datap, dataalloc));
        }
      memcpy(datap + datalen, line + 1, headerlen);
      *(datap + datalen + headerlen) = 0;
      datalen += headerlen + 1;


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
          datap = static_cast<char *>(xrealloc(datap, dataalloc));
        }
      uint64_t datalen_seqlen = datalen;
      memcpy(datap + datalen, & length, sizeof(unsigned int));
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
                  nt_buffer |= ((static_cast<uint64_t>(m))-1) << (2 * nt_bufferlen);
                  length++;
                  nt_bufferlen++;

                  if (nt_bufferlen == nt_buffersize)
                    {
                      while (datalen + sizeof(nt_buffer) > dataalloc)
                        {
                          dataalloc += memchunk;
                          datap = static_cast<char *>(xrealloc(datap, dataalloc));
                        }

                      memcpy(datap + datalen, & nt_buffer, sizeof(nt_buffer));
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

          /* check length of longest sequence */
          if (length > max_sequence_length) {
            fatal(error_prefix, "Sequences longer than 67,108,861 symbols are not supported.");
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

      memcpy(datap + datalen_seqlen, & length, sizeof(unsigned int));

      if (length == 0)
        {
          fatal(error_prefix, "Empty sequence found on line ", lineno - 1, ".");
        }

      nucleotides += length;

      if (length > longest) {
        longest = length;
      }


      /* save remaining padded 64-bit value with nt's, if any */

      if (nt_bufferlen > 0)
        {
          while (datalen + sizeof(nt_buffer) > dataalloc)
            {
              dataalloc += memchunk;
              datap = static_cast<char *>(xrealloc(datap, dataalloc));
            }

          memcpy(datap + datalen, & nt_buffer, sizeof(nt_buffer));
          datalen += sizeof(nt_buffer);

          nt_buffer = 0;
          nt_bufferlen = 0;  // that value is never read again, all tests pass without it
        }

      sequences++;

      if (is_regular) {
        progress_update(filepos);
      }
    }
  progress_done();

  fclose(input_fp);



  /* create indices */

  seqindex = new seqinfo_t[sequences];
  seqinfo_t * seqindex_p {seqindex};



  char * pl {datap};
  progress_init("Indexing database:", sequences);
  for(auto i = 0ULL; i < sequences; i++)
    {
      /* get line number */
      pl += sizeof(unsigned int);

      /* get header */
      seqindex_p->header = pl;
      seqindex_p->headerlen = static_cast<int>(strlen(seqindex_p->header));
      pl += seqindex_p->headerlen + 1;

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

  // user report
  fprintf(logfile, "Database info:     %" PRIu64 " nt", db_getnucleotidecount());
  fprintf(logfile, " in %u sequences,", db_getsequencecount());
  fprintf(logfile, " longest %u nt\n", db_getlongestsequence());
}


auto db_gethash(uint64_t seqno) -> uint64_t
{
  return seqindex[seqno].seqhash;
}


auto db_getsequence(uint64_t seqno) -> char *
{
  return seqindex[seqno].seq;
}


void db_getsequenceandlength(uint64_t seqno,
                             char ** address,
                             unsigned int * length)
{
  *address = seqindex[seqno].seq;
  *length = seqindex[seqno].seqlen;
}


auto db_getsequencelen(uint64_t seqno) -> unsigned int
{
  return seqindex[seqno].seqlen;
}


auto db_getheader(uint64_t seqno) -> char *
{
  return seqindex[seqno].header;
}


auto db_getabundance(uint64_t seqno) -> uint64_t
{
  return seqindex[seqno].abundance;
}


void db_free()
{
  if (datap != nullptr) {
    xfree(datap);
  }
  datap = nullptr;
  delete [] seqindex;
  seqindex = nullptr;
}


auto db_fprintseq(std::FILE * fastaout_fp, const unsigned int seqno) -> void
{
  const unsigned int len {db_getsequencelen(seqno)};
  char * const seqptr {db_getsequence(seqno)};
  static std::vector<char> buffer(db_getlongestsequence() + 1, '\0');

  // decode to nucleotides (A, C, G and T)
  for(auto i = 0U; i < len; i++) {
    buffer[i] = sym_nt[1 + nt_extract(seqptr, i)];
  }
  buffer[len] = '\0';

  fprintf(fastaout_fp, "%.*s\n", len, buffer.data());
}
