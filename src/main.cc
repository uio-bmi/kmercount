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


/* OPTIONS */

std::string progname;  // unused variable?

struct Parameters p;
std::string opt_log;
int64_t opt_threads;

/* fine names and command line options */

std::FILE * outfile {nullptr};
std::FILE * logfile {stderr};  // cstdio stderr macro is expanded to type std::FILE*

constexpr int n_options {26};
std::array<int, n_options> used_options {{0}};  // set int values to zero by default

char short_options[] = "hl:o:t:v"; /* unused: abcdefgijkmnpqrsuwxyz*/

static struct option long_options[] =
  {
   {"help",                  no_argument,       nullptr, 'h' },
   {"log",                   required_argument, nullptr, 'l' },
   {"output-file",           required_argument, nullptr, 'o' },
   {"threads",               required_argument, nullptr, 't' },
   {"version",               no_argument,       nullptr, 'v' },
   {nullptr,                 0,                 nullptr, 0 }
  };


const std::vector<std::string> header_message
  {"kmercount ", program_version,
   "\n"};


const std::vector<std::string> args_usage_message
  /*0         1         2         3         4         5         6         7          */
  /*01234567890123456789012345678901234567890123456789012345678901234567890123456789 */
  {"Usage: kmercount [OPTIONS] KMERFILE [FASTAFILE]\n",
   "\n",
   "General options:\n",
   " -h, --help                          display this help and exit\n",
   " -t, --threads INTEGER               number of threads to use (1)\n",
   " -v, --version                       display version information and exit\n",
   "\n",
   "Input/output options:\n",
   " -l, --log FILENAME                  log to file, not to stderr\n",
   " -o, --output-file FILENAME          output result to file (stdout)\n",
   "\n",
#ifndef __WIN32
   "\n",
   "See 'man kmercount' for more details.\n",
#endif
   "\n"
  };

auto args_long(char * str, const char * option) -> int64_t;
void args_show();
void show(const std::vector<std::string> & message);
void args_init(int argc, char **argv, std::array<int, n_options> & used_options);
void args_check(std::array<int, n_options> & used_options);
void open_files();
void close_files();


auto args_long(char * str, const char * option) -> int64_t
{
  static constexpr int base_value {10};
  char * endptr {nullptr};
  const int64_t temp = strtol(str, & endptr, base_value);
  if (*endptr != 0)
    {
      fatal(error_prefix, "Invalid numeric argument for option ", option, ".\n\n",
            "Frequent causes are:\n",
            " - a missing space between an argument and the next option,\n",
            " - a long option name not starting with a double dash\n",
            "   (the program accepts '--help' or '-h', but not '-help')\n\n",
            "Please run again with '--help' for more details.");
    }
  return temp;
}


void args_show()
{
#ifdef __x86_64__
  cpu_features_detect(p);
  cpu_features_test(p);
  cpu_features_show(p);
#endif

  fprintf(logfile, "Kmer file:         %s\n", p.kmer_filename.c_str());
  fprintf(logfile, "Sequence file:     %s\n", p.seq_filename.c_str());
  fprintf(logfile, "Output file:       %s\n", p.opt_output_file.c_str());
  fprintf(logfile, "Threads:           %" PRId64 "\n", opt_threads);
  fprintf(logfile, "\n");
}


void show(const std::vector<std::string> & message)
{
  for (const auto & m : message) {
    fprintf(logfile, "%s", m.c_str());
  }
}


void args_init(int argc, char **argv, std::array<int, n_options> & used_options)
{
  static constexpr unsigned int threads_default {1};

  progname = argv[0];
  opt_threads = threads_default;
  opterr = 1;  // unused variable? get_opt option?

  int c {0};

  while (true)
  {
    int option_index {0};
    c = getopt_long(argc, argv, short_options, long_options, &option_index);

    if (c == -1) {
      break;
    }

    /* check if any option is specified more than once */

    if ((c >= 'a') && (c <= 'z'))
      {
        auto optindex = static_cast<unsigned int>(c - 'a');  // c - 'a' cannot be negative
        if (used_options[optindex] == 1)
          {
            int longoptindex {0};
            while (long_options[longoptindex].name != nullptr)
              {
                if (long_options[longoptindex].val == c) {
                  break;
                }
                longoptindex++;
              }
            fatal(error_prefix, "Option -", static_cast<char>(c),
                  " or --", long_options[longoptindex].name,
                  " specified more than once.");
          }
        used_options[optindex] = 1;
      }

    switch(c)
      {
      case 'h':
        /* help */
        p.opt_help = true;
        break;

      case 'l':
        /* log */
        opt_log = optarg;
        break;

      case 'o':
        /* output-file */
        p.opt_output_file = optarg;
        break;

      case 't':
        /* threads */
        opt_threads = args_long(optarg, "-t or --threads");
        break;

      case 'v':
        /* version */
        p.opt_version = true;
        break;

      default:
        show(header_message);
        show(args_usage_message);
        fatal();
    }
  }

  if (optind < argc)
    {
      if (optind + 1 < argc)
	{
	  p.kmer_filename = argv[optind];
	  p.seq_filename = argv[optind + 1];
	}
      else
	{
	  p.kmer_filename = argv[optind];
	}
    }
  else
    fatal("At least one filename must be specified (kmer file)");
}


void args_check(std::array<int, n_options> & used_options) {
  static constexpr unsigned int max_threads {256};
  // meaning of the used_options values

  (void) used_options;

  if ((opt_threads < 1) || (opt_threads > max_threads))
    {
      fatal(error_prefix, "Illegal number of threads specified with "
            "-t or --threads, must be in the range 1 to ", max_threads, ".");
    }

  if (p.opt_version) {
    show(header_message);
    std::exit(EXIT_SUCCESS);
  }

  if (p.opt_help) {
    show(header_message);
    show(args_usage_message);
    std::exit(EXIT_SUCCESS);
  }
}


void open_files()
{
  // special case (always '-')??
  outfile = fopen_output(p.opt_output_file.c_str());
  if (outfile == nullptr) {
    fatal(error_prefix, "Unable to open output file for writing.");
  }

  /* open files */

  if (! opt_log.empty())
    {
      logfile = fopen_output(opt_log.c_str());
      if (logfile == nullptr) {
        fatal(error_prefix, "Unable to open log file for writing.");
      }
    }
}


auto close_files() -> void {
  const std::vector<std::FILE *> file_handles
    {outfile, logfile};
  for (auto * const file_handle : file_handles) {
    if (file_handle != nullptr) {
      fclose(file_handle);
    }
  }
}


auto main(int argc, char** argv) -> int
{
  args_init(argc, argv, used_options);
  args_check(used_options);
  open_files();
  show(header_message);
  args_show();
  kmercount(p.kmer_filename.c_str(), p.seq_filename.c_str());
  close_files();
}
