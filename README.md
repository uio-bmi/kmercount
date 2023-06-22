# Kmercounter

Kmercounter (`kmercounter`) is a command line tool to count the number of occurences of a set of selected kmers in a set of sequences. A k-mer is defined as a continuous sequence of `k` nucleotides. Kmercounter will work with kmers of up to 32 nucleotides.

Kmercounter uses an efficient hash function, Bloom filter and hash table to perform the counting rapidly.


## Compilation and installation

The code is C++11 standard compliant and should compile easily using
`make` and a modern C++ compiler (e.g. GNU GCC or LLVM Clang). Run
`make clean`, `make`, `make install` in the main
folder to clean, build, and install the tool. There are no
dependencies except for the C and C++ standard libraries.


## General options

Here are the Kmercount options:

```
Kmercount 0.0.1

Usage: kmercount [OPTIONS] KMERFILENAME [SEQUENCEFILENAME]

General options:
 -h, --help                 display this help and exit
 -k, --kmer-length INTEGER  kmer length [1-32] (31)
 -t, --threads INTEGER      number of threads to use [1-256] (1)
 -v, --version              display version information and exit

Input/output options:
 -l, --log FILENAME         log to file (stderr)
 -o, --output FILENAME      output result to file (stdout)
```

The input file with the selected kmers must be specified as the first
positional argument. It must be a FASTA formatted file with sequences
of the correct length (k). The FASTA headers are ignored.

The input file with the sequences to scan for kmers may be specified
as the second postional argument. If not specified, or specified as
`-`, the program will read from standard input. The input must be in
FASTA format. The headers are ignored.

Use the `-h` or `--help` option to show some help information.

The kmer length may be specified with the `-k` or `--kmerlength`
option. The length must be in the range from 1 to 32. The default kmer
length is 31.

The number of parallel threads requested may be specified with the
`-t` or `--threads` option. However, the code is not multi-threaded,
yet.

Run the program with `-v` or `--version` for version information.

While the program is running it will print some status and progress
information to standard error (stderr) unless a log file has been
specified with the `-l` or `--log` option. Error messages and warnings
will also be written here.

The results will be written to standard output (stdout) unless a file
name has been specified with the `-o` or `--output` option.


## General information

All input sequences must be nucleotide sequences. The letters
`ACGTNacgtn` are accepted, but be aware that `N` and `n` are treated
as `A`.

The kmer sequences should be distinct. If they are not distinct, a
warning will be given.


## Example

Command line:

```
kmercount --kmer-length 31 --output counts.tsv kmers50m.fasta human30M.fasta
```

Output to terminal:

```
Kmercount 0.0.1

Kmer file:         kmers50m.fasta
Sequence file:     human30M.fasta
Kmer length:       31
Output file:       counts.tsv
Threads:           1

Reading kmer file
Reading sequences: 100%  
Indexing database: 100%  
Database info:     1550000000 nt in 50000000 sequences, longest 31 nt
Indexing kmers:    100%  
Unique kmers:      50000000

Reading sequence file
Reading sequences: 100%  
Indexing database: 100%  
Database info:     3030000000 nt in 30000000 sequences, longest 101 nt
Counting matches:  100%  

Writing results:   100%  
Matching kmers:    93468
Total matches:     1903965
```

Output to `counts.tsv` (only first 10 lines shown):

```
AAATGGGCTAAATGCTCCAATTAAAAGACAC	45
GTAATCCCAGCACTTTGGAAGGCCGAGGCAG	15
GACAGAGAGCCAAATCATGAGTGAACTCCCA	15
AGATCACTTGAGGCCAGGAGTTCAAGACCAG	15
TTGGTGCAGAGCTGAGTTCAATTCCTGGGTA	15
AAATACAAAAATTAGCCAGGTGTGGTGGCAC	60
CACATCACACTTATTCCAAAATTGACCACAT	15
TGAGGTTTTGGGACTCAGACTGGCTTCCTTG	15
GGATCTCACTCTGTTGCCCAGGCTAGAGTGC	15
CTTGAGGTGTCAGTGGCAGATAGAGATGCTG	15
```
