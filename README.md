# Kmercount

Kmercount (`kmercount`) is a simple command line tool to quickly count
the number of occurences of a set of selected kmers in a set of
sequences. A k-mer is defined as a continuous sequence of `k`
nucleotides. Kmercount will work with kmers of up to 32 nucleotides.

Kmercount uses an efficient hash function, Bloom filter and hash
table to perform the counting rapidly.


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
 -t, --threads INTEGER      number of threads to use [1-1] (1)
 -v, --version              display version information and exit

Input/output options:
 -l, --log FILENAME         log to file (stderr)
 -o, --output FILENAME      output result to file (stdout)
```

Use the `-h` or `--help` option to show some help information.

Run the program with `-v` or `--version` for version information.

The input file with the selected kmers must be specified as the first
positional argument. It must be a FASTA formatted file with sequences
of the correct length (k). The FASTA headers are ignored.

The input file with the sequences to scan for kmers may be specified
as the second postional argument. If not specified, or specified as
`-`, the program will read from standard input. The input must be in
FASTA format. The headers are ignored.

The kmer length may be specified with the `-k` or `--kmerlength`
option. The length must be in the range from 1 to 32. The default kmer
length is 31.

The number of parallel threads requested may be specified with the
`-t` or `--threads` option. However, the code is not multi-threaded,
yet.

While the program is running it will print some status and progress
information to standard error (stderr) unless a log file has been
specified with the `-l` or `--log` option. Error messages and warnings
will also be written here.

The results will be written to standard output (stdout) unless a file
name has been specified with the `-o` or `--output` option. The output
is a plain text file with tab-separated values. The first column
contains the kmer sequences, while the second column contains the
counts. The kmers are sorted by descending number of occurences.


## General information

All input sequences must be nucleotide sequences. The letters
`ACGTNacgtn` are accepted, but be aware that `N` and `n` are treated
as `A`.

The kmer sequences should be distinct. The number of distinct (unique)
kmers will be shown.


## Example

Command line:

```
kmercount -k 31 -o counts.tsv kmers.fasta seq.fasta
```

The `kmers.fasta` file:

```
>kmer1
AAGAAATGAGAAGTAATCAGAAAACCACTTA
>kmer2
AGAAATGAGAAGTAATCAGAAAACCACTTAA
>kmer3
GAAATGAGAAGTAATCAGAAAACCACTTAAG
>kmer4
AAATGAGAAGTAATCAGAAAACCACTTAAGG
```


The `seq.fasta` file:

```
>seq
AAGAAATGAGAAGTAATCAGAAAACCACTTAAGG
```


Output to terminal:

```
Kmercount 0.0.1

Kmer file:         kmers.fasta
Sequence file:     seq.fasta
Kmer length:       31
Output file:       counts.tsv
Threads:           1

Reading kmer file
Reading sequences: 100%  
Indexing database: 100% 
Database info:     124 nt in 4 sequences, longest 31 nt
Indexing kmers:    100% 
Unique kmers:      4

Reading sequence file
Reading sequences: 100%  
Indexing database: 100%
Database info:     34 nt in 1 sequences, longest 34 nt
Counting matches:  100%  

Sorting results:   100%
Writing results:   100% 
Matching kmers:    4
Total matches:     4
```

Output to `counts.tsv` (only first 10 lines shown):

```
AAATGAGAAGTAATCAGAAAACCACTTAAGG	1
AGAAATGAGAAGTAATCAGAAAACCACTTAA	1
AAGAAATGAGAAGTAATCAGAAAACCACTTA	1
GAAATGAGAAGTAATCAGAAAACCACTTAAG	1
```
