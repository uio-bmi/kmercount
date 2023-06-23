#!/bin/sh

if ! [ -e ../src/kmercount ] ; then
    echo The kmercount binary is missing
    echo Test failed.
    exit 1
fi

../src/kmercount -k 31 kmers.fasta seq.fasta -l kmercount.log -o counts.tsv

if diff -q counts.tsv expected.tsv; then
    echo Test completed successfully.
else
    echo Test failed.
    exit 1
fi
