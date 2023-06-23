# Makefile for Kmercount

ifndef PREFIX
	PREFIX=/usr/local
endif

all : kmercount

kmercount:
	make -C src kmercount

test: kmercount
	make -C test

install: kmercount test
	/usr/bin/install -d $(PREFIX)/bin
	/usr/bin/install -c src/kmercount $(PREFIX)/bin/kmercount

clean:
	make -C src clean
	make -C test clean
	rm -f *~
