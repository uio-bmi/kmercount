# Makefile for kmercount

PROG=bin/kmercount

kmercount : $(PROG)

$(PROG) :
	make -C src kmercount

install : $(PROG) $(MAN)
	/usr/bin/install -c $(PROG) '/usr/local/bin'

clean :
	make -C src clean
