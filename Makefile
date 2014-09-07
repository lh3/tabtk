CC=gcc
PROG=tabtk
CFLAGS=-g -Wall -O2 -Wno-unused-function

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

tabtk:tabtk.o kseq.h kstring.h
		$(CC) $(CFLAGS) tabtk.c -o $@ -lz -lm

clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) -- *.c)

# DO NOT DELETE THIS LINE -- make depend depends on it.

regexp9.o: regexp9.h
tabtk.o: kvec.h kstring.h ksort.h kseq.h
