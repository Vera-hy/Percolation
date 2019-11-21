MF=	Makefile

#CC=	mpicc -cc=icc
CC= mpicc -cc=gcc
CFLAGS=	-O3
LIB=	-lm

LFLAGS= $(CFLAGS)

PAREXE=	percolate_par

SEQEXE= percolate_seq

PARINC= \
	percolate.h \
	pinit_data.h    \
	pdistri_pro.h   \
	pupdate_squares.h   \
	pcollect_data.h \
	test_perc.h

SEQINC= \
    percolate.h

PARSRC= \
	percolate.c \
	percwrite.c \
	uni.c   \
	arralloc.c  \
	pinit_data.c    \
	pdistri_pro.c   \
	pupdate_squares.c   \
	pcollect_data.c \
	test_perc.c

SEQSRC= \
    percolate.c

#
# Sequential program
#

.SUFFIXES:
.SUFFIXES: .c .o

SEQOBJ=	$(SEQSRC:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -c $<

percolate_seq:	$(SEQEXE)

$(SEQOBJ):	$(SEQINC)

$(SEQEXE):	$(SEQOBJ)
	$(CC) $(LFLAGS) -o $@ $(SEQOBJ) $(LIB)

$(SEQOBJ):	$(MF)

clean_seq:
	rm -f $(SEQEXE) $(SEQOBJ) core map.pgm

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .c .o

PAROBJ=	$(PARSRC:.c=.o)

#.c.o:
#	$(CC) $(CFLAGS) -c $<

percolate_p:	$(PAREXE)

$(PAROBJ):	$(PARINC)

$(PAREXE):	$(PAROBJ)
	$(CC) $(LFLAGS) -o $@ $(PAROBJ) $(LIB)

$(PAROBJ):	$(MF)

clean:
	rm -f $(PAREXE) $(PAROBJ) core map.pgm
