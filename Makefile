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
    perclib.h    \
	mplib.h

SEQINC= \
    percolate.h \
    perclib.h

PARSRC= \
	percolate.c \
	percwrite.c \
	uni.c   \
	arralloc.c  \
	par_init_data.c    \
	par_distri_pro.c   \
	par_update_squares.c   \
	par_collect_data.c \
	test_perc.c \
	mplib.c

SEQSRC= \
    percolate.c \
    seq_init_data.c    \
    seq_distri_pro.c   \
    seq_update_squares.c   \
    seq_collect_data.c \
    percwrite.c \
    uni.c   \
    arralloc.c  \
    test_perc.c


#
# Sequential version
#

.SUFFIXES:
.SUFFIXES: .c .o

SEQOBJ=	$(SEQSRC:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -c $<

percolate_s:	$(SEQEXE)

$(SEQOBJ):	$(SEQINC)

$(SEQEXE):	$(SEQOBJ)
	$(CC) $(LFLAGS) -o $@ $(SEQOBJ) $(LIB)

$(SEQOBJ):	$(MF)

clean_seq:
	rm -f $(SEQEXE) $(SEQOBJ) core map.pgm

#
# Parallel version
#

.SUFFIXES:
.SUFFIXES: .c .o

PAROBJ=	$(PARSRC:.c=.o)

percolate_p:	$(PAREXE)

$(PAROBJ):	$(PARINC)

$(PAREXE):	$(PAROBJ)
	$(CC) $(LFLAGS) -o $@ $(PAROBJ) $(LIB)

$(PAROBJ):	$(MF)

clean_par:
	rm -f $(PAREXE) $(PAROBJ) core map.pgm
