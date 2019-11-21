MF=	Makefile

#CC=	mpicc -cc=icc
CC= mpicc -cc=gcc
CFLAGS=	-O3
LIB=	-lm

LFLAGS= $(CFLAGS)

EXE=	percolate

INC= \
	percolate.h \
	pinit_data.h    \
	pdistri_pro.h   \
	pupdate_squares.h   \
	pcollect_data.h \
	test_perc.h

SRC= \
	percolate.c \
	percwrite.c \
	uni.c   \
	arralloc.c  \
	pinit_data.c    \
	pdistri_pro.c   \
	pupdate_squares.c   \
	pcollect_data.c \
	test_perc.c

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .c .o

OBJ=	$(SRC:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -c $<

percolate_para:	$(EXE)

$(OBJ):	$(INC)

$(EXE):	$(OBJ)
	$(CC) $(LFLAGS) -o $@ $(OBJ) $(LIB)

$(OBJ):	$(MF)

clean:
	rm -f $(EXE) $(OBJ) core map.pgm
