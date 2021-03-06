CC=/usr/bin/gcc
MPI=/usr/bin/mpicc
LIB= -lm -lgmp -lmpfr
OPT= -g -O3 -funroll-loops -fexpensive-optimizations -Wall 
ARGS=$(OPT) 
OBJFILES=alphaCertified.o certify_float.o certify_rational.o classify.o classify_over.o eval.o io.o io_float.o io_rational.o loadSettings.o LUdecomp.o misc.o newton.o newtonOnly.o norm.o output.o refine.o sqrt.o

all : $(OBJFILES) $(POBJFILES) cadenza ;

cadenza : cadenza.c $(OBJFILES) ;
	$(CC) $(ARGS) -o cadenza cadenza.c $(OBJFILES) $(LIB)

.c.o : 
	$(CC) $(ARGS) -c $*.c 

clean : 
	rm -f cadenza $(OBJFILES)
reset :
	rm -f approxSolns constantValues distinctSolns isApproxSoln isDistinctSoln isRealSoln nonrealDistinctSolns realDistinctSolns redundantSolns refinedPoints unknownPoints summary
