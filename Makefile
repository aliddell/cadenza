CC=gcc
MPI=`which mpicc`
LIB= -lm -lmpfr -lgmp
OPT= -g -O3 -funroll-loops -fexpensive-optimizations -Wall 
ARGS=$(OPT) 
OBJFILES=alphaCertified.o certify_float.o certify_rational.o classify.o classify_over.o eval.o io.o io_float.o io_rational.o loadSettings.o LUdecomp.o misc.o newton.o newtonOnly.o norm.o output.o refine.o sqrt.o

all : $(OBJFILES) $(POBJFILES) blueharvest ;

blueharvest : blueharvest.c $(OBJFILES) ;
	$(MPI) $(ARGS) -o blueharvest blueharvest.c $(OBJFILES) $(LIB)

.c.o : 
	$(MPI) $(ARGS) -c $*.c 

clean : ;
	rm -f blueharvest $(OBJFILES)
