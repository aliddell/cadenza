CC=gcc
MPI=`which mpicc`
LIB= -lm -lmpfr -lgmp
OPT= -O3 -funroll-loops -fexpensive-optimizations -Wall 
ARGS=$(OPT) 
#OBJFILES=misc.o eval.o norm.o LUdecomp.o newton.o classify.o classify_over.o loadSettings.o output.o refine.o sqrt.o newtonOnly.o 
OBJFILES=io.o

#all : $(OBJFILES) $(POBJFILES) blueharvest ;
all : blueharvest ;

#blueharvest : blueharvest.c $(OBJFILES) ;
blueharvest : blueharvest.c ;
	#$(CC) $(ARGS) -o blueharvest blueharvest.c $(OBJFILES) $(LIB)
	$(MPI) $(ARGS) -o blueharvest blueharvest.c $(LIB)

.c.o : 
	$(CC) $(ARGS) -c $*.c 

clean : ;
	#rm -f blueharvest $(OBJFILES)
	rm -f blueharvest
