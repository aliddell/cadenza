CC=gcc
MPI=`which mpicc`
LIB= -lm -lmpfr -lgmp
OPT= -g -O3 -funroll-loops -fexpensive-optimizations -Wall 
ARGS=$(OPT) 
OBJFILES=io.o

all : $(OBJFILES) $(POBJFILES) blueharvest ;

blueharvest : blueharvest.c $(OBJFILES) ;
	$(MPI) $(ARGS) -o blueharvest blueharvest.c $(OBJFILES) $(LIB)

.c.o : 
	$(MPI) $(ARGS) -c $*.c 

clean : ;
	rm -f blueharvest $(OBJFILES)
