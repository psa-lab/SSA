#########################################################################
#   Tharuna Niranjan							#
#   Kuhn Laboratory							#
#   Dept. of Biochemistry						#
#   Michigan State University						#
#									#
#	Makefile for SSA automation code		5.8.96		#
#########################################################################

CC = cc 

#Add -DDEBUG to CFLAGS for displaying debugging statements (optional)
CFLAGS = -g 
 
LIBS = -lm

install:ssa
	rm -f core *.o
	
all:ssa
	chmod 755 ssa
	@ echo "ssa uptodate"

SSA_FILES=ssa.o templatematch.o rotlsqfit.o rotate.o matvec.o eigen.o emsg.o

SRC_FILES=ssa.c templatematch.c rotlsqfit.c rotate.c matvec.c eigen.c emsg.c

ssa: $(SSA_FILES)
	$(CC) $(CFLAGS) -o ssa $(SSA_FILES) $(LIBS)

clean:
	rm -f core *.o
	@echo "Deleted the OBJ files"

ssa.o: ssa.c templatematch.c ssa.h
	$(CC) -c $(CFLAGS) ssa.c

templatematch.o: templatematch.c ssa.h
	$(CC) -c $(CFLAGS) templatematch.c

rotlsqfit.o: rotlsqfit.c rotate.c matvec.c graph.h mymacs.h
	$(CC) -c $(CFLAGS) rotlsqfit.c

rotate.o: rotate.c eigen.c
	$(CC) -c $(CFLAGS) rotate.c

matvec.o:matvec.c 
	$(CC) -c $(CFLAGS) matvec.c

eigen.o:eigen.c emsg.c
	$(CC) -c $(CFLAGS) eigen.c

emsg.o:emsg.c
	$(CC) -c $(CFLAGS) emsg.c
