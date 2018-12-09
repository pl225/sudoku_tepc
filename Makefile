EXEFILE      = sudoku
VERSION      = -D_VERSION=\"BENCHMARK_0\"
APPLICATION  = -D_APPLICATION=\"$(EXEFILE)\"
CPUCC     = mpicc
CPUCPP    = g++
DEFS      = $(APPLICATION)  $(VERSION)

INCLUDES  =	-I. -I/usr/local/share/papi/testlib/

LIBDIR   = -L/usr/local/lib/
         
LIBS     =  -lm  -lpapi

LINK     =  $(LIBDIR) $(LIBS) 

DEFS      += 
CPPFLAGS  +=  -O3 -std=c11 -march=native -Wall -fopenmp -lpthread


all:	 main
	$(CPUCC)	$(DEFS) $(INCLUDES) $(CPPFLAGS) sudoku.o $(LINK) \
	/usr/local/share/papi/testlib/libtestlib.a \
	/usr/local/lib/libpapi.a \
	-o $(EXEFILE) 

main:
	$(CPUCC) $(DEFS) $(INCLUDES) $(CPPFLAGS) -c sudoku.c
	
clean:
	rm *.o; rm $(EXEFILE)

files:
	rm *.txt; rm *.dat; rm *.o