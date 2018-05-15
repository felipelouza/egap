CC = gcc
CFLAGS = -D_FILE_OFFSET_BITS=64 -O -g -Wall -Winline -std=c99 -fms-extensions -pthread
CFFLAGS = -D_FILE_OFFSET_BITS=64 -O3 -g -Wall -Winline -std=c99 -fms-extensions -pthread

# -ldl is required by malloc_count 
LFLAGS = -lm -ldl 

HEADERS = *.h

CFILES = gap.c util.c io.c mergegap.c mergehm.c alphabet.c malloc_count/malloc_count.c threads.c multiround.c
CFILES0 = gap.c util.c io.c mergegap.c mergehm.c alphabet.c malloc_count/malloc_count0.c threads.c multiround.c


EXECS = gap gap1 gap2 gap4

# targets not producing a file declared phony
.PHONY: all tools clean tarfile

all: $(EXECS) tools

# BWTs/LCPs merging (no malloc_count: had some conflicts with -O)
gap:  $(CFILES0) $(HEADERS)
	$(CC) $(CFLAGS) $(CFILES0) -lm -DBSIZE=2 -ogap

# fast nodebug versions of the gap algorithm 
gap1:  $(CFILES) $(HEADERS)
	$(CC) $(CFFLAGS) $(CFILES) $(LFLAGS) -DNDEBUG -DBSIZE=1 -ogap1

gap2:  $(CFILES) $(HEADERS)
	$(CC) $(CFFLAGS) $(CFILES) $(LFLAGS) -DNDEBUG -DBSIZE=2 -ogap2

gap4:  $(CFILES) $(HEADERS)
	$(CC) $(CFFLAGS) $(CFILES) $(LFLAGS) -DNDEBUG -DBSIZE=4 -ogap4

# executables in tools directory for phases 1 and 3
tools:
	make -C tools

tarfile:
	tar -zcf egap.tgz readme.txt eGap Makefile *.[ch] malloc_count/*.[ch]\
           tools/*.[ch] tools/Makefile tools/*/*.[ch]

clean:
	\rm -f $(EXECS)
	make clean -C tools

####

INPUT=dataset/pacbio.fasta	#input data
RAM=10		#memory in MB 
K=0 			#number of strings

run:
	./tools/gsacak-64 $(INPUT) $(K) -bm $(RAM)
	./gap2 -vEla -A128 -g128 $(INPUT) 
	./tools/mergelcp $(INPUT) 5 2

####

