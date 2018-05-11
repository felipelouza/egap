CC = gcc
CFLAGS += -Wall 
CFLAGS += -D_FILE_OFFSET_BITS=64 -m64 -O3 -fomit-frame-pointer -Wno-char-subscripts 

LFLAGS = -lm -lrt -ldl
DEBUG	= 0

DEFINES = -DDEBUG=$(DEBUG) 
CFLAGS += $(DEFINES)

##

# object files for gsaca
LIBOBJ = \
	lib/utils.o\
	lib/file.o\
	lib/suffix_array.o\
	lib/lcp_array.o\
	src/gsacak.o\
	src/malloc_count.o


# object files for gsaca-64
LIBOBJ64 = \
	lib/utils-64.o\
	lib/file-64.o\
	lib/suffix_array-64.o\
	lib/lcp_array-64.o\
	src/gsacak-64.o\
	src/malloc_count.o


# object files for mergelcp
MERGEOBJ = \
	lib/utils.o\
	src/malloc_count.o\
	heap.o

EXECS = gsacak gsacak-64 mergelcp

all: ${EXECS}

# rule for 64 bits object files 
%-64.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $< -DM64

gsacak: main_gsacak.c ${LIBOBJ} 
	$(CC) $(CFLAGS) -o $@ $< ${LIBOBJ} $(LFLAGS) -DM64=0

gsacak-64: main_gsacak.c ${LIBOBJ64} 
	$(CC) $(CFLAGS) -o $@ $< ${LIBOBJ64} $(LFLAGS) -DM64=1

heap.o: heap.c heap.h
	$(CC) $(CFLAGS) -c -o $@ $<

mergelcp: mergelcp.c ${MERGEOBJ}
	$(CC) $(CFLAGS) -o $@ $< ${MERGEOBJ} $(LFLAGS)

clean:
	\rm -f *.o lib/*.o src/*.o ${EXECS}
