#ifndef IO_H_INCLUDED
#define IO_H_INCLUDED

#include "config.h"


// used to access BWTs in external memory 
void open_bw_files(g_data *g);
void rewind_bw_files(g_data *g);
void close_bw_files(g_data *g);

FILE *gap_tmpfile(char* path);
void huge_pwrite(int fd, const void *buf, size_t count, off_t offset);
void huge_pread(int fd, void *buf, size_t count, off_t offset);
void cwriter_put(cwriter *w, palette c);
void cwriter_skip(cwriter *w, uint64_t s);
void cwriter_close(cwriter *w);
off_t cwriter_tell(cwriter *w);
void cwriter_init(cwriter *w, int fd, size_t size, off_t o);


// by default a bitfile buffer is 8 file buffers
#define Bitfile_bufsize_bytes (8*BUFSIZ)

// structure supporting sequential read&write on a file representing a sequence of bits
typedef struct{
  int fd;       // file descriptor
  off_t offset; // offset inside file descriptor (in bytes)
  uint8_t *buffer; 
  size_t filesize;  // max number of bits stored in file
  int size;         // actual buffer size  (in bits)
  int cur;          // current position (in bits)
} bitfile;

void bitfile_flush(bitfile *b);
void bitfile_create(bitfile *b, size_t size, char *path);
void bitfile_destroy(bitfile *b);
void bitfile_rewind(bitfile *b);
bool bitfile_read_or_write(bitfile *b, bool new);
void bitfile_skip(bitfile *b, uint64_t s);
off_t bitfile_tell(bitfile *b);








#endif
