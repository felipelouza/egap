// collection of functions for accessing data in external memory 
#include "util.h"
#include "io.h"

// open a file for each input BWT, file pointers are stored to bwf[] 
void open_bw_files(g_data *g) {
  assert(g->extMem);
  g->bwf = malloc(g->numBwt*sizeof(FILE *));
  if(g->bwf==NULL) die(__func__);
  for(int i=0; i< g->numBwt; i++) {
    g->bwf[i] = fopen(g->bwfname,"r");
    if(!g->bwf[i]) die(__func__);
  }
}

// use bws[] to make bwf[i] point at the beginning of bws[i]  
void rewind_bw_files(g_data *g) {
  assert(g->extMem);
  for(int i=0; i< g->numBwt; i++) {
    int e = fseek(g->bwf[i],sizeof(symbol)*(g->bws[i]-g->bws[0]+g->symb_offset),SEEK_SET);
    if(e!=0) die(__func__);
  }
}

// use bws[] to make bwf[i] point at the beginning of bws[i]  
void close_bw_files(g_data *g) {
  assert(g->extMem);
  for(int i=0; i< g->numBwt; i++) {
    int e = fclose(g->bwf[i]);
    if(e!=0) die(__func__);
  }
  free(g->bwf);
}

// creation of temporary files for irrelevant blocks
// the file is not visible since it is deleted after creation
FILE *gap_tmpfile(char* path)
{
  // create local copy of template
  char s[strlen(path)+11];
  sprintf(s,"%s.tmpXXXXXX",path);
  assert(strlen(s)==strlen(path) + 10);
  // get file descriptor for tmp file 
  int fd = mkstemp(s);
  if(fd == -1) die("gap_tmpfile: Tempfile creation failed (1)");
  // get the FILE * (we need buffering)
  FILE *f = fdopen(fd,"w+"); 
  if(f==NULL)  die("gap_tmpfile: Tempfile creation failed (2)");
  // unlink file so it is deleted as soon as it is closed
  int e = unlink(s);
  if(e!=0)     die("gap_tmpfile: Tempfile creation failed (3)");
  return f;
}


// read/write huge blocks to file using multiple pread/pwrite calls 

void huge_pwrite(int fd, const void *vbuf, size_t count, off_t offset)
{
  const char *buf = (const char *) vbuf;
  while(count>0) {
    ssize_t w = pwrite(fd,buf,count,offset);
    if(w==0) die(__func__);
    count -= w;
    buf += w;
    offset += w;
  }
}
    
void huge_pread(int fd, void *vbuf, size_t count, off_t offset)
{
  char *buf = (char *) vbuf;
  while(count>0) {
    ssize_t w = pread(fd,buf,count,offset);
    if(w==0) die(__func__);
    count -= w;
    buf += w;
    offset += w;
  }
}


// ----- color writer functions 
static void cwriter_flush(cwriter *w)
{
  if(w->cur>0) {
    huge_pwrite(w->fd,w->buffer,w->cur*sizeof(palette),w->offset);
    w->offset += w->cur*sizeof(palette);
    w->cur=0;
  }
}

void cwriter_put(cwriter *w, palette c)
{
  if(w->cur==w->size) cwriter_flush(w);
  assert(w->cur < w->size);
  w->buffer[w->cur++] = c;
} 

void cwriter_skip(cwriter *w, uint64_t s) {
  cwriter_flush(w);
  w->offset += s*sizeof(palette); 
}

void cwriter_close(cwriter *w) {
  cwriter_flush(w);
  free(w->buffer);
}

off_t cwriter_tell(cwriter *w) {
  return w->offset + (w->cur*sizeof(palette));
}

void cwriter_init(cwriter *w, int fd, size_t size, off_t o) {
  assert(size>0);
  w->buffer = malloc(size*sizeof(palette));
  if(!w->buffer) die(__func__);
  w->size = size;
  w->cur = 0;
  w->offset = o;
  w->fd = fd;
}

// --- bit file structure and related functions ---

// save b->cur bits to b->fd. correspondingly advance b->offset
static void bitfile_save(bitfile *b, bool endfile)
{
  if(b->cur>0) {
    if(b->cur%8!=0 && !endfile) die("Illegal bitfile save");
    size_t bytes = (b->cur+7)/8;
    assert(bytes*8 <= b->size);
    huge_pwrite(b->fd,b->buffer,bytes,b->offset);
    b->offset += bytes;
    b->cur=b->size=0;
  }
}

// save to file the bits currently in buffer
// used only at the end of a reading/writing cycle 
void bitfile_flush(bitfile *b)
{
  bitfile_save(b, true);
  assert(b->offset==(b->filesize+7)/8); 
}

// init a bitfile: opening an anonymous file and filling it with zeros 
void bitfile_create(bitfile *b, size_t size, char *path) {
  // create local copy of template
  char s[strlen(path)+11];
  sprintf(s,"%s.bitXXXXXX",path);
  assert(strlen(s)==strlen(path) + 10);
  // get file descriptor for tmp file fill it zith 0s and delete file  
  b->fd = mkstemp(s);
  if(b->fd == -1) die("bitfile_create: Tempfile creation failed");
  int e = ftruncate(b->fd,(size+7)/8); // fill with size 0 bits
  if(e!=0)        die("bitfile_create: Tempfile ftruncate failed");
  e = unlink(s);
  if(e!=0)        die("bitfile_create: Tempfile unlink failed");
  b->filesize = size;  // total size of the bitfile (in bits) 
  b->buffer = malloc(Bitfile_bufsize_bytes);
  if(!b->buffer)  die("bitfile_create: malloc error");
  b->size = 0;         // current size of the buffer in bits
  b->cur = 0;          // bit index inside buffer
  b->offset = 0;       // offset in bytes in the file 
}

// destroy a bitfile closing the corresponding file and freeing the buffer 
void bitfile_destroy(bitfile *b) {
  free(b->buffer);
  b->buffer= NULL;
  if(close(b->fd)!=0) die(__func__);
}

// set virtual pointer at the beginning of the file 
void bitfile_rewind(bitfile *b)
{
  b->size = 0;         // current size of the buffer in bits
  b->cur = 0;          // bit index inside buffer
  b->offset = 0;       // offset in bytes in the file   
}

// read the next bit b from bitfile
// change it to b|new and return b 
bool bitfile_read_or_write(bitfile *b, bool new)
{
  // make sure there is a bit to read in the buffer, possibily reading from disk 
  if(b->cur >= b->size) {
    assert(b->cur==b->size);
    bitfile_save(b,false);
    assert(b->cur==0);
    assert(b->size==0);
    ssize_t n = pread(b->fd,b->buffer,Bitfile_bufsize_bytes,b->offset);
    if(n<=0) die("Unable to read bitfile data (bitfile_read_or_write)");
    b->size = n*8; // number of bits available for reading
    // note we did not change offset since we are going to write at this offset 
  }
  assert(b->cur<b->size);
  int i = b->cur/8;
  int j = b->cur%8;
  b->cur++;
  bool old = b->buffer[i] & (1<<j); // get old bit value
  if(new) b->buffer[i] |= (1<<j);
  return old;
}

// skip an assigned number of bits 
void bitfile_skip(bitfile *b, uint64_t s) {
  if(b->cur+s <= b->size) { // easy case: we stay in the current buffer
    b->cur+= s;
    return;
  }
  // complete current byte
  int delta = (8 - (b->cur%8)) % 8;
  s -= delta;
  b->cur += delta;
  assert(b->cur%8==0 && b->cur <= b->size&& s>0);
  bitfile_save(b,false); // advance offset to the current b->cur
  // skip as many full bytes as possible
  b->offset += s/8;       
  s = s%8;                // we are left with < 8 bits
  ssize_t n = pread(b->fd,b->buffer,Bitfile_bufsize_bytes,b->offset); // read a full buffer 
  if(n<0 ) die("Error reading bitfile data (bitfile_skip)");
  else if(n==0 && s>0) die("Unable to read bitfile data (bitfile_skip)");
  b->size = n*8; // number of bits available for reading
  b->cur = s;    // virtually skip the remaining bits 
}
 
// return index of the next bit to be read
off_t bitfile_tell(bitfile *b) {
  return 8*b->offset + b->cur;
}
