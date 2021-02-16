#include "util.h"
#include "io.h"
#include "alphabet.h"
#include "mergegap.h"


// prototype from blocks.h
int tba_get(uint64_t *a,customInt i);
void remap_bwts(g_data *g);


void array_clear(customInt* array, customInt size, customInt v) {
  for (customInt i = 0; i < size; ++i ) array[i] = v;
}


/**
 * read BWTs from a single file, the length of each input sequence is in the .size file  
 * alloc and init the fields in g:
 *   numBwt
 *   mergeLen
 *   bwtLen[]
 *   bws[]
 * 
 * Return false if the processing was not started since there is a single
 * multiBWT and there are not LCP or DB info values to compute (ie there
 * is nothing to do), otherwise return true   
 * */
bool readBWTsingle(char *path, g_data *g)
{
  char filename[Filename_size];
  customInt n=0;
  
  // open .size file and get number of (multi) BWT
  snprintf(filename,Filename_size,"%s.%s",path,LEN_EXT);
  FILE *f = fopen(filename,"r");
  if(f==NULL) die(__func__);
  if(fseek(f, 0, SEEK_END) != 0) {
    perror("fseek error");
    die(__func__);
  }
  off_t flen = ftello(f);
  rewind(f);
  if(flen < 0) {
    perror("ftello error");
    die(__func__);
  }
  else if(flen%8!=0) 
    die("Invalid format of ."LEN_EXT" file");
  g->numBwt = flen/8;
  assert(g->numBwt>0);

  // single BWT with no lcp or dbInfo computation: nothing to do
  if(g->numBwt==1 && !g->lcpCompute && g->dbOrder==0){ 
    if(g->outputSA) { // if there is a single (multi)-bwt the DA does not change
      char tmp1[Filename_size];
      char tmp2[Filename_size];
      snprintf(tmp1,Filename_size,"%s.%d.%s",path,g->outputSA, SA_BL_EXT);
      snprintf(tmp2,Filename_size,"%s.%d.%s",path,g->outputSA, SA_EXT);
      if(rename(tmp1,tmp2)!=0)
        die("Cannot rename SA file");
    }
    if(g->outputDA) { // if there is a single (multi)-bwt the DA does not change
      char tmp1[Filename_size];
      char tmp2[Filename_size];
      snprintf(tmp1,Filename_size,"%s.%d.%s",path,g->outputDA, DA_BL_EXT);
      snprintf(tmp2,Filename_size,"%s.%d.%s",path,g->outputDA, DA_EXT);
      if(rename(tmp1,tmp2)!=0)
        die("Cannot rename DA file");
    }
    if(g->outputQS) { // if there is a single (multi)-bwt the QS does not change
      char tmp1[Filename_size];
      char tmp2[Filename_size];
      snprintf(tmp1,Filename_size,"%s.%s",path,QS_BL_EXT);
      snprintf(tmp2,Filename_size,"%s.%s",path,QS_EXT);
      if(rename(tmp1,tmp2)!=0)
        die("Cannot rename QS file");
    }
    size_t size;
    size_t r = fread(&size, 8, 1, f);// sizes are stored in 64 bits 
    if(r!=1) die(__func__);
    g->mergeLen = size;
    fclose(f);
    return false;
  }

  if(g->outputDA) { 
    g->bwtDocs = (customInt *) malloc (g->numBwt * sizeof(customInt));
    snprintf(filename,Filename_size,"%s.%s",path,DOCS_EXT);
    FILE *f = fopen(filename,"r");
    if(f==NULL) die(__func__);
    size_t sum=0;
    for (int i = 0; i < g->numBwt; ++i) {
      customInt docs;
      size_t r = fread(&docs, 8, 1, f);// docs are stored in 64 bits 
      if(r!=1) die(__func__);
      g->bwtDocs[i] = sum;
      sum+=docs;
      //printf("--> %zu\n",  docs);
    }
    fclose(f);
    f = fopen(filename,"wb");
    size_t e = fwrite(&sum,8,1,f);
    if(e!=1) die(__func__);
    fclose(f);
  }

  // we have the number of (multi)BWTs: allocate g->bwtLen and g->bws  
  g->bwtLen = (customInt *) malloc (g->numBwt * sizeof(customInt));
  g->bws = (symbol **) malloc (g->numBwt * sizeof(symbol *));
  if(!g->bwtLen || !g->bws) die(__func__);
  // read and store sizes of multibwt in g->bwtLen 
  for (int i = 0; i < g->numBwt; ++i) {
    size_t size;
    size_t r = fread(&size, 8, 1, f);// sizes are stored in 64 bits 
    if(r!=1) die(__func__);
    if(size>MAX_OUTPUT_SIZE) {
      fprintf(stderr,"Bwt len %zu larger than maximum sequence size: %llu\n",
      size,MAX_OUTPUT_SIZE);
      die(__func__);
    }
    g->bwtLen[i] = size;
    n += size;
  }
  g->mergeLen = n;
  fclose(f);
  
  // We now have g->mergeLen: open and read concatenated bwt file 
  snprintf(g->bwfname,Filename_size,"%s.%s",path,BWT_EXT);
  int fd = open(g->bwfname,O_RDWR);
  if(fd == -1) die(__func__);
  // allocate/mmap memory for BWT concatenation and read whole file
  if(g->mmapBWT || g->extMem) { // mmap input file to RAM modifications are written back to the input file 
    g->bws[0] = mmap(NULL,n*sizeof(symbol),PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
    if(g->bws[0] == MAP_FAILED) die(__func__);
    #ifdef USE_MMAP_ADVISE
    // the first access is sequential to remap alphabet
    madvise(g->bws[0], g->mergeLen*sizeof(symbol), MADV_SEQUENTIAL);
    #endif
  }
  else { // allocate memory and read from file descriptor
    uint8_t *b = malloc(n*sizeof(symbol));
    if(!b) die(__func__);  
    huge_pread(fd,b,sizeof(symbol)*n,0);
    g->bws[0] = (symbol *) b;
  } 
  if(close(fd)!=0) die(__func__);

  // init g->bws[1..]: we do this even for g->extMem, since g->bws[i]-g->bws[0] is the starting point in the file of the i-th bwt
  for (int i = 0; i < g->numBwt -1; ++i)  
    g->bws[i+1] = g->bws[i] + g->bwtLen[i]; 
    
  // remap BWTs and also init g->sizeOfAlpha etc   
  remap_bwts(g);    
  if(g->extMem) { // if we are working in external memory don't keep the BWTs mmaped 
    int e = munmap(g->bws[0],g->mergeLen*sizeof(symbol));
    if(e) die("main (unmap bws)");
  }
  return true;
}

void array_copy(customInt *dest, customInt *src, customInt size){
  for (customInt i = 0; i < size; ++i) {
    dest[i] = src[i];
  }
}

// debug only funcion 
void check_g_data(g_data *g)
{
  for(int i=0;i<g->numBwt-1;i++) {
    assert(g->bws[i]+g->bwtLen[i]==g->bws[i+1]);
  }
  if(g->smallAlpha) 
    for(int i=0;i<g->numBwt;i++) {
      customInt tot=0;
      for(int j=0;j<g->sizeOfAlpha;j++) {
        tot += g->bwtOcc[i][j];
      }
      assert(tot==g->bwtLen[i]);
    }
  customInt mergeLen=0;
  for(int i=0;i<g->numBwt;i++) mergeLen += g->bwtLen[i];
  assert(mergeLen==g->mergeLen);
}


void open_unsortedLCP_files(g_data *g)
{
  char filename[Filename_size];
  snprintf(filename,Filename_size,"%s.pair.lcp",g->outPath);
  g->unsortedLcp = fopen(filename,"wb");
  if(g->unsortedLcp==NULL) {perror(filename); die(__func__);}
  snprintf(filename,Filename_size,"%s.size.lcp",g->outPath);
  g->unsortedLcp_size = fopen(filename,"wb");
  if(g->unsortedLcp_size==NULL) {perror(filename); die(__func__);}
  // write total size of the LCP array to first position of .size.lcp file 
  size_t e = fwrite(&(g->mergeLen),sizeof(customInt),1,g->unsortedLcp_size);  
  if(e!=1) die(__func__);
}

// write a lcp/position pair to g->unsortedLcp
// uses BSIZE bytes for the LCP value POS_SIZE bytes (currently 5) for the position
// everything is valid for little endian only! 
void writeLcp(customInt k, uint32_t lcp, g_data *g)
{
  assert(g->unsortedLcp!=NULL);
  assert(lcp <= MAX_LCP_SIZE);
  if(k>=(1ULL<<(8*POS_SIZE))) { // at most pos_size bytes for the position
    fprintf(stderr,"%d bytes per position are not enough. Increase POS_SIZE and recompile\n",POS_SIZE); 
    exit(EXIT_FAILURE);
  } 
  size_t e = fwrite(&lcp,BSIZE,1,g->unsortedLcp);
  if(e!=1) die(__func__);
  e = fwrite(&k,POS_SIZE,1,g->unsortedLcp);
  if(e!=1) die(__func__);
}

// write a lcp/position EOF to g->unsortedLcp, and the size of the current sorted block to g->unsortedLcp_size
void writeLcp_EOF(uint64_t size, g_data *g)
{
  assert(g->unsortedLcp!=NULL && g->unsortedLcp_size!=NULL);

  // write eof to .pair.lcp file
  uint64_t b = ~0ULL; // all 1's 
  size_t e = fwrite(&b,BSIZE,1,g->unsortedLcp);  
  if(e!=1) die(__func__);
  e = fwrite(&b,POS_SIZE,1,g->unsortedLcp);  
  if(e!=1) die(__func__);
  // write size of complete LCP segment to .size.lcp file 
  e = fwrite(&size,8,1,g->unsortedLcp_size);  
  if(e!=1) die(__func__);
}


void close_unsortedLCP_files(g_data *g)
{
   fclose(g->unsortedLcp); 
   fclose(g->unsortedLcp_size);
}



// allocate memory for the LCP values and mmap there the input LCPs
// only use when lcpMerge==true
// init g->lcps[i] = starting position of i-th input sequence LCP values
// g->lcps is initialized with a map-shared copy of the input LCP file
// so new values overwrite old ones and there is no need to copy them to the output file 
void initLCPmem(g_data *g)
{
  char filename[Filename_size];
  
  assert(g->lcpMerge);
  // allocate g->lcps and g->lcp[0..]
  g->lcps = (lcpInt **) malloc(g->numBwt*sizeof(lcpInt *));
  if(!g->lcps) die(__func__);
  
  snprintf(filename,Filename_size,"%s.%s",g->lcpinPath,LCP_EXT);
  int fd = open(filename,O_RDWR);
  if(fd == -1) die(__func__);
  // mmap lcp values, init lcps[] with a copy of the LCP values  
  g->lcps[0] = mmap(NULL,g->mergeLen*sizeof(lcpInt),PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
  if(g->lcps[0] == MAP_FAILED) die(__func__);
  #ifdef USE_MMAP_ADVISE
  madvise(g->lcps[0], g->mergeLen*sizeof(lcpInt), MADV_DONTNEED);// lcp values are not need untill the end of the bwt-merge phase 
  #endif
  if(close(fd)!=0) die(__func__);

  // init lcps[1..]
  for(int i=0;i<g->numBwt-1;i++)
    g->lcps[i+1] = g->lcps[i] + g->bwtLen[i];
}


static FILE *openDAFile(g_data *g)
{
  char filename[Filename_size];

  if(g->verbose>1) puts("Writing Document Array");
  assert(g->outputDA);
  snprintf(filename,Filename_size,"%s.%d.%s",g->outPath,g->outputDA,DA_EXT);
  FILE *f = fopen(filename,"wb");
  if(f==NULL) die("Error opening Document array file");
  return f;
}

static FILE *openSAFile(g_data *g)
{
  char filename[Filename_size];

  if(g->verbose>1) puts("Writing Suffix Array");
  assert(g->outputSA);
  snprintf(filename,Filename_size,"%s.%d.%s",g->outPath,g->outputSA,SA_EXT);
  FILE *f = fopen(filename,"wb");
  if(f==NULL) die("Error opening Suffix array file");
  return f;
}

static FILE *openQSFile(g_data *g)
{
  char filename[Filename_size];

  if(g->verbose>1) puts("Writing QS");
  assert(g->outputQS);
  snprintf(filename,Filename_size,"%s.%s",g->outPath,QS_EXT);
  FILE *f = fopen(filename,"wb");
  if(f==NULL) die("Error opening QS file");
  return f;
}

// use g->mergeColor to merge bwt (and LCP) values 
// used by gap gap16 and hm
void mergeBWTandLCP(g_data *g, bool lastRound)
{
  if(sizeof(symbol)>sizeof(palette)) die("Sorry, alphabet is too large (mergeBWTandLCP");
  assert(!g->lcpMerge || g->blockBeginsAt!=NULL); // if lcpMerge we need g->blockBeginsAt
  FILE *daOutFile=NULL;

  check_g_data(g);
  array_clear(g->inCnt,g->numBwt,0);
  if(g->outputDA && lastRound)
    daOutFile = openDAFile(g);
  if(g->extMem) {
    open_bw_files(g);// open files in read mode
    rewind_bw_files(g);   // set file pointers at the beginning of each BWT
    // mmap merge values
    int fd = open(g->merge_fname,O_RDWR);
    if(fd == -1) die(__func__);
    g->mergeColor = mmap(NULL,g->mergeLen*sizeof(palette),PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
    if(g->mergeColor == MAP_FAILED) die(__func__);
    #ifdef USE_MMAP_ADVISE
    madvise(g->mergeColor, g->mergeLen*sizeof(palette), MADV_SEQUENTIAL);
    #endif
    if(close(fd)!=0) die(__func__);
  }
  symbol *bwtout = (symbol *) g->mergeColor;   // merged BWT stored in mergeColor
  for (customInt i = 0; i < g->mergeLen; ++i) {
    int currentColor = g->mergeColor[i];
    assert(currentColor < g->numBwt);
    // check that all LCP values have been obtained
    if(g->lcpCompute && lastRound)
      assert( tba_get(g->bitB,i)==3 ); 
    // if requested output merge array
    if(g->outputDA && lastRound) 
      if(fputc(currentColor, daOutFile)==EOF)
        die("mergeBWTandLCP: Error writing to Document Array file");   
    // save new BWT char overwriting mergeColor[i]
    if(g->extMem) {
      int e = fread(bwtout+i,sizeof(symbol),1,g->bwf[currentColor]);
      if(e!=1) die(__func__);
    }      
    else bwtout[i] = g->bws[currentColor][g->inCnt[currentColor]];
    if(lastRound) bwtout[i] = alpha_enlarge(bwtout[i]); 
    
    if(g->lcpMerge) { // compute lcp value and store it to blockBeginsAt[]
      lcpInt tmp = g->lcps[currentColor][g->inCnt[currentColor]]; // lcp retrieved from file (suffixes from the same string)      
      if(g->blockBeginsAt[i]>0)
        g->blockBeginsAt[i] -= 1;  // lcp computed by the merging algorithm 
      else 
        g->blockBeginsAt[i] = tmp; // lcp retrieved from file (suffixes from the same string) 
    }
    g->inCnt[currentColor]++; // one more char read from currentColor BWT/LCP
  }  
  // final check on the merging 
  for(int i=0;i<g->numBwt;i++) assert(g->inCnt[i]==g->bwtLen[i]);
  if(g->extMem) close_bw_files(g);

  // merging done: copy back to g->bws[0] or to file
  if(g->extMem) {
    // printf("Copy back to file. Bytes: %ld, offset %ld\n",g->mergeLen, g->symb_offset);
    int fd = open(g->bwfname,O_WRONLY);
    if(fd == -1) die(__func__);
    huge_pwrite(fd, bwtout,sizeof(symbol)*g->mergeLen,sizeof(symbol)*g->symb_offset);
    if(close(fd)!=0) die(__func__);
    // unmap mergeColor
    fd = munmap(g->mergeColor,g->mergeLen*sizeof(palette));
    if(fd == -1) die(__func__);
    g->mergeColor=NULL;  
  }
  else 
    memcpy(g->bws[0],bwtout,g->mergeLen*sizeof(symbol));
  // close document array file 
  if(g->outputDA && lastRound)
      if(fclose(daOutFile)!=0) die("mergeBWTandLCP: Error closing Document Array file");   
  //copy merged values back to g->lcps[0]
  if(g->lcpMerge)
      memcpy(g->lcps[0],g->blockBeginsAt,g->mergeLen*sizeof(lcpInt)); // copy lcp values back to g->lcps[0]
}


// use colors in g->merge_fname to merge bwt values 
// the merged bwt is written back to file g->bwfname
// used only by gap128ext in merge128ext.h
void mergeBWT128ext(g_data *g, bool lastRound)
{
  if(sizeof(symbol)>sizeof(palette)) die("Sorry, alphabet is too large (mergeBWT128ext");
  assert(!g->lcpMerge && g->extMem && g->numBwt<=128);
  FILE *daOutFile=NULL;
  FILE *saOutFile=NULL;
  FILE *qsOutFile=NULL;
  check_g_data(g);
  array_clear(g->inCnt,g->numBwt,0);
  if(g->outputSA && lastRound)
    saOutFile = openSAFile(g);
  if(g->outputDA && lastRound)
    daOutFile = openDAFile(g);  
  if(g->outputQS && lastRound)
    qsOutFile = openQSFile(g);

  open_bw_files(g);     // open files in read mode
  rewind_bw_files(g);   // set file pointers at the beginning of each BWT

  if(g->outputSA && lastRound){
    snprintf(g->safname,Filename_size,"%s.%d.%s",g->outPath,g->outputSA, SA_BL_EXT);
    open_sa_files(g);
    rewind_sa_files(g);   // set file pointers at the beginning of each SA
  }
  if(g->outputDA && lastRound){
    snprintf(g->dafname,Filename_size,"%s.%d.%s",g->outPath,g->outputDA, DA_BL_EXT);
    open_da_files(g);
    rewind_da_files(g);   // set file pointers at the beginning of each DA
  }
  if(g->outputQS && lastRound){
    snprintf(g->qsfname,Filename_size,"%s.%s",g->outPath, QS_BL_EXT);
    open_qs_files(g);
    rewind_qs_files(g);   // set file pointers at the beginning of each DA
  }

  // mmap merge values so we can read/write simultaneously 
  int fd = open(g->merge_fname,O_RDWR);
  if(fd == -1) die(__func__);
  g->mergeColor = mmap(NULL,g->mergeLen*sizeof(palette),PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);
  if(g->mergeColor == MAP_FAILED) die(__func__);
  #ifdef USE_MMAP_ADVISE
  madvise(g->mergeColor, g->mergeLen*sizeof(palette), MADV_SEQUENTIAL);
  #endif
  if(close(fd)!=0) die(__func__);

  symbol *bwtout = (symbol *) g->mergeColor;   // merged BWT stored in mergeColor
  for (customInt i = 0; i < g->mergeLen; ++i) {
    int currentColor = g->mergeColor[i] & 0x7F; // delete additional bit
    assert(currentColor < g->numBwt);
    // if requested output merge array
    if(g->outputSA && lastRound){ 
      int sa_value=0;
      int e = fread(&sa_value, g->outputSA, 1, g->saf[currentColor]);
      if(e!=1) die(__func__);
      if(fwrite(&sa_value, g->outputSA, 1, saOutFile)==EOF)
        die("mergeBWT128ext: Error writing to Suffix Array file");   
    } 
    if(g->outputDA && lastRound){ 
      int da_value=0;
      int e = fread(&da_value, g->outputDA, 1, g->daf[currentColor]);
      if(e!=1) die(__func__);
      da_value+=g->bwtDocs[currentColor];
      //if(fputc(currentColor, daOutFile)==EOF)
      if(fwrite(&da_value, g->outputDA, 1, daOutFile)==EOF)
        die("mergeBWT128ext: Error writing to Document Array file");   
      //printf("%d ==> %d\n", currentColor, da_value);
    } 
    if(g->outputQS && lastRound){ 
      symbol qs_value=0;
      int e = fread(&qs_value, sizeof(symbol), 1, g->qsf[currentColor]);
      if(e!=1) die(__func__);
      if(fwrite(&qs_value, sizeof(symbol), 1, qsOutFile)==EOF)
        die("mergeBWT128ext: Error writing to QS file");   
    } 
    // save new BWT char overwriting mergeColor[i]
    int e = fread(bwtout+i,sizeof(symbol),1,g->bwf[currentColor]);
    if(e!=1) die(__func__);
    if(lastRound) bwtout[i] = alpha_enlarge(bwtout[i]);     
    g->inCnt[currentColor]++; // one more char read from currentColor BWT
  }  
  close_bw_files(g);
  if(g->outputSA  && lastRound) close_sa_files(g);
  if(g->outputDA  && lastRound) close_da_files(g);
  if(g->outputQS  && lastRound) close_qs_files(g);
  // final check on the merging 
  for(int i=0;i<g->numBwt;i++) assert(g->inCnt[i]==g->bwtLen[i]);

  // merging done: copy back to file
  fd = open(g->bwfname,O_WRONLY);
  if(fd == -1) die(__func__);
  huge_pwrite(fd, bwtout,sizeof(symbol)*g->mergeLen,sizeof(symbol)*g->symb_offset);
  if(close(fd)!=0) die(__func__);
  // unmap mergeColor
  fd = munmap(g->mergeColor,g->mergeLen*sizeof(palette));
  if(fd == -1) die(__func__);
  g->mergeColor=NULL;  
  // close suffix array file 
  if(g->outputSA && lastRound){
    if(fclose(saOutFile)!=0) die("mergeBWT128ext: Error closing Suffix Array file");   
    remove(g->safname);
  }
  // close document array file 
  if(g->outputDA && lastRound){
    if(fclose(daOutFile)!=0) die("mergeBWT128ext: Error closing Document Array file");   
    remove(g->dafname);
  }
  // close QS file 
  if(g->outputQS && lastRound){
    if(fclose(qsOutFile)!=0) die("mergeBWT128ext: Error closing QS file");   
    remove(g->qsfname);
  }
}


// merge BWT using colors stored in 7 bits of g->mergeColor16
// merged BWT is stored back to g->mergeColor
// used only by gap128 in merge128.h
void mergeBWT128(g_data *g, bool lastRound)
{
  if(sizeof(symbol)>2) die("Sorry, alphabet is too large (mergeBWT16)");
  symbol *bwtout = (symbol *) g->mergeColor16;   // merged BWT stored in mergeColor16
  assert(!g->lcpMerge);
  FILE *daOutFile=NULL;

  check_g_data(g);
  array_clear(g->inCnt,g->numBwt,0);
  if(g->outputDA && lastRound)
    daOutFile = openDAFile(g);  

  int shift = 7; (void) shift;  // shift is not used outside assertions 
  for (customInt i = 0; i < g->mergeLen; ++i) {
    int currentColor = g->mergeColor16[i] & 0x7F;
    assert(currentColor == ((g->mergeColor16[i]>>shift) & 0x7F)); // we must have mergeColor==newMergeColor
    assert(currentColor < g->numBwt);
    // check that all LCP values have been obtained
    if(g->lcpCompute && lastRound)
      assert( ((g->mergeColor16[i]>>2*shift) & 3)==3);
    // if requested output merge array
    if(g->outputDA && lastRound) 
      if(fputc(currentColor, daOutFile)==EOF)
        die("mergeBWT128: Error writing to Document Array file");   
    // save new BWT char overwriting mergeColor16[i]
    bwtout[i] = g->bws[currentColor][g->inCnt[currentColor]];
    g->inCnt[currentColor]++; // one more char read from currentColor BWT
  }
  // merging done: copy back to g->bws[0]
  if(lastRound)
    for(customInt i=0;i<g->mergeLen;i++)
      g->bws[0][i] = alpha_enlarge(bwtout[i]); // if last round remap while copying 
  else  
    memcpy(g->bws[0],bwtout,g->mergeLen*sizeof(symbol));  
  // final check
  for(int i=0;i<g->numBwt;i++) assert(g->inCnt[i]==g->bwtLen[i]);
  // close document array file 
  if(g->outputDA && lastRound)
      if(fclose(daOutFile)!=0) die("mergeBWT128: Error closing Document Array file");   
}


// merge BWT using colors stored in 3 bits of g->mergeColor
// merged BWT is stored first to g->mergeColor and then back to g->bws[0]
// used only by gap8() 
void mergeBWT8(g_data *g, bool lastRound)
{
  if(sizeof(symbol)>sizeof(palette)) die("Sorry, alphabet is too large (mergeBWT8)");
  symbol *bwtout = (symbol *) g->mergeColor;   // merged BWT stored in mergeColor
  assert(!g->lcpMerge);
  FILE *daOutFile=NULL;
  FILE *saOutFile=NULL;
  FILE *qsOutFile=NULL;

  check_g_data(g);
  array_clear(g->inCnt,g->numBwt,0); // clear counter inside each bwt
  if(g->outputSA && lastRound){
    saOutFile = openSAFile(g);    
    snprintf(g->safname,Filename_size,"%s.%d.%s",g->outPath,g->outputSA, SA_BL_EXT);
    open_sa_files(g);
    rewind_sa_files(g);   // set file pointers at the beginning of each SA
  }
  if(g->outputDA && lastRound){
    daOutFile = openDAFile(g);    
    snprintf(g->dafname,Filename_size,"%s.%d.%s",g->outPath,g->outputDA, DA_BL_EXT);
    open_da_files(g);
    rewind_da_files(g);   // set file pointers at the beginning of each DA
  }
  if(g->outputQS && lastRound){
    qsOutFile = openQSFile(g);    
    snprintf(g->qsfname,Filename_size,"%s.%s",g->outPath,QS_BL_EXT);
    open_qs_files(g);
    rewind_qs_files(g);   // set file pointers at the beginning of each QS
  }
  if(g->extMem) {
    open_bw_files(g);     // open files in read/write mode
    rewind_bw_files(g);   // set file pointers at the beginning of each BWT 
  }
  int shift = 3; (void) shift;  // shift is not used outside assertions
  for (customInt i = 0; i < g->mergeLen; ++i) {
    int currentColor = g->mergeColor[i] & 0x7;
    assert(currentColor == ((g->mergeColor[i]>>shift) & 0x7)); // we must have mergeColor==newMergeColor
    assert(currentColor < g->numBwt);
    // check that all LCP values have been obtained
    if(g->lcpCompute && lastRound)
      assert( ((g->mergeColor[i]>>2*shift) & 3)==3);    
    // if requested output merge array
    if(g->outputSA && lastRound){ 
      int sa_value=0;
      int e = fread(&sa_value, g->outputSA, 1, g->saf[currentColor]);
      if(e!=1) die(__func__);
      if(fwrite(&sa_value, g->outputSA, 1, saOutFile)==EOF)
        die("mergeBWT128ext: Error writing to Sufix Array file");   
    } 
    if(g->outputDA && lastRound){ 
      int da_value=0;
      int e = fread(&da_value, g->outputDA, 1, g->daf[currentColor]);
      if(e!=1) die(__func__);
      //if(fputc(currentColor, daOutFile)==EOF)
      da_value+=g->bwtDocs[currentColor];
      if(fwrite(&da_value, g->outputDA, 1, daOutFile)==EOF)
        die("mergeBWT128ext: Error writing to Document Array file");   
      //printf("%d ==> %d\n", currentColor, da_value);
    } 
    if(g->outputQS && lastRound){ 
      symbol qs_value=0;
      int e = fread(&qs_value, sizeof(symbol), 1, g->qsf[currentColor]);
      if(e!=1) die(__func__);
      if(fwrite(&qs_value, sizeof(symbol), 1, qsOutFile)==EOF)
        die("mergeBWT128ext: Error writing to QS file");   
    } 
    // save new BWT char overwriting mergeColor[i]
    if(g->extMem) {
      int e = fread(bwtout+i,sizeof(symbol),1,g->bwf[currentColor]);
      if(e!=1) die(__func__);
    }      
    else bwtout[i] = g->bws[currentColor][g->inCnt[currentColor]];
    if(lastRound) bwtout[i] = alpha_enlarge(bwtout[i]); 
    g->inCnt[currentColor]++; // one more char read from currentColor BWT
  }
  // final check on the merging 
  for(int i=0;i<g->numBwt;i++) assert(g->inCnt[i]==g->bwtLen[i]);
  if(g->extMem) close_bw_files(g);
  // close suffix array file 
  if(g->outputSA && lastRound){
    if(fclose(saOutFile)!=0) die("mergeBWT8: Error closing Suffix Array file");       
    remove(g->safname);
  }
  // close document array file 
  if(g->outputDA && lastRound){
    if(fclose(daOutFile)!=0) die("mergeBWT8: Error closing Document Array file");       
    remove(g->dafname);
  }
  // close suffix array file 
  if(g->outputQS && lastRound){
    if(fclose(qsOutFile)!=0) die("mergeBWT8: Error closing QS file");       
    remove(g->qsfname);
  }

  // merging done: copy back to g->bws[0] or to file
  if(g->extMem) {
    if(g->verbose>1) printf("Copy back to file. Bytes: %zu, offset %zu\n",g->mergeLen, g->symb_offset);
    int fd = open(g->bwfname,O_WRONLY);
    if(fd == -1) die(__func__);
    huge_pwrite(fd, bwtout,sizeof(symbol)*g->mergeLen,sizeof(symbol)*g->symb_offset);
    if(close(fd)!=0) die(__func__);
  }
  else 
    memcpy(g->bws[0],bwtout,g->mergeLen*sizeof(symbol));
}


// merge BWTs and if required compute LCPs using B and Z stored in g->array32
// the merged BWT is stored back to g->bws
// LCP values if computed are written to the outPath.BSIZE.lcp file 
// Note: we are not doing a lcpMerge so we are not using g->lcps (this is why we write directly to the LCP file)  
// Note: if it is not the last round the LCP values are not used
// used only by gap256()
void mergeBWTandLCP256(g_data *g, bool lastRound)
{
  if(sizeof(symbol)>sizeof(g->array32[0])) die("Sorry, alphabet is too large (mergeBWTandLCP256)");
  symbol *bwtout = (symbol *) g->array32;   // merged BWT stored in g->array32
  assert(!g->lcpMerge);
  FILE *daOutFile=NULL;
  FILE *lcpfile = NULL;
  int shift = 8;

  // if last round open file for LCP values if requested
  if(lastRound && g->lcpCompute) {
    char filename[Filename_size];
    snprintf(filename,Filename_size,"%s.%s",g->outPath,LCP_EXT);
    lcpfile = fopen(filename,"wb");
    if(lcpfile==NULL) {perror("Unable to open LCP output file"); die(__func__);}
  }
  
  check_g_data(g);
  array_clear(g->inCnt,g->numBwt,0);
  if(g->outputDA && lastRound)
    daOutFile = openDAFile(g);  
  for (customInt i = 0; i < g->mergeLen; ++i) {
    int currentColor = (g->array32[i]) & 0xFF;
    assert(currentColor == ((g->array32[i]>>shift) & 0xFF));// merge and newMerge should be identical 
    assert(currentColor < g->numBwt);
    // if requested output merge array
    if(g->outputDA && lastRound) 
      if(fputc(currentColor, daOutFile)==EOF)
        die("mergeBWTandLCP256: Error writing to Document Array file");       
    // get next BWT char
    symbol c = g->bws[currentColor][g->inCnt[currentColor]];
    // write LCP to file if requested
    if(lcpfile!=NULL) {
      uint64_t curB = (g->array32[i]>>(2*shift)) & 0xFFFF;
      if(curB==0) die("Illegal value in B-array");
      curB -= 1; // obtain true LCP
      if(curB > MAX_LCP_SIZE) die("LCP value too large for requested lcp size");  
      // only valid for little endian architectures!!!!
      size_t e = fwrite(&curB,BSIZE,1,lcpfile);
      if(e!=1) {perror("Error writing LCP values"); die(__func__);}
    }
    // save new BWT char overwriting array32[i]
    bwtout[i] =  c; 
    g->inCnt[currentColor]++; // one more char read from currentColor BWT
  }
  // merging done
  if(lcpfile!=NULL) fclose(lcpfile);
  // copy the merged bwt back to bws[0]
  if(lastRound)
    for(customInt i=0;i<g->mergeLen;i++)
      g->bws[0][i] = alpha_enlarge(bwtout[i]); // if last round remap while copying 
  else  
    memcpy(g->bws[0],bwtout,g->mergeLen*sizeof(symbol));
  // final check
  for(int i=0;i<g->numBwt;i++) assert(g->inCnt[i]==g->bwtLen[i]);
  // close document array file 
  if(g->outputDA && lastRound)
      if(fclose(daOutFile)!=0) die("mergeBWTandLCP256: Error closing Document Array file");   

}


// allocate or mmap arrays Z and newZ
// only used by mergegap and mergehm (the latter does not support extermnal memory)
void alloc_merge_arrays(g_data *g) {
  if(g->extMem) { // notice extMem overrules mmap
    int e,fd;
    e = asprintf(&(g->merge_fname),"%s.mrg0_XXXXXX",g->outPath);
    if(e<0) die("merge temp file name creation failed (1)");
    fd = mkstemp(g->merge_fname);
    if(fd == -1) die("merge temp file creation failed (1)");
    if(close(fd)!=0) die("merge temp file close failed (1)");  // we don't need it now
    e = asprintf(&(g->newmerge_fname),"%s.mrg1_XXXXXX",g->outPath);
    if(e<0) die("merge temp file name creation failed (2)");
    fd = mkstemp(g->newmerge_fname);
    if(fd == -1) die("merge temp file creation failed (2)");
    if(close(fd)!=0) die("merge temp file close failed (2)");  // we don't need it now
    g->mergeColor = g->newMergeColor = NULL;
  }
  else if(g->mmapZ) {
    g->mergeColor = mmap(NULL,g->mergeLen*sizeof(palette),PROT_READ|PROT_WRITE,MAP_PRIVATE|MAP_ANONYMOUS,-1,0);
    if(g->mergeColor== MAP_FAILED) die(__func__);
    g->newMergeColor =  mmap(NULL,g->mergeLen*sizeof(palette),PROT_READ|PROT_WRITE,MAP_PRIVATE|MAP_ANONYMOUS,-1,0);
    if(g->newMergeColor== MAP_FAILED) die(__func__);
  }
  else {
    g->mergeColor =     malloc(g->mergeLen*sizeof(palette)); 
    g->newMergeColor =  malloc(g->mergeLen*sizeof(palette));
     if(!g->mergeColor || !g->newMergeColor) die(__func__);
  }
  // make sure these are not used 
  g->mergeColor16=NULL;
  g->array32 = NULL;
}

void free_merge_arrays(g_data *g) {
  if(g->extMem) {
    int e = unlink(g->merge_fname);
    if(e!=0) perror("Error deleting temporary merge file (1)");
    e = unlink(g->newmerge_fname);
    if(e!=0) perror("Error deleting temporary merge file (2)");
    free(g->merge_fname); // deallocate file names 
    free(g->newmerge_fname);
  }
  else if(g->mmapZ) {
    int e = munmap(g->mergeColor,g->mergeLen*sizeof(palette));
    if(e!=0) die("gap (unmap mergeColor)");
    e = munmap(g->newMergeColor,g->mergeLen*sizeof(palette));
    if(e!=0) die("gap (unmap newMergeColor)");  }
  else {
    free(g->newMergeColor);
    free(g->mergeColor);
  }
  g->newMergeColor = NULL;
  g->mergeColor = NULL;
}


// allocate or mmap array Z and not newZ (since it shares the same space as Z)
// usied inmerge16.h and merge8.h
void alloc_merge_array(g_data *g) {
  if(g->mmapZ) {
    g->mergeColor = mmap(NULL,g->mergeLen*sizeof(palette),PROT_READ|PROT_WRITE,MAP_PRIVATE|MAP_ANONYMOUS,-1,0);
    if(g->mergeColor== MAP_FAILED) die(__func__);
  }
  else {
    g->mergeColor =  calloc(g->mergeLen,sizeof(palette)); 
     if(!g->mergeColor) die(__func__);
  }
  // make sure these are not used 
  g->newMergeColor = NULL;
  g->mergeColor16=NULL;
  g->array32 = NULL;
  g->fmergeColor = NULL;
  g->fnewMergeColor = NULL;
}

void free_merge_array(g_data *g) {
  assert(g->newMergeColor==NULL);
  if(g->mmapZ) {
    int e = munmap(g->mergeColor,g->mergeLen*sizeof(palette));
    if(e!=0) die("gap (unmap mergeColor)");
  }
  else {
    free(g->mergeColor);
  }
  g->mergeColor = NULL;
}

// allocate or mmap array B: zero intialized (we use this) 
void alloc0_B_array(g_data *g) {
  if(g->mmapB) {
    // zero intialized 
    g->blockBeginsAt = mmap(NULL,g->mergeLen*sizeof(lcpInt),PROT_READ|PROT_WRITE,MAP_PRIVATE|MAP_ANONYMOUS,-1,0);
    if(g->blockBeginsAt== MAP_FAILED) die(__func__);
  }
  else {
    g->blockBeginsAt = (lcpInt *) calloc(g->mergeLen,sizeof(lcpInt));
    if(!g->blockBeginsAt) die(__func__);
  }
}

void free_B_array(g_data *g) {
  if(g->mmapB) {
    int e = munmap(g->blockBeginsAt,g->mergeLen*sizeof(lcpInt));
    if(e) die(__func__);
  }
  else free(g->blockBeginsAt); // if used free storage for lcp values
  g->blockBeginsAt = NULL;
}


// same as above for 16 bit merge array used in merge128.h 
// allocate or mmap array Z and newZ
// array is initialized to 0 
void alloc_merge_array16(g_data *g) {
  if(g->mmapZ) {
    g->mergeColor16 = mmap(NULL,g->mergeLen*sizeof(uint16_t),PROT_READ|PROT_WRITE,MAP_PRIVATE|MAP_ANONYMOUS,-1,0);
    if(g->mergeColor16 == MAP_FAILED) die(__func__);
  }
  else {
    g->mergeColor16 =  calloc(g->mergeLen,sizeof(uint16_t)); 
     if(!g->mergeColor16) die(__func__);
  }
  g->mergeColor = g->newMergeColor = NULL;
  g->array32 = NULL;
  g->fmergeColor = NULL;
  g->fnewMergeColor = NULL;
}

void free_merge_array16(g_data *g) {
  assert(g->mergeColor==NULL && g->newMergeColor==NULL && g->array32==NULL);
  if(g->mmapZ) {
    int e = munmap(g->mergeColor16,g->mergeLen*sizeof(uint16_t));
    if(e!=0) die("gap (unmap mergeColor16)");
  }
  else {
    free(g->mergeColor16);
  }
  g->mergeColor16 = NULL; // make sure it is no longer used
}

// same as above for 32-bit array 
// allocate or mmap array Z and not newZ
// array is initialized to 0 
void alloc_array32(g_data *g) {
  if(g->mmapZ) {
    g->array32 = mmap(NULL,g->mergeLen*sizeof(uint32_t),PROT_READ|PROT_WRITE,MAP_PRIVATE|MAP_ANONYMOUS,-1,0);
    if(g->array32 == MAP_FAILED) die(__func__);
  }
  else {
    g->array32 =  calloc(g->mergeLen,sizeof(uint32_t)); 
    if(!g->array32) die(__func__);
  }
  g->mergeColor = g->newMergeColor = NULL;
  g->mergeColor16 = NULL; // make sure it is not used
}

void free_array32(g_data *g) {
  assert(g->mergeColor==NULL && g->newMergeColor==NULL && g->mergeColor16 == NULL);
  if(g->mmapZ) {
    int e = munmap(g->array32,g->mergeLen*sizeof(uint32_t));
    if(e!=0) die("gap (unmap array32)");
  }
  else {
    free(g->array32);
  }
  g->array32 = NULL; // make sure it is no longer used
}

 
// error message and exit
void die(const char* where)  {
  printf("Error at %s: %s.\n",where,errno ? strerror(errno) : "errno not set");
  exit(errno?errno:1);
}

