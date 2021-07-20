#ifndef CONFIG_H
#define CONFIG_H
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <time.h>
#include <pthread.h>
#include <semaphore.h>
#include <limits.h>
#include <inttypes.h>
#ifdef __linux__
#include <linux/limits.h>
#endif

// comment to prevent the use of madvise
#define USE_MMAP_ADVISE 1

// size of buffer for building file names
#ifndef PATH_MAX
#define PATH_MAX 4096
#endif
#define Filename_size PATH_MAX

// size of all buffers for producer/consumer threads
#define Threads_buf_size 20
// size of each buffer for external memory newMerge array
#define COLOR_WBUFFER_SIZE (1024*1024)


// type used to represent an input symbol
typedef unsigned char symbol;
#define SIZE_OF_ALPHABET 256
#define BWT_EOF 0x00
#define ILLEGAL_SYMBOL (-1)

//type used to represent occurrences in small solid blocks
typedef uint16_t smallSolidInt;
#define SMALLSOLID_FORMAT "%"PRIu16
#define SMALLSOLID_LIMIT UINT16_MAX

// type used to represent positions in the BWT's 
// define the maximum (merged) BWT size 
typedef uint64_t customInt;
#define CUSTOM_FORMAT "%"PRIu64
#define MAX_OUTPUT_SIZE 0xFFFFFFFFFFFFFFFEULL

// number of bytes used to represent a position in the BWT's for irrelevant blocks and mergeLcp 
#define POS_SIZE 5

// type used to represent LCP values (B array)
#ifndef BSIZE
#define BSIZE 2
#endif
#if BSIZE==1
typedef uint8_t lcpInt;
#define LCP_FORMAT "%"PRIu8
#define MAX_LCP_SIZE 0xFFULL
#define LCP_EXT "1.lcp"
#define HM_LCP_EXT "1.lcpHM"
#elif BSIZE==2
typedef uint16_t lcpInt;
#define LCP_FORMAT "%"PRIu16
#define MAX_LCP_SIZE 0xFFFFULL
#define LCP_EXT "2.lcp"
#define HM_LCP_EXT "2.lcpHM"
#else
typedef uint32_t lcpInt;
#define LCP_FORMAT "%"PRIu32
#define MAX_LCP_SIZE 0xFFFFFFFFULL
#define LCP_EXT "4.lcp"
#define HM_LCP_EXT "4.lcpHM"
#endif

// type used to represent the Merge (aka Z) array
// defines the maximum number of input BWTs for a single round
// Note that, for LCP computation, we need an open file for each BWT, 
// so the max number of BWTs can be limited in practice by the maximum 
// number of open files allowed by the system  (eg 1024 on my linux box)
#define ZSIZE 1
#if ZSIZE==1
typedef uint8_t palette;
#define PALETTE_FORMAT "%"PRIu8
#define MAX_NUMBER_OF_BWTS 0x100ULL
#elif ZSIZE==2
typedef uint16_t palette;
#define PALETTE_FORMAT "%"PRIu16
#define MAX_NUMBER_OF_BWTS 0x10000ULL
#else
typedef uint32_t palette;
#define PALETTE_FORMAT "%"PRIu32
#define MAX_NUMBER_OF_BWTS 0x100000000ULL
#endif


#define BWT_EXT "bwt"
#define LEN_EXT "size"
#define DOCS_EXT "docs"
#define HM_BWT_EXT "bwtHM"
#define DA_EXT "da"
#define COLORS_EXT "cda" //colored Document array
#define DA_BL_EXT "da_bl"
#define SA_EXT "sa"
#define SA_BL_EXT "sa_bl"
#define QS_EXT "bwt.qs"
#define QS_BL_EXT "bwt.qs_bl"

// structure for writing in a segment of a newMergeColor
typedef struct{
  int fd;  // file descriptor
  off_t offset; // offset inside file descriptor (in bytes)
  palette *buffer; 
  int size; // buffer size  (in palette units)
  int cur;  // current position (in palette units)
} cwriter;


 
typedef struct {
  // input: shared among threads
  symbol **bws;            // bws[0] ... bws[numBwt-1] are input bwt
  customInt *bwtLen;       // bwtLen[i] is size of bwt[i]
  customInt *bwtDocs;      // bwtDocs[i] is the starting DOC-id of each chunk
  customInt **bwtOcc;      // occ of each symbol in each bwt (used ony if smallAlpha==true)
  int numBwt;              // number of bwt to be merged
  customInt mergeLen;      // sum_i bwtLen[i]
  customInt symb_offset;   // offset with respect to the global bwt/lcp (used to read in single lcp file)
  int sizeOfAlpha;         // size of the alphabet   
  bool bwtOnly;            // do not compute LCP
  bool lcpMerge;           // merge LCP values. LCP of input sequences must exist and are overwritten  
  bool lcpCompute;         // compute LCP from scratch
  bool mwXMerge;           // use external multiway mergesort when computing LCP from scratch
  int dbOrder;             // if > 1 output info useful for order-k dbGraph construction (only with -A 128) 
  int outputDA;            // if > 0 output Merge array (=Document Array) for last iteration using outputDA bytes per symbol
  int outputColors;        // if > 0 output Merge array (=Color Array) for last iteration using outputColors bytes per symbol
  int outputSA;            // if > 0 output Merge array (=Suffix Array) for last iteration using outputSA bytes per symbol
  int outputQS;            // if 1 output Merge array (=QS) for last iteration using 1 bytes per symbol
  FILE *unsortedLcp;       // if !NULL file containing unsorted LCP values
  FILE *unsortedLcp_size;  // if !NULL file containing the size of sorted blocks in unsorted_Lcp
  bool smallAlpha;         // the alphabet is small
  bool extMem;             // if true run in external memory
  int algorithm;           // preferred algorithm to use: gap8 gap16 gap128 gap256, if!=8,16,128,256 then use gap  
  bool mmapZ;              // mmap Z arrays
  bool mmapB;              // mmap B array
  bool mmapBWT;            // mmap BWT arrays
  char *lcpinPath;         // base path for lcp files
  char *outPath;           // path for output and  temporary files 
  lcpInt **lcps;           // if lcpMem, lcps[0] ... lcps[numBwt-1] are input lcp values
  int solid_limit;         // irrelevant blocks become solid after this size
  int verbose;
  union {                   // in Gap use one or the other according to bwtOnly 
    lcpInt *blockBeginsAt;  // B array and then merged LCP. output but allocatd by caller  (shared as well)
    uint64_t *bitB;         // B array represented in 2 bits private 
  };
  // external memory access 
  char bwfname[Filename_size];   // filename of the input bwt file 
  char dafname[Filename_size];   // filename of the input DA file 
  char safname[Filename_size];   // filename of the input SA file 
  char qsfname[Filename_size];   // filename of the input QS file 
  FILE **bwf;              // bwf[0] ... bwf[numBwt-1] are pointer inside the input bwt file  
  FILE **daf;              // daf[0] ... daf[numBwt-1] are pointer inside the input DA file  
  FILE **saf;              // saf[0] ... saf[numBwt-1] are pointer inside the input SA file  
  FILE **qsf;              // qsf[0] ... qsf[numBwt-1] are pointer inside the input QS file  
  FILE *fmergeColor;       // mergecolor file   
  cwriter *fnewMergeColor; // newmergecolor files (one per symbol) 
  char *merge_fname;       // name of merge file
  char *newmerge_fname;    // name of new merge file  
  // glocal private
  palette *mergeColor;     // Z
  palette *newMergeColor;  // newZ 
  uint16_t *mergeColor16;  // Z and newZ + 2 bits array combined in a single uint16_t 
  uint32_t *array32;       // Z and newZ + B array combined in a single uint32_t 
  customInt *firstColumn;  // compact first column array
  customInt *F;            // positions inside the newMerge (array F in the pseudocode)   
  customInt *inCnt;        // counters inside each bwt (k_i in the pseudo code)
} g_data; 



// generic data structure for a  1 producer N consumers scheme
typedef struct {
  int pindex;             // index in buffer for producer
  int cindex;             // index in buffer for consumers  
  pthread_mutex_t cmutex; // consumer mutex 
  sem_t *free_slots;       // free slots in the buffer
  sem_t *ready;            // data ready to be processed
  void *buffer;           // generic data buffer
  int buf_size;           // size of buffer
  // from here fields are more or less specific for merge 
  int remaining;          // remaining jobs before a sincronization 
  pthread_mutex_t remutex;// mutex for access to remaining
  pthread_cond_t recond;  // condition variable for remaining 
  union {
    bool hm;                // used for hm vs gap
    g_data *data;
  };
} pc_system;


#endif
