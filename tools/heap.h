#ifndef HEAP_H
#define HEAP_H

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <limits.h>
#include <string.h>
#include <assert.h>


typedef unsigned __int128 uint128_t;

static const uint64_t mask[9] = { 
0,
0xFF,             // 1 byte
0xFFFF,           // 2 bytes
0xFFFFFF,         // 3 bytes
0xFFFFFFFF,       // 4 ... 
0xFFFFFFFFFF,     // 5 
0xFFFFFFFFFFFF,   // 6  
0xFFFFFFFFFFFFFF, // 7
0xFFFFFFFFFFFFFFFF// 8
};

//#define pos(i) ((i&mask[h->pos_size]))
//#define lcp(i) ((i>>((h->pos_size)*8)))

#define pos(i) ((uint64_t)(i>>((h->lcp_size)*8)))
#define lcp(i) ((uint64_t)(i&mask[h->lcp_size]))

//sorting key is <pos>
#define key(i) ((uint128_t)((h->heap[i]->buffer[h->heap[i]->idx])))
#define MAX_KEY ((uint128_t)(~0ULL)&mask[h->lcp_size+h->pos_size])

#define heap_size(h) ((h)->size)

#ifndef DEBUG
#define DEBUG 1
#endif

#ifndef INPUT_SIZE
#define INPUT_SIZE 1024/sizeof(pair) //1K
#endif

#ifndef OUTPUT_BUFFER
#define OUTPUT_BUFFER 1 
#endif

#if OUTPUT_BUFFER
#define OUTPUT_SIZE 4096/sizeof(pair) //4K
#endif

//typedef struct pair{
//int pos, lcp;
//} pair;

//typedef uint64_t pair;
typedef uint128_t pair;

// struct representing a file of [position|lcpValue] pairs
typedef struct heap_node{
  //uint64_t key;         // next key to be sorted
  //FILE *f_in;           // file containing the other pairs
  //size_t  seek;
  //uint64_t tot_size;    // total number of lcp values in file
  						// useful only for an extra check
  pair *buffer;	
  int idx; 
  size_t seek;
  
} heap_node;


typedef struct heap {
  heap_node** heap;
  
  int size;  
  char file_name[500];
  
  FILE* f_in;
  
  #if OUTPUT_BUFFER
    pair *out_buffer;
    int out_idx;
  #endif
  
  int pos_size;
  int lcp_size;

} heap;

/**********************************************************************/

heap* heap_alloc(int n, char* file_name, int level, int pos_size, int lcp_size);
void heap_free(heap* h, FILE *f_out, int level);

int heap_insert(heap *h, size_t pos);
pair heap_delete_min(heap *h);

void heap_write(heap *h, FILE* f_out, pair tmp, int level);


#endif
