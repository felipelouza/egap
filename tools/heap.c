#include <errno.h>
#include <math.h>
#include <stdlib.h>

#include "heap.h"

/********************************************************************************
  Allocates an empty heap for at most n nodes.

  Returns a pointer to the heap on success.  Returns null if no space is left in
  memory and sets errno to ENOMEM.
********************************************************************************/
heap* heap_alloc(int n, char* file_name, int level, int pos_size, int lcp_size) {
   
  heap* h = (heap*) malloc(sizeof(heap));

  if (!h) {
    errno = ENOMEM;
    return 0;
  }
  
  h->heap = (heap_node**) malloc((n+2)*sizeof(heap_node*));
  if (!h->heap) {
    free(h);
    errno = ENOMEM;
    return 0;
  }
  
  int i;
  for(i=0; i<n; i++)
    h->heap[i] = NULL;
   
  h->size = 0;
  
  //input file for all heap nodes
  strcpy(h->file_name, file_name);
 
  h->f_in = fopen(h->file_name, "rb"); 
  if (!h->f_in) perror ("file_open"); 

  #if OUTPUT_BUFFER
    h->out_buffer = (pair*) malloc((OUTPUT_SIZE+1)*sizeof(pair));
    if(!h->out_buffer) perror("malloc(heap_alloc)");  
    h->out_idx = 0;
  #endif

  h->pos_size = pos_size;
  h->lcp_size = lcp_size;

  return h;
}

#if OUTPUT_BUFFER

void write_buffer(heap *h, FILE *f_out, int level) {

  int i;
  if(level){
    for(i=0; i<h->out_idx; i++){
      fwrite(&h->out_buffer[i], h->pos_size+h->lcp_size, 1, f_out);
    }
  }
  else{
    for(i=0; i<h->out_idx; i++){
      uint64_t lcp = lcp(h->out_buffer[i]);
      fwrite(&lcp, h->lcp_size, 1, f_out);
    }
  }
  h->out_idx = 0;
}

#endif

/********************************************************************************
  Releases a heap from memory.
********************************************************************************/
void heap_free(heap *h, FILE *f_out, int level) {

  int i;
  for (i=0; i<h->size; i++){
    //fclose(h->heap[i]->f_in);
    free(h->heap[i]->buffer);
    free(h->heap[i]);
  }

  #if OUTPUT_BUFFER 
    if(h->out_idx) write_buffer(h, f_out, level);
    free(h->out_buffer);
  #endif
  
  fclose(h->f_in);
      
  free(h->heap);
  free(h);
}

/********************************************************************************
   Moves the node at position i up in the heap.
********************************************************************************/
void heapfy_up(heap *h, int i) {

  int p = (int) floor((i-1)/2);

  while (i>0 && key(i) < key(p)) {  
      
    heap_node *auxn = h->heap[i];
    h->heap[i] = h->heap[p];
    h->heap[p] = auxn;

    i = p;
    p = (int) floor((i-1)/2);
  }
}

/********************************************************************************
  Moves the node at position i down in the heap.
********************************************************************************/
void heapfy_down(heap *h, int i) {

  int last = h->size-1;

  // Selects the child s of i with smaller cost:
  int l = 2*i+1;
  int r = l+1;
  int s;
  if (r <= last)
    s = key(l) < key(r) ? l : r;    
  else if (l == last)
    s = l;
  else
    return;

  // Swaps down:
  while (key(i) > key(s)) {   
    
    heap_node *aux = h->heap[i];
    h->heap[i] = h->heap[s];
    h->heap[s] = aux;

    i = s;
    l = 2*i+1;
    r = l+1;

    if (r <= last)
      s = key(l) < key(r) ? l : r;      
    else if (l == last)
      s = l;
    else
      return;
  }
}
/********************************************************************************
  Reads the file into the buffer.

********************************************************************************/

void  read_buffer(heap *h, heap_node *node){

//  fseek(h->f_in, node->seek*sizeof(pair), SEEK_SET);
    fseek(h->f_in, node->seek*(h->pos_size+h->lcp_size), SEEK_SET);

//  fread(node->buffer, h->pos_size+h->lcp_size, INPUT_SIZE, h->f_in); 
    int i;
    for(i=0; i<INPUT_SIZE; i++){
      node->buffer[i]=0;
      int e = fread(&node->buffer[i], h->pos_size+h->lcp_size, 1, h->f_in);
      if(e!=1) {perror(__func__); exit(1);}
      if(node->buffer[i]==MAX_KEY) break; 
    }

    node->seek += INPUT_SIZE;
}

/********************************************************************************
  Inserts a node in the heap. Open the file and load the buffer.

  Returns 0 on success.  If no space is left in memory then returns ENOMEM and
  sets errno to ENOMEM.

********************************************************************************/
int heap_insert(heap *h, size_t pos) {

  if (!h->heap[h->size])
    h->heap[h->size] = malloc(sizeof(heap_node));

  if (!h->heap[h->size]) 
    return errno = ENOMEM;

  heap_node *node = h->heap[h->size];

  //alloc buffer <pos, lcp>
  node->buffer = (pair*) malloc(INPUT_SIZE*sizeof(pair));
  if (!node->buffer) 
    return errno = ENOMEM;
  node->idx = 0;
  
  //load buffer
  node->seek = pos;
  read_buffer(h, node);

  #if DEBUG
    int i;
    for(i=0; i<INPUT_SIZE; i++) 
      printf("<%lu, %lu (%lu)> ", lcp(node->buffer[i]), pos(node->buffer[i]), node->buffer[i]);
    printf("**\n");
  #endif
  
  heapfy_up(h,h->size++);

  return 0;
}

/********************************************************************************
  Removes the key with minimum cost <pos, lcp> and returns it.

  The caller must ensure that the heap is not empty.
********************************************************************************/
pair heap_delete_min(heap *h) {

  heap_node *node = h->heap[0];

  pair tmp=0;
  tmp= node->buffer[node->idx];

  if(++node->idx == INPUT_SIZE){ //load from disk
    node->idx = 0; 
    read_buffer(h, node);
  }

  heapfy_down(h,0);
  
  return tmp;
}

/********************************************************************************
  Writes the pair <pos, lcp> if level == 1, otherwise it writes <lcp>.

********************************************************************************/
void heap_write(heap *h, FILE* f_out, pair value, int level){
  
  if(level){
    #if OUTPUT_BUFFER
      h->out_buffer[h->out_idx++]= value;
      if(h->out_idx==OUTPUT_SIZE){
        write_buffer(h, f_out, level);
        //fwrite(h->out_pair_buffer, sizeof(pair), h->out_idx, f_out);  
        //h->out_idx=0;
      }   
    #else
      //uint64_t pos = pos(value);
      //uint64_t lcp = lcp(value);
      //fwrite(&lcp, h->lcp_size, 1, f_out);
      //fwrite(&pos, h->pos_size, 1, f_out);
      fwrite(&value, h->pos_size+h->lcp_size, 1, f_out);
    #endif
  }
  else{
    #if OUTPUT_BUFFER   
      //h->out_lcp_buffer[h->out_idx] = value.lcp;
      h->out_buffer[h->out_idx++] = value;
      if(h->out_idx==OUTPUT_SIZE){
        write_buffer(h, f_out, level);
      }
    #else
      uint64_t lcp = lcp(value);
      fwrite(&lcp, h->lcp_size, 1, f_out);  
    #endif
  }
  
}
