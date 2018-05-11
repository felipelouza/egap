#include "alphabet.h"
#include "util.h"

// ===== functions to reduce/enlarge alphabet

// map from input alphabet -> restricted alphabet
static int alphabet_map[SIZE_OF_ALPHABET];
// map from input restritcted alphabet -> alphabet
static int restricted_unmap[SIZE_OF_ALPHABET];


// compute frequencies of symbols in b[n]
// and store them in freq (which should be of len SIZE_OF_ALPHABET)
void init_freq(symbol *b, customInt n, customInt *freq)
{
  for(int i=0;i<SIZE_OF_ALPHABET;i++)
    freq[i]=0;
  for(customInt i=0;i<n;i++)
    freq[b[i]]++;
}
#if 0
void init_freq_ext(char *fname, customInt n, customInt *freq)
{
  for(int i=0;i<SIZE_OF_ALPHABET;i++)
    freq[i]=0;
  FILE *f = fopen(fname,"r");  
  symbol b;
  for(customInt i=0;i<n;i++) {
    int e = fread(&b,sizeof(b),1,f);
    if(e!=1) die(__func__);
    freq[b]++;
  }
  fclose(f);
}
#endif

void init_freq_no0(symbol *b, customInt n, customInt *freq)
{
  for(customInt i=0;i<n;i++)
    freq[b[i]]++;
}


/**
 * Compute alphabet map and its inverse for a string b[n]
 * \return   number of distinct symbols in b
 * */ 
int init_alpha_maps(g_data *g)
{
  
  // init maps with illegal value -1
  for(int i=0;i<SIZE_OF_ALPHABET;i++)
    alphabet_map[i]=restricted_unmap[i] = ILLEGAL_SYMBOL;

  // compute frequencies of the input symbols
  customInt freq[SIZE_OF_ALPHABET];
  customInt n = g->mergeLen;
  init_freq(g->bws[0],n,freq); // complete scan of bws[0][0..n]
  // compute maps
  int j=0;
  for(int i=0;i<SIZE_OF_ALPHABET;i++)
    if(freq[i]!=0) {
      alphabet_map[i] = j;
      restricted_unmap[j++] = i;
    }
  return j;
}

int alpha_reduce(int c)
{
  assert(alphabet_map[c]!= ILLEGAL_SYMBOL);
  return alphabet_map[c];
}

int alpha_enlarge(int c)
{
  assert(restricted_unmap[c]!= ILLEGAL_SYMBOL);
  return restricted_unmap[c];
}

// remap the string b according to alpha_reduce
// if freq!=null update freq[] according to the remapped values
// freq[] must be appropriately allocated and initialized
static void remap_string(symbol *b, customInt n, customInt *freq)
{
  for(customInt i=0;i<n;i++) {
    b[i] = alpha_reduce(b[i]);
    if(freq) freq[b[i]]++;
  }  
}

// remap BWTs and init g>bwtOcc[i] and g->sizeOfAlpha
void remap_bwts(g_data *g)
{
  int sizeOfAlphabet; // size of alphabet
  
  // remap alphabet, remap bwts and compute occ 
  sizeOfAlphabet = init_alpha_maps(g);
  assert(sizeOfAlphabet>1 && sizeOfAlphabet<= SIZE_OF_ALPHABET);
  if(g->verbose>0) printf("Alphabet size: %d\n", sizeOfAlphabet);
  
  // remap; also compute Occs for each BWT if smallAlpha
  if(g->smallAlpha) {
    g->bwtOcc = malloc(g->numBwt*sizeof(customInt *));
    if(!g->bwtOcc) die(__func__);
    g->bwtOcc[0] = calloc(g->numBwt*sizeOfAlphabet,sizeof(customInt));
   if(!g->bwtOcc[0]) die(__func__);
  }  
  for(int i=0;i<g->numBwt;i++) {
    customInt *bwtOcc = NULL;
    if(g->smallAlpha) { // if the alphabet is small keep bwtOcc[i]
      if(i>0) g->bwtOcc[i] = g->bwtOcc[i-1] + sizeOfAlphabet;
      bwtOcc = g->bwtOcc[i];
    }
    remap_string(g->bws[i],g->bwtLen[i],bwtOcc);
  }
  g->sizeOfAlpha = sizeOfAlphabet;
}

