#include "util.h"
#include "mergehm.h"
#if MALLOC_COUNT_FLAG
  #include "malloc_count/malloc_count.h"
#endif

// customInt should represent sizeOfMerge

// input 
//   g->numBwt # bwts to be merged
//   g->bws[]  g->bws[i] is i-th BWT, has size g->bwtLen[i]

// output
//   g->bwtOut & (if !g->BWTOnly) g->lcpOut

static void init_arrays(g_data *g)
{
  customInt i=0; // position inside Z newZ and B 
  for(int j=0;j<g->sizeOfAlpha;j++) {
    if(!g->bwtOnly) g->blockBeginsAt[i]=1;
    g->firstColumn[j] = i;                    // symbol j starts at position i
    for(int b=0;b<g->numBwt;b++)         // scan all occs of symbol j
      for(customInt t=0;t<g->bwtOcc[b][j];t++) {
        if(j==0) { // zero chars are all different, Z are newZ do not change 
          assert(j==0);     // no sqeezing: only zero char is zero
          if(!g->bwtOnly) g->blockBeginsAt[i]=1;
          g->newMergeColor[i]=b;          //for 0-chars write b to both Z and newZ 
        }
        g->mergeColor[i++] = b;
      }
  }
  assert(i==g->mergeLen);
  // extra check on mergeColor, can be commented out
  #ifndef NDEBUG
  customInt cnt[MAX_NUMBER_OF_BWTS] = {0};
  for(i=0;i<g->mergeLen;i++) cnt[g->mergeColor[i]]++;
  bool stop=false;
  for(int i=0; i<g->numBwt; i++)
    if(cnt[i]!=g->bwtLen[i]) {
      printf("INIT %d cnt:"CUSTOM_FORMAT" len:"CUSTOM_FORMAT"\n",i,cnt[i],g->bwtLen[i]);
      stop=true;
    }
  assert(!stop);
  #endif
}

// init Z, newZ and B array computing g->bwtOcc[i][j] and then discarding it
static void init_arrays_largealpha(g_data *g)
{  
  // compute bwtOcc on the spot with a complete scan of input BWTs
  assert(g->bwtOcc==NULL);
  g->bwtOcc = malloc(g->numBwt*sizeof(customInt *));
  if(!g->bwtOcc) die(__func__);
  for(int i=0;i<g->numBwt;i++) {
    g->bwtOcc[i] = calloc(g->sizeOfAlpha,sizeof(customInt));
    if(!g->bwtOcc[i]) die(__func__);
    init_freq_no0(g->bws[i],g->bwtLen[i],g->bwtOcc[i]); 
  }
  init_arrays(g);
  for(int i=0;i<g->numBwt;i++)
    free(g->bwtOcc[i]);
  free(g->bwtOcc);
  g->bwtOcc=NULL;
    
  #ifndef NDEBUG
  // extra check on mergeColor
  customInt cnx[MAX_NUMBER_OF_BWTS] = {0};
  for(int i=0;i<g->mergeLen;i++) cnx[g->mergeColor[i]]++;
  bool stop=false;
  for(int i=0; i<g->numBwt; i++)
    if(cnx[i]!=g->bwtLen[i]) {
      printf("INIT %d cnx:"CUSTOM_FORMAT" len:"CUSTOM_FORMAT"\n",i,cnx[i],g->bwtLen[i]);
      stop=true;
    }
  assert(!stop);
  #endif
}



/* execute a single pass of HM algo 
 * input:
 *   sizeOfAlphabet
 *   initialCharPosition[sizeOfAlphabet] (F array)
 *   sizeOfMerge  (size of the merged BWT) 
 *   mergeColor[sizeOfMerge] (current Z array)
 *   
 * parameters
 *   position (copy of F array)
 *     
 * return
 *   true if the merge vector has changed false otherwise
 * 
 * If only the bwt is required (bwtOnly==true) another iteration is required
 * if the merge vector has changed. If the vector does not change in one
 * iteration, it will not change in subsequent iterations so we are sure 
 * the bwt is the correct one.
 * 
 * It also the lcp is required, we cannot stop if simply the merge vector has
 * not changed, because the bwt is final but the lcp can still change.
 * Consider for example the case s_0 = AATCATC$_0 s_1 = TATCATC$_1. After 4 
 * iterations the bwt is  correctly computed adnthe merge vaector does not change
 * However, we need three more iterations to compute the lcp between 
 * ATCATC$_0 and ATCATC$_1. So we stop the iterations when for two iterations all
 * blocks are monotone. This is done using the boolean var nonMonotoneBlocks
 * and by returning the value 1 when there are only monothone blocks. If we
 * are ot interested in lco values, we return 2 when the merge vector does
 * not change to simplify the holtMcMillan function. 
 * 
 * Note: for the detection of non-monothone blocks we had to get rid
 * of the larget_lcp boolean var that could be uses to avoid some
 * reduntant checking on the B arrays  
 * 
 * Note: we need two distinct arrays merge and mergeNew. 
 * consider the case seq 0 contains only B followed by C and sequence 1 
 * contains only B followed by A. Initially the B region is 0^k 1^h 
 * and at the end it should be 1^h 0^k. However, when the first iteration reaches
 * the B region the 1's have overrwitten the 0's  
 */
static int addCharToPrefix_HM(customInt prefixLength, g_data *g) 
{
  assert(prefixLength<= MAX_LCP_SIZE);
  bool mergeChanged = false; // the Z vector has changed in this iteration
  bool nonMonotoneBlocks = false; // non monotone blocks found in this iteration
  // copy F array
  customInt position[g->sizeOfAlpha];
  array_copy(position, g->firstColumn, g->sizeOfAlpha); //initialize char positions
  // pointer inside each BWT (k_0 & k_1  in the pseudocode)
  array_clear(g->inCnt,g->numBwt,0);
  // blockID: current block for each alphabet char. Used only if blockBeginsAt!=NULL
  customInt blockID[g->sizeOfAlpha];
  array_clear(blockID,g->sizeOfAlpha, g->mergeLen); // mergeLen is an invalid id 
  customInt id = 0;

  // main loop
  for (customInt k = 0; k < g->mergeLen; ++k) {
    if (g->verbose>2) {printf("k="CUSTOM_FORMAT" ",k);}
    int currentColor = g->mergeColor[k]; //  b in pseudocode
    int currentChar = g->bws[currentColor][g->inCnt[currentColor]++]; // c in pseudocode

    // if we are computing the LCP check if we are entering next block
    if (!g->bwtOnly && g->blockBeginsAt[k] && (g->blockBeginsAt[k] < prefixLength)) 
      id = k; 

    // write color in new Z array, except 0 chars
    if(currentChar!=0) {
      customInt positionToUpdate = position[currentChar]++;
      g->newMergeColor[positionToUpdate] = currentColor;
      if(g->mergeColor[positionToUpdate]!=currentColor)
        mergeChanged=true;
      if (!g->bwtOnly) {
        // create new block?
        if(blockID[currentChar] != id) {
            if(g->blockBeginsAt[positionToUpdate]==0) {// only 0 values in B should be overwritten
              g->blockBeginsAt[positionToUpdate] = prefixLength;
            }
          blockID[currentChar] = id;
        }
        // the previous block continues, check if it is non-monotone
        else { 
          assert(positionToUpdate>0);
          nonMonotoneBlocks = nonMonotoneBlocks||(g->newMergeColor[positionToUpdate-1]!=currentColor);
        }
      }
    }
    else // extra check for has0 chars, can be commented out 
      assert(currentChar==0);
  } // end main loop 
  
  // check we have seen all BWT positions 
  for(int i=0; i<g->numBwt; i++)
    assert(g->inCnt[i]==g->bwtLen[i]);

  // extra check on newMergeColor
  #ifndef NDEBUG
  bool stop=false;
  customInt cnt[MAX_NUMBER_OF_BWTS] = {0};
  for(customInt i=0;i<g->mergeLen;i++) 
    cnt[g->newMergeColor[i]]++;
  for(int i=0; i<g->numBwt; i++)
    if(cnt[i]!=g->bwtLen[i]) {
      printf("NEWMERGE %d cnt:"CUSTOM_FORMAT" len:"CUSTOM_FORMAT"\n",i,cnt[i],g->bwtLen[i]);
      stop=true;
    }
  assert(!stop);
  #endif

  // swap Merge and newMerge
  palette *tmp=g->mergeColor; g->mergeColor=g->newMergeColor; g->newMergeColor=tmp;
  if(g->bwtOnly && !mergeChanged)
    return 2; // no lcp and bwt not changed we can stop immediately
  if(!g->bwtOnly  && nonMonotoneBlocks==false)
    return 1; // al blocks are monothone, stop at next iteration 
  return  0; 
}

/**
 * Main entry point for the Holt-McMillan algorithm
 * 
 * Note: allocate
 *   mergeColor[mergeLen]  
 *   newMergeColor[mergeLen]  
 *   blockBeginsAt[mergeLen] 
 *   firstColumn[alphaSize]
 *   inCnt[numBwt]
 * */
void holtMcMillan(g_data *g, bool lastRound) {
  // init local global vars
  check_g_data(g);
  assert(g->numBwt <= MAX_NUMBER_OF_BWTS); 
  int stop = 0; 
  
  // clear array B if necessary
  if(g->bwtOnly) assert(g->blockBeginsAt==NULL);
  else alloc0_B_array(g); // alloc and clear blockBeginsAt array

  // allocate Z (merge) Znew 
  alloc_merge_arrays(g);
  // allocate firstColumn and inCnt
  g->inCnt = malloc(g->numBwt*sizeof(customInt)); 
  g->firstColumn = malloc(g->sizeOfAlpha*sizeof(customInt));
  if(!g->firstColumn || !g->inCnt) die(__func__);
  // init the above arrays
  if(g->smallAlpha) init_arrays(g);
  else init_arrays_largealpha(g);

  customInt prefixLength = 1;  
  do { // main loop. Add something to do nothing if there is a single BWT as in Gap? 
    prefixLength+= 1;
    if(prefixLength>MAX_LCP_SIZE && !g->bwtOnly) {fprintf(stderr,"LCP too large\n");die(__func__);}
    stop += addCharToPrefix_HM(prefixLength, g);
    if(g->verbose>2 || (g->verbose>1 && lastRound)) { 
      printf("Lcp: "CUSTOM_FORMAT" Stop=%d   ",prefixLength-1,stop);
      #if MALLOC_COUNT_FLAG
        printf("Mem: %zu peak, %zu current, %.2lf/%.2lf bytes/symb\n", malloc_count_peak(),
             malloc_count_current(), (double)malloc_count_peak()/g->mergeLen,
             (double)malloc_count_current()/g->mergeLen);
      #endif
    }
  } while(stop<2);
  
  // computation complete, do the merging. The following call writes the
  // (possibly remapped) merged BWT back to g->bws[0]; and if lcpMerge==true the merged LCP to g->lcps[0]  
  mergeBWTandLCP(g,lastRound);

  free(g->firstColumn); // last four arrays deallocated
  free(g->inCnt);
  free_merge_arrays(g);
  if(!g->bwtOnly) free_B_array(g);
}
