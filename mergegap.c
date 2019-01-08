#include "mergegap.h"
#include "util.h"
#include "io.h"
#include "alphabet.h"
#include "malloc_count/malloc_count.h"

#define min(a,b) ((a)<(b) ? (a) : (b))
#define max(a,b) ((a)>(b) ? (a) : (b))


// the gap routine supports:
//   lcpMerge, lcpCompute
// external memory storage for the input BWTs and the merge arrays 
// if !lcpMerge the RAM usage is only for the 2n bit arrays  


// --- summary of the available construction algorithms ---

// assuming a symbol takes 1 byte, we merge at most 256 Bwt, each LCP value takes 2 bytes
// these are the space usage of the different versions.   

// For only merging the BWT the space in bytes usage of the different versions:
//   gap:       n (BWTs) 2n (Z) n/4 (B) = 3.25 n [extMem: n/4]
//   gap128:    n (BWTs) 2n (Z+B)       = 3 n    [extMem 3n could be reduced to 2n]
//   gap16:     n (BWTs) n (Z) n/4 (B)  = 2.25 n [extMem idem, could be reduced to 1.25 n]
//   gap8:      n (BWTs) n (Z+B)        = 2 n    [extMem n]

// for merging the BWT and the LCP (external memory not supported) 
//   gap:       n (BWTs) 2n (Z) 2n (BlockB) = 5 n
//   gap16:     n (BWTs) n (Z)  2n (BlockB) = 4 n
//   

// for merging the BWT and computing the LCP with the compute from scratch procedure 
//   gap:       3.25 n [extMem: n/4]
//   gap8:      2 n    [extMem n]
//   gap128:    3 n    [extMem 3n could be reduced to 2n]
//   gap256:    n (BWTs) 4n (Z+BlockB)  = 5 n    [extMem 5n, could be reduced to 4n]
//              Note: gap256 is used only to avoid the lcp merge step (option -x)
//                    but there is no good reason to do that.  

// gap128ext computes bwt and possibly LCP completely in external memory 


// meaning of the g->algorithm parameter:
// N = 8,16,256 --> use gapN algorithm if possible
// N=128 --> use gap128ext if g->extMem of gap128 otherwise (again if possible)
// N=0   --> use the old best fit strategy minimizing the amount of RAM 
// any other value --> use gap (to force the use of gap use -A 256 without -x) 


// input from variables stored in g 
//   mergeLen
//   sizeOfAlpha
//   numOfBwt  
//   others...

// output
//   the final merge overwritten to bws[0]
//   if lcpMerge==true LCP values are stored to lcps[]
//   if lcpCompute==true and lastRound, lcp values are stored to the .pair files

// used to access BWTs in external memory 
static void open_merge_files(g_data *g);
static void close_merge_files(g_data *g);
static void fwrite_color(int b, FILE *f);
static  int fread_color(FILE *f);




#include "blocks.h"
#include "merge8.h"        // ext:BWTs lcpCompute !lcpMerge
#include "merge16.h"       // !lcpCompute lcpMerge !ext
#include "merge128.h"      // lcpCompute !lcpMerge !ext
#include "merge128ext.h"   // lcpCompute !lcpMerge ext:BWTs:Z:B
#include "merge256.h"      // lcpCompute (without mergesort) !lcpMerge !ext: do not use for bwtOnly


/**
 * Using the number of occs of each symbol in each bwt (stored in bwtOcc)
 * init the array Z (mergeColor) and B (blockBeginsAt) at the value
 * they should have after the first iteration of the Gap algorithm:
 *   blockBeginsAt[i]=1 if i is the first occurrence in the first
 *                      column F of a new symbol or F[i]=0
 *                      (0 occurrences are assumed to be all different)
 *   in each region of the F column with the same symbol j in Z
 *   we have: #occ(j) in bwt[0], #occ(j) in bwt(1), and so on
 * Since the region corresponding to 0 does not change and 0 has a 
 * special update rule, we init that region also in Znew (newMergeColor)
 * and we never modify that region in the algorithm.   
 * The array firstColumn (compact representation of F) is also initialized
 *  
  * Note, when computing only the BWT instead of B we init bitB which
 * uses 2 bits per entry to encode the values: 
 *   never set->00,  recently set->01 or 10,  set at least 2 iterations before->11
 * during this initialization we write 01 for recently set entries,
 * therefore in the first iteration the mask for access to bitB should be 10 (eg 2)  
* */    
 
// init Z, newZ B, and first Column array using g->bwtOcc[i][j]
static void init_arrays(g_data *g)
{
  customInt i=0; // position inside Z newZ and B  
  FILE *fnewmerge=NULL, *fmerge=NULL;
  if(g->extMem) {
    fnewmerge = fopen(g->newmerge_fname,"wb");
    if(!fnewmerge) die("mergegap:init_arrays:fnewmerge open");
    fmerge = fopen(g->merge_fname,"wb");
    if(!fmerge) die("mergegap:init_arrays:fmerge open");
  }
  for(int j=0;j<g->sizeOfAlpha;j++) {
    if(!g->lcpMerge) tba_or_m(g->bitB,i,1);
    else g->blockBeginsAt[i]=1;  // start of symbol j, correct lcp is 0
    g->firstColumn[j] = i;  // symbol j starts at position i
    for(int b=0;b<g->numBwt;b++) {
      for(customInt t=0;t<g->bwtOcc[b][j];t++) { 
        if(j==0) { // zero chars are all different, Z, Znew never change 
          if(!g->lcpMerge) tba_or_m(g->bitB,i,1);
          else g->blockBeginsAt[i]=1;
           //for 0-chars write b also to newmerge aka newZ 
          if(g->extMem)  fwrite_color(b,fnewmerge);
          else           g->newMergeColor[i]=b;          
        }
        // all colors are written to merge aka Z
        if(g->extMem) fwrite_color(b,fmerge);
        else g->mergeColor[i] = b;
        i++;
      } // end for t
    } // end for b 
  } // end for j
  assert(i==g->mergeLen); 
  // extra check on mergeColor, can be commented out
  if(g->extMem) {
    assert(ftell(fmerge)==g->mergeLen*sizeof(palette));
    if(fclose(fmerge)!=0) die("init_arrays:fmerge close"); 
    if(fclose(fnewmerge)!=0) die("mergegap:init_arrays:fnewmerge close"); 
  }
  else {
    #ifndef NDEBUG
    customInt cnt[MAX_NUMBER_OF_BWTS] = {0};
    for(i=0;i<g->mergeLen;i++) cnt[g->mergeColor[i]]++;
    for(int i=0; i<g->numBwt; i++)
      assert(cnt[i]==g->bwtLen[i]);
    #endif
  }
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
}


// single iteration of the Gap algorithm
// input is head of the irrelevant lists (fin and fout) and an empty liquid block
// return true if all sequencen has become irrelevant.  
// when working in external memory assumes that g->merge_fname and g->newmerge_fname 
// hold the name of the files containing the merge and newmerge array
// during the iteration colors are read from merge sequentially and written to newmerge 
// in positions corresponding to the nonzero characters. At the end of the iteration, 
// the file names for merge and newmerge are swapped so that the one with name merge points to
// the current merge array 
static bool addCharToPrefix(solidBlockFile *solidHead, liquidBlock *liquid, customInt prefixLength, bool *mergeChanged, const int round, g_data *g) {
  assert(liquid->empty);
  liquid->beginsAt = liquid->endsAt = 0;  
  for(int i=0;i<liquid->occ_size;i++) assert(liquid->occ[i]==0);
  // copy first column to F 
  array_copy(g->F, g->firstColumn, g->sizeOfAlpha); //initialize char positions
  // pointer inside each BWT (k_0 & k_1 in the pseudocode)
  array_clear(g->inCnt,g->numBwt,0);
  if(g->extMem) {
    rewind_bw_files(g);  // set file pointers at the beginning of each BWT 
    open_merge_files(g); // open merge file for reading and newmerge files for writing
  }
  // id for each character, init with an invalid id
  customInt blockID[g->sizeOfAlpha];
  array_clear(blockID,g->sizeOfAlpha, g->mergeLen);  // mergeLen is an invalid id
  customInt id = 0, k; 
  int m = (round%2==1) ? 1 : 2;  // mask for the bitB array
  uint64_t lcpWritten =0;

  protoBlock cblock = {.mono = false};
  solidBlock *next = readBlock(solidHead); // first block
  solidBlock *last = NULL;                 // previous block

  for (k = 0; k < g->mergeLen; ) { 
    assert(next==NULL || k <= next->beginsAt);  // we did not pass next block
    assert(last==NULL || last->nextBlock == next); // last is the immediately preceeding block

    // check if we are entering a block, and if the block is at least 2 iterations old
    bool start_block, last_block_recent=true; // for k=0 a new block starts, so last is properly initialized
    if(g->lcpMerge) {
      start_block = (g->blockBeginsAt[k]>0) && (g->blockBeginsAt[k] < prefixLength);
      if(start_block) last_block_recent = g->blockBeginsAt[k]>=prefixLength-1;
    }
    else 
      start_block = tba_block_test_set(g->bitB,k,m,&last_block_recent);    
    if(start_block && last_block_recent && g->lcpCompute)
      {writeLcp(k,prefixLength-2,g); lcpWritten++;} // save lcp value found in previous iteration
    if (start_block) {
      // if the block we just left is a singleton we add it to liquid that remains active
      if(!last_block_recent && cblock.mono==true && (cblock.beginsAt==k-1)) {
        cblock.endsAt = k;           // solidifiable singleton block just ended
        add_singleton2liquid(&cblock, liquid);
      }
      // if the block we left is not recent, monochrome, we are not computing LCPs and not using extMem add it  
      else if(!last_block_recent && cblock.mono==true && !g->lcpCompute && !g->extMem) {
        cblock.endsAt = k;           // solidifiable monochorome block just ended
        add_proto2liquid(&cblock,liquid);  // add proto to liquid that remains active 
      }
      else { // proto block cannot be added, close current liquid
        if(!liquid->empty)
          last = finalize_liquid(last,liquid,next,solidHead); // this is the only point where a new block is created
        assert(liquid->empty);
        liquid->beginsAt=liquid->endsAt=k; // start empty liquid block
      }
      assert(liquid->endsAt==k);
      // block ending at k considered, now look forward 
      if(next!=NULL && next->beginsAt==k) { // entering an irrelevant block
        skip(next, g);                      // skip block
        k = next->endsAt;                   // update k
        // merge liquid with next block and possibly previous 
        if(last==NULL || last->endsAt!=liquid->beginsAt) {
            if(!liquid->empty) merge_liquid(liquid,next,solidHead); // simple merge
            if(last!=NULL) writeBlock(last,solidHead);  // save current last  
            last = next;   // advance last           
        }
        else  //three way merge, next is freed, last does not change
          merge_sls(last,liquid,next,solidHead); // only point where a solid block can be destroyed  
        assert(liquid->empty);
        liquid->beginsAt=liquid->endsAt=k; // start empty liquid block
        next = readBlock(solidHead);       // next has become last, update next (was: next = last->nextBlock; )
        last->nextBlock = next;
        assert(k==last->endsAt);
        cblock.mono = false;        // prevent re-adding the just skipped block
        continue;                   // resume from the end of the block 
      }
      // we are entering a relevant block, unless it is a recent one make it a candidate for solidification  
      if( !last_block_recent ) {
        cblock.beginsAt = k; cblock.mono=true; 
        if(!g->extMem) { // monochrome blocks of length>1 are not used in external memory  
          cblock.color = g->mergeColor[k];
          cblock.start = &g->bws[cblock.color][g->inCnt[cblock.color]]; // bwt-position of first char in block
        }
      }
      else cblock.mono = false; // not a candidate for solid block, wait next iteration 
      if(last_block_recent) 
         id = k;    // id of the new block 
    }
    // processing a char in a relevant block
    int currentColor=0;   // b in pseudocode
    int currentChar=0;  // c in pseudocode: next char in current BWT 
    if(g->extMem) {
      assert(ftell(g->fmergeColor)==k*sizeof(palette));
      currentColor = fread_color(g->fmergeColor);
      assert(ftell(g->bwf[currentColor])==(g->inCnt[currentColor]+(g->bws[currentColor]-g->bws[0])+g->symb_offset)*sizeof(symbol)); 
      int e = fread(&currentChar,sizeof(symbol),1,g->bwf[currentColor]);
      if(e!=1) die("mergegap:addCharToPrefix:bwt[color]Read"); 
      g->inCnt[currentColor]++;
    }
    else {
      currentColor = g->mergeColor[k];
      currentChar  = g->bws[currentColor][g->inCnt[currentColor]++]; // c in pseudocode
    }
    // add currentChar/Color to proto block  
    cblock.lastChar =  currentChar;   // save lastchar, only useful for singleton blocks
    cblock.lastColor =  currentColor; // save lastcolor, only useful for singleton blocks
    if(!g->extMem)
      if(currentColor != cblock.color) cblock.mono = false;       // block is not monochrome
    // write color in new Z array, except 0 chars
    if(currentChar!=0) {
      customInt positionToUpdate = g->F[currentChar]++;
      if(g->extMem) {
        cwriter_put(&g->fnewMergeColor[currentChar],currentColor);
        assert(cwriter_tell(&g->fnewMergeColor[currentChar])==g->F[currentChar]*sizeof(palette));
        *mergeChanged=true; // this could be a problem....
      }
      else {
        g->newMergeColor[positionToUpdate] = currentColor;
        if(g->bwtOnly && !*mergeChanged && g->mergeColor[positionToUpdate]!=currentColor)
          *mergeChanged=true;  // remember there is a difference from the previous iteration
      }
      // create new block?
      if (blockID[currentChar] != id) {
        if(!g->lcpMerge) { // no lcp just mark B array 
          if(last_block_recent) tba_mark_if0(g->bitB,positionToUpdate,m);
        }
        else  // update lcp
          if(last_block_recent && g->blockBeginsAt[positionToUpdate]==0) // only 0 values in B are overwritten
            g->blockBeginsAt[positionToUpdate] = prefixLength;
        blockID[currentChar] = id; // update block id, always!
      }
    }
    k++;
  } // end main loop
  assert(k==g->mergeLen);
  assert(next==NULL); 
  if(cblock.mono==true && cblock.beginsAt==k-1) {
    cblock.endsAt = k;           // solidifiable singleton block just ended
    add_singleton2liquid(&cblock, liquid);
    assert(!liquid->empty);
  }
  else if(cblock.mono==true && !g->lcpCompute && !g->extMem) {
    cblock.endsAt = k;
    add_proto2liquid(&cblock,liquid);  // add proto to liquid that remains active
    assert(!liquid->empty);
  }
  if(!liquid->empty) 
    last = finalize_liquid(last,liquid,NULL,solidHead); // a new block could be created 
  assert(liquid->empty);
  // add EOF value to lcp file and entry to .size file 
  if(g->lcpCompute && lcpWritten>0)
    writeLcp_EOF(++lcpWritten, g);
  // check if all sequence has become irrelevant
  bool everything_irrelevant = false;
  if(last!=NULL && last->beginsAt==0 && last->endsAt==g->mergeLen)
    everything_irrelevant = true;
  // save last block   
  if(last!=NULL) writeBlock(last,solidHead);
  // check we have read all chars from all BWT's 
  for(int i=0; i<g->numBwt; i++)
    assert(g->inCnt[i]==g->bwtLen[i]);
    
  // swap merge and newMerge
  if(g->extMem) { // close merge files and swap file names 
    assert(ftell(g->fmergeColor)==g->mergeLen*sizeof(palette));
    assert(g->F[g->sizeOfAlpha-1]==g->mergeLen);
    assert(cwriter_tell(&g->fnewMergeColor[g->sizeOfAlpha-1])==g->mergeLen*sizeof(palette));
    close_merge_files(g);
    char *tmp=g->merge_fname; g->merge_fname = g->newmerge_fname; g->newmerge_fname = tmp;
  } else { // swap arrays 
    palette *tmp=g->mergeColor; g->mergeColor=g->newMergeColor; g->newMergeColor=tmp;
  }
  return everything_irrelevant;
}


// entry point for the gap bwt/lcp merging procedures (including gap8 gap16 etc)
// if we are only interested in BWT merge, blockBeginsAt is replaced by a bit array 
void gap(g_data *g, bool lastRound) {
  
  // maybe there is nothing to do
  if(g->numBwt<2 && !g->lcpCompute && g->dbOrder==0) {
    assert(g->numBwt==1);
    if(g->verbose>0) puts("Single BWT/LCP merging: nothing to do!");
    return;
  }
  // try preferred algorithm
  if(g->algorithm==8 && g->numBwt <=8 && !g->lcpMerge)
    return gap8(g,lastRound);
  else if(g->algorithm==16 && g->numBwt <=16 && !g->extMem &&!g->lcpCompute)
    return gap16(g,lastRound);    
  else if(g->algorithm==128 && g->numBwt <=128 && !g->lcpMerge)
    return g->extMem ? gap128ext(g,lastRound) : gap128(g,lastRound);    
  else if(g->algorithm==256 && g->numBwt <=256 && g->lcpCompute && !g->mwXMerge)
    return gap256(g,lastRound);
  else if(!g->algorithm) {      
    // use a best fit strategy
    if(g->numBwt<=128 && g->extMem)
      return gap128ext(g,lastRound); // gap128ext is the best extMem also since uses o(n) RAM 
    // case of lcpCompute or bwtOnly not in external memory for at most 128 bwts
    if(!g->lcpMerge && g->numBwt <= 128 && !g->extMem) 
      return gap128(g,lastRound);  // possibly use gap8 or gap16      
    // case <=16 with lcpMerge not in external memory
    if(g->numBwt<=16 && g->lcpMerge && !g->extMem) return gap16(g,lastRound);
  }

  // everything else is handled here!
  if(g->verbose>0) {
    if(g->lcpCompute) puts("BWT Merge + LCP compute with gap");
    else if(g->lcpMerge) puts("BWT/LCP Merge with gap");  
    else puts("BWT only merging with gap");
  }
  assert(g->numBwt <= MAX_NUMBER_OF_BWTS);
  if(g->lcpCompute) {       // we compute LCP values only if we are at the last round 
    assert(!g->bwtOnly && lastRound);
    open_unsortedLCP_files(g);
  }
  else assert(g->bwtOnly || g->lcpMerge);
  
  // init local global vars
  check_g_data(g);
  if(g->extMem) open_bw_files(g);
  // allocate and clear bit/int array B 
  if(!g->lcpMerge) g->bitB = tba_alloc(g->mergeLen, g->mmapB);
  else alloc0_B_array(g); // alloc and clear blockBeginsAt array
        
  // allocate Z (merge) Znew 
  alloc_merge_arrays(g);
  // allocate other useful arrays
  g->inCnt = malloc(g->numBwt*sizeof(customInt));
  g->firstColumn = malloc(g->sizeOfAlpha*sizeof(customInt));
  g->F = malloc(g->sizeOfAlpha*sizeof(customInt));
  if(!g->inCnt || !g->firstColumn || !g->F)  die(__func__);
  // init the above arrays
  if(g->smallAlpha) init_arrays(g);
  else init_arrays_largealpha(g);
  // we are now ready to give mmap advise
  #ifdef USE_MMAP_ADVISE
  if(g->mmapZ) { // advise on g->mergeColor
    for (int i = 0; i < g->sizeOfAlpha-1; ++i)
      madvise(g->mergeColor + g->firstColumn[i], (g->firstColumn[i+1]-g->firstColumn[i])*sizeof(palette), MADV_SEQUENTIAL);
    madvise(g->mergeColor + g->firstColumn[g->sizeOfAlpha-1], (g->mergeLen-g->firstColumn[g->sizeOfAlpha-1])*sizeof(palette), MADV_SEQUENTIAL);
  }
  #endif

  // init liquid block (containing list of allocated mem)
  liquidBlock *liquid = liquid_new(g); 
  // init list (on disk) of irrelevant blocks, initially empty 
  solidBlockFile *ibList = ibHead_new(g);

  // main loop
  customInt prefixLength = 1;      
  int lcpSize = POS_SIZE + BSIZE;   // number of bytes for each pos,lcp pair: see writeLcp()
  int round=0;
  bool merge_completed; 
  do {
    prefixLength+= 1;
    if(prefixLength>MAX_LCP_SIZE && g->lcpMerge) {fprintf(stderr,"LCP too large (use --lbytes=4)\n");exit(EXIT_FAILURE);}
    if(g->lcpCompute && prefixLength-2>MAX_LCP_SIZE) {fprintf(stderr,"LCP too large (2) (use --lbytes=4)\n");exit(EXIT_FAILURE);}
    bool mergeChanged = false; // the Z vector has changed in this iteration (used when g->bwtOnly)
    ibList->fout = gap_tmpfile(g->outPath);
    merge_completed=addCharToPrefix(ibList,liquid,prefixLength,&mergeChanged,round,g);
    if (g->verbose>1 && lastRound) {
      printf("Lcp: "CUSTOM_FORMAT". Memory peak/current: %.2lf/%.2lf bytes/symbol. ibList: %ju\n", 
           prefixLength-1, (double)malloc_count_peak()/g->mergeLen,
           (double)malloc_count_current()/g->mergeLen, (uintmax_t) ftello(ibList->fout));
      // also EOF are written to unsortedLcp file so percentages are not accurate     
      if(g->unsortedLcp) printf("   unsorted lcp values: %ju (%.2lf%%)\n",  
      (uintmax_t) ftello(g->unsortedLcp)/lcpSize, (double) 100*ftello(g->unsortedLcp)/(lcpSize*g->mergeLen));
    }
    round  = 1 - round; // change round parity
    if(g->bwtOnly && !mergeChanged) {
      if(g->verbose>1) puts("Gap bwt-only early termination");
      fclose(ibList->fout);
      break;
    }
    if(ibList->fin!=NULL) fclose(ibList->fin);
    rewind(ibList->fout);
    ibList->fin = ibList->fout;
  } while(!merge_completed);  // end main loop
  if(ibList->fin!=NULL) fclose(ibList->fin);

  if (g->verbose>0) {
    if(lastRound)
      printf("Merge completed (%d bwts). Mem: %zu peak, %zu current, %.2lf/%.2lf bytes/symbol\n", g->numBwt, malloc_count_peak(),
           malloc_count_current(), (double)malloc_count_peak()/g->mergeLen,
           (double)malloc_count_current()/g->mergeLen);
    else if(g->verbose>1)
      printf("Merge completed (%d bwts). Mem: %zu peak, %zu current\n", g->numBwt, malloc_count_peak(),
           malloc_count_current());
  }
  liquid_free(liquid);
  ibHead_free(ibList);
  if(g->extMem) close_bw_files(g);

  // computation complete, do the merging. The following call writes the
  // (possibly remapped) merged BWT back to g->bws[0]; and if lcpMerge==true the merged LCP to g->lcps[0] 
  mergeBWTandLCP(g,lastRound);
  if(g->lcpCompute) {
    assert(lastRound);
    // close lcp file (and merge them?)
    close_unsortedLCP_files(g);
    if(g->verbose>0) printf("Remind to run mergelcp to obtain the final LCP array\n"); 
  }

  // free B array 
  if(!g->lcpMerge) tba_free(g->bitB, g->mergeLen, g->mmapB);
  else free_B_array(g); 
  free(g->F); // last five arrays deallocated
  free(g->firstColumn); 
  free(g->inCnt);
  free_merge_arrays(g);
}


// ----- access to color files on disk: declared static because they are not used elsewhere

// open file fmergeColor for reading sequentially from 0 to mergeLen
// open one file pointer pointing at fnewMerge for each symbol except 0
// the fnewmerge pointers are also set to the correct position according to g->firstColumn 
static void open_merge_files(g_data *g) {
  assert(g->extMem);
  // mergeColor file for reading (Z in pseudocode)
  g->fmergeColor = fopen(g->merge_fname,"rb");
  if(!g->fmergeColor) die("merge_open");
  #ifndef NDEBUG
  customInt c[256] = {0}; // init to zero
  for(customInt i = 0; i < g->mergeLen; i++) {
    int col=fread_color(g->fmergeColor);
    col &= 0x7F;
    assert(col>=0 && col<g->numBwt);
    c[col]++;
  }
  for(int i=0;i<g->numBwt;i++) 
    assert(c[i]==g->bwtLen[i]);
  rewind(g->fmergeColor);
  #endif
  // alphaSize-1 newMergeColor cwriters for writing (newZ in pseudocode)
  int fd = open(g->newmerge_fname,O_WRONLY);
  if(fd == -1) die("new_merge_open");
  g->fnewMergeColor = malloc(g->sizeOfAlpha*sizeof(cwriter));
  if(g->fnewMergeColor==NULL) die("new_merge_alloc");
  g->fnewMergeColor[0].fd = -1; // invalid file descriptor
  for(int i=1; i< g->sizeOfAlpha; i++)
    cwriter_init(&g->fnewMergeColor[i],fd,COLOR_WBUFFER_SIZE, g->firstColumn[i]*sizeof(palette));
}

// use bws[] to make bwf[i] point at the beginning of bws[i]  
static void close_merge_files(g_data *g) {
  assert(g->extMem);
  // close new merge files
  for(int i=1; i< g->sizeOfAlpha; i++) 
    cwriter_close(&g->fnewMergeColor[i]);
  close(g->fnewMergeColor[1].fd); // close file   
  free(g->fnewMergeColor);
  // close merge file
  int e = fclose(g->fmergeColor);
  if(e!=0) die("merge_close");
}

// write a single color to f (that should be merge or newmerge)
static void fwrite_color(int b, FILE *f) {
  int e = fwrite(&b,sizeof(palette),1,f);
  if(e!=1) die(__func__);
}

static int fread_color(FILE *f) {
  int b=0;
  int e = fread(&b,sizeof(palette),1,f);
  if(e!=1) die(__func__);
  return b;
}
