// mergegap with numBwt <= 256; merge newMerge and blockBeginsAt share 
// the same array g->array32. The rational is that these arrays
// have the same access pattern. 
// supports bwtOnly (inefficient, better gap or gap128) and
//          lcpCompute (without external mergesort)
// does not support lcpMerge and extMem 

/**
 * Using the number of occs of each symbol in each bwt (stored in bwtOcc)
 * init the array Z (mergeColor) and B (blockBeginsAt) at the value
 * they should have after the first iteration of the Gap algorithm:
 *   blockBeginsAt[i]=1 if i is the first occurrence in the first
 *                      column F of a new symbol or F[i]=0
 *                      (0 occurrences are assumed to be all different)
 *   in each region of the F column with the same symbol j in Z
 *   we have: #occ(j) in bwt[0], #occ(j) in bwt(1), and so on
 * Since the region corresponding to 0 does not chage and 0 has a 
 * special udate rule, we init that region also in Znew (newMergeColor)
 * and we never modify thet region in the algorithm.   
 * The array firstColumn (compact representation of F) is also initialized
 * 
 * Note, when computing only the BWT instead of B we init bitB which
 * uses 2 bits per entry to encode the values: 
 *   never set->00,  recently set->01 or 10,  set at least 2 iterations before->11
 * during this initialization we write 01 for recently set entries,
 * therefore in the first iteration the mask for access to bitB should be 10 (eg 2)  
 * */     


// Note g->array32 stores merge, newMerge, blockBeginsAt in
// 16 + 8 + 8 + bit format 


// macros to access colors and B-values inside g->array32
// one color is in the last 8 bits, the other in bits 8-15
// round is 0 or 1. the current color is in bits [round*8,(round+1)*8-1]
// B-values are in bits 16-31
#define get_mergeColor(k,round) ((round)? (g->array32[(k)]>>8)&0xFF : (g->array32[(k)]&0xFF))
#define set_mergeColor(k,c,round) (g->array32[(k)] = (round)? \
        ( (g->array32[(k)]&0xFFFFFF00) | (c)) : ((g->array32[(k)]&0xFFFF00FF) | ((c)<<8))  )

#define get_blockBeginsAt(k) ((g->array32[(k)]>>16)&0xFFFF)
#define set_blockBeginsAt(k,c) (g->array32[(k)] = \
        ( (g->array32[(k)]& 0xFFFF) | ((c)<<16) ))


// init Z, newZ B, and first Column array using g->bwtOcc[i][j]
static void init_arrays256(g_data *g)
{
  assert(g->numBwt<=256);
  customInt i=0; // position inside Z newZ and B  

  for(int j=0;j<g->sizeOfAlpha;j++) {
    set_blockBeginsAt(i,1);   // start of symbol j, correct lcp is 0
    g->firstColumn[j] = i;    // symbol j starts at position i    
    for(int b=0;b<g->numBwt;b++) {
      for(customInt t=0;t<g->bwtOcc[b][j];t++) {
        if(j==0)  // zero chars are all different, Z are newZ do not change 
          set_blockBeginsAt(i,1);  
        g->array32[i++] |= b*257;      // write b in both cur and next
      } // end for t
    } // end for b 
  } // end for j
  assert(i==g->mergeLen);
  // extra check on mergeColor, can be commented out
  #ifndef NDEBUG
  customInt cnt[MAX_NUMBER_OF_BWTS] = {0};
  for(i=0;i<g->mergeLen;i++) cnt[g->array32[i]&0xFF]++;
  bool stop=false;
  for(int i=0; i<g->numBwt; i++)
    if(cnt[i]!=g->bwtLen[i]) {
      printf("INIT %d cnt:"CUSTOM_FORMAT" len:"CUSTOM_FORMAT"\n",i,cnt[i],g->bwtLen[i]);
      stop=true;
    }
  assert(!stop);
  #endif
}


// init Z, newZ and B array without using g->bwtOcc[i][j]
static void init_arrays256_largealpha(g_data *g)
{
  assert(!g->extMem);                 // if extMem we do not want to reread the text
  // compute bwtOcc on the spot with a complete scan of input BWTs
  assert(g->bwtOcc==NULL);
  g->bwtOcc = malloc(g->numBwt*sizeof(customInt *));
  if(!g->bwtOcc) die(__func__);
  for(int i=0;i<g->numBwt;i++) {
    g->bwtOcc[i] = calloc(g->sizeOfAlpha,sizeof(customInt));
    if(!g->bwtOcc[i]) die(__func__);
    init_freq_no0(g->bws[i],g->bwtLen[i],g->bwtOcc[i]);
  }
  init_arrays256(g);
  for(int i=0;i<g->numBwt;i++)
    free(g->bwtOcc[i]);
  free(g->bwtOcc);
  g->bwtOcc=NULL;
}


// single iteration of the Gap algorithm
// input is head of the irrelevant lists (fin and fout) and an empty liquid block
// return true if the whole sequence has become irrelevant.
static bool addCharToPrefix256(solidBlockFile *solidHead, liquidBlock *liquid, uint32_t prefixLength, bool *mergeChanged, const int round, g_data *g) {
  assert(prefixLength <= 65535);// lengths are stored in 16 bits in g->array32
  assert(liquid->empty);
  liquid->beginsAt = liquid->endsAt = 0;  
  for(int i=0;i<liquid->occ_size;i++) assert(liquid->occ[i]==0);
  // copy first column to F 
  array_copy(g->F, g->firstColumn, g->sizeOfAlpha); //initialize char positions
  // pointer inside each BWT (k_0 & k_1 in the pseudocode)
  array_clear(g->inCnt,g->numBwt,0);
  // id for each character, init with an invalid id
  customInt blockID[g->sizeOfAlpha];
  array_clear(blockID,g->sizeOfAlpha, g->mergeLen);  // mergeLen is an invalid id
  customInt id = 0, k; 

  protoBlock cblock = {.mono = false};
  solidBlock *next = readBlock(solidHead); // first block
  solidBlock *last = NULL;
  
  for (k = 0; k < g->mergeLen; ) { 
    assert(next==NULL || k <= next->beginsAt);  // we did not pass next block
    assert(last==NULL || last->nextBlock == next); // last is the immediately preceeding block

    // check if we are entering a block, and if the block is at least 2 iterations old
    bool start_block, last_block_recent=true; // for k=0 a new block starts, so last id properly initialized
    int bk = get_blockBeginsAt(k);
    start_block = (bk>0) && (bk < prefixLength);// if bk==prefixLength the value has been written during this iteration 
    if(start_block) last_block_recent = (bk >=prefixLength-1);
    if (start_block) {
      // check if the block we just left is not recent and monochrome and if(unsortedLcp) singleton  
      if( !last_block_recent && cblock.mono && ((cblock.beginsAt==k-1)||g->bwtOnly)) {
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
        else  //three way merge, next is freed last does not change
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
        cblock.beginsAt = k; cblock.mono=true; cblock.color =  get_mergeColor(k,round);     // g->mergeColor[k];
        cblock.start = &g->bws[cblock.color][g->inCnt[cblock.color]]; // bwt-position of first char in block
      }
      else cblock.mono = false; // not a candidate for solid block, wait next iteration 
      if(last_block_recent) 
         id = k;    // id of the new block 
    }
    // processing a char in a relevant block
    int currentColor = get_mergeColor(k,round);   // g->mergeColor[k] b in pseudocode
    int currentChar = g->bws[currentColor][g->inCnt[currentColor]++]; // c in pseudocode
    if(currentColor != cblock.color) cblock.mono = false;       // block is not monochrome
    // write color in new Z array, except 0 chars
    if(currentChar!=0) {
      customInt positionToUpdate = g->F[currentChar]++;
      set_mergeColor(positionToUpdate,currentColor,round); // g->newMergeColor[positionToUpdate] = currentColor;
      if(g->bwtOnly && !*mergeChanged && get_mergeColor(positionToUpdate,round)!=currentColor)
        *mergeChanged=true;  // remember there is a difference from the previous iteration 
      // create new block?
      if (blockID[currentChar] != id) {
        if(last_block_recent && get_blockBeginsAt(positionToUpdate)==0) // only 0 values in B are overwritten
           set_blockBeginsAt(positionToUpdate, prefixLength);
        blockID[currentChar] = id; // update block id, always!
      }
    }
    k++;
  } // end main loop
  assert(next==NULL); 
  if(cblock.mono==true) { // final mono block can be created 
    cblock.endsAt = k;
    add_proto2liquid(&cblock,liquid);  // add proto to liquid that remains active
    assert(!liquid->empty);
  }
  if(!liquid->empty) 
    last = finalize_liquid(last,liquid,NULL,solidHead); // a new block could be created 
  assert(liquid->empty);
  // check if all sequence has become irrelevant
  bool everything_irrelevant = false;
  if(last!=NULL && last->beginsAt==0 && last->endsAt==g->mergeLen)
    everything_irrelevant = true;
  // save last block   
  if(last!=NULL) writeBlock(last,solidHead);
  // check we have read all chars from all BWT's 
  for(int i=0; i<g->numBwt; i++)
    assert(g->inCnt[i]==g->bwtLen[i]);

  return everything_irrelevant;
}


// entry point for the gap bwt/lcp merging procedure
// even when we are only interested in the BWT we use blockBeginsAt (here embedded in array32) 
// to keep track of blocks. 
// This function is specialized in computing the LCP from scratch, but without using 
// external memory mergesort (at the cost of using a 16-bit B array instead of a 2-bit array as in gap128)
// Note that gap128 is probably a better alternative unless B can be stored in RAM
// when the function returns the merged BWT is in g->bws[0] while the LCP values are in the output file 
void gap256(g_data *g, bool lastRound) 
{
  assert(g->lcpCompute && !g->mwXMerge); // does not make sense to be here unless we want compute LCP from scratch  with multiway mergesort 
  // extra check we can really use 8 bits for Z and newZ and we are not interested in merging 
  assert(g->numBwt<=256 &&  !g->lcpMerge);
  if(g->verbose>0) puts("BWT merging with gap256");
  assert(sizeof(lcpInt) <= 2);
  assert(!g->lcpMerge);
  assert(!g->extMem);
  if(g->lcpCompute)       // we compute LCP values only if we are at the last round 
    assert(!g->bwtOnly && lastRound);
  else assert(g->bwtOnly); 

  // init local global vars
  check_g_data(g);
  // allocate Z (merge) Znew and B and clear them
  alloc_array32(g);
  
  // allocate other useful arrays
  g->inCnt = malloc(g->numBwt*sizeof(customInt)); 
  g->firstColumn = malloc(g->sizeOfAlpha*sizeof(customInt));
  g->F = malloc(g->sizeOfAlpha*sizeof(customInt));
  if(!g->inCnt || !g->firstColumn || !g->F)  die(__func__);
  // init the above arrays
  if(g->smallAlpha) init_arrays256(g);
  else init_arrays256_largealpha(g);  
  uint32_t prefixLength = 1;
  
  #ifdef USE_MMAP_ADVISE
  // advise on g->array32
  for (int i = 0; i < g->sizeOfAlpha-1; ++i) {
    madvise(g->array32 + g->firstColumn[i], (g->firstColumn[i+1]-g->firstColumn[i])*4, MADV_SEQUENTIAL);
  }
  madvise(g->array32 + g->firstColumn[g->sizeOfAlpha-1], (g->mergeLen-g->firstColumn[g->sizeOfAlpha-1])*4, MADV_SEQUENTIAL);
  #endif

  // init liquid block (containing list of allocated mem)
  liquidBlock *liquid = liquid_new(g); 
  // init list (on disk) of irrelevant blocks, initially empty 
  solidBlockFile *ibList = ibHead_new(g);
  // main loop
  int round=0;
  if(g->numBwt>1) {
    bool merge_completed; 
    do {
      prefixLength+= 1;
      if(prefixLength> 0xFFFF ) {fprintf(stderr,"prefixLength too large: %u\n", prefixLength);die(__func__);}
      bool mergeChanged = false; // the Z vector has changed in this iteration (used when g->bwtOnly)
      ibList->fout = gap_tmpfile(g->outPath);
      merge_completed=addCharToPrefix256(ibList,liquid,prefixLength,&mergeChanged,round,g);
      #if MALLOC_COUNT_FLAG
        if (g->verbose>1 && lastRound) {
          printf("Lcp: %u. Memory peak/current: %.2lf/%.2lf bytes/symbol. ibList: %ju\n", 
               prefixLength-1, (double)malloc_count_peak()/g->mergeLen,
               (double)malloc_count_current()/g->mergeLen, (uintmax_t) ftello(ibList->fout));
        }
      #endif
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
  }
  if (g->verbose>0) {
    #if MALLOC_COUNT_FLAG
      if(lastRound)
        printf("Merge256 completed (%d bwts). Mem: %zu peak, %zu current, %.2lf/%.2lf bytes/symbol\n", g->numBwt, malloc_count_peak(),
             malloc_count_current(), (double)malloc_count_peak()/g->mergeLen,
             (double)malloc_count_current()/g->mergeLen);
      else if(g->verbose>1)
        printf("Merge256 completed (%d bwts). Mem: %zu peak, %zu current\n", g->numBwt, malloc_count_peak(),
             malloc_count_current());
    #else
      printf("Merge256 completed (%d bwts).\n", g->numBwt);
    #endif
  }
  liquid_free(liquid);
  ibHead_free(ibList);

  // computation complete, do the merging. Write the merged BWT to g->array32
  // the LCP if requested is in the high 16 bits of g->array32 and is written to the output file  
  mergeBWTandLCP256(g,lastRound);

  free(g->F); // last five arrays deallocated
  free(g->firstColumn); 
  free(g->inCnt);
  free_array32(g);
  assert(g->blockBeginsAt==NULL); // lcp values are not here!
}

#undef get_blockBeginsAt
#undef set_blockBeginsAt
#undef get_mergeColor
#undef set_mergeColor
