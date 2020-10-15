
// mergegap with numBwt <= 8 and B represented as a 2bit array (lcpMerge==false)
// since the sizes of Z+newZ+B = 3+3+2 = 8 we use g->mergeColor to store them all  
// supports bwtOnly & lcpCompute NOT lcpMerge
// external memory: BWTs 

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
 * special update rule, we init that region also in Znew (newMergeColor)
 * and we never modify thet region in the algorithm.   
 * The array firstColumn (compact representation of F) is also initialized
 *  
 * Note, when computing only the BWT instead of B we init bitB which
 * uses 2 bits per entry to encode the values: 
 *   never set->00,  recently set->01 or 10,  set at least 2 iterations before->11
 * during this initialization we write 01 for recently set entries,
 * therefore in the first iteration the mask for access to bitB should be 10 (eg 2)  
 * */     
// init Z, newZ B, and first Column array using g->bwtOcc[i][j]
static void init_arrays8(g_data *g)
{
  assert(g->numBwt<=8);
  customInt i=0; // position inside Z newZ and B  

  for(int j=0;j<g->sizeOfAlpha;j++) {
    tba_or_m8(g->mergeColor,i,1);
    g->firstColumn[j] = i;  // symbol j starts at position i
    for(int b=0;b<g->numBwt;b++) {
      for(customInt t=0;t<g->bwtOcc[b][j];t++) {
        if(j==0)   // zero chars are all different, Z are newZ do not change 
          tba_or_m8(g->mergeColor,i,1);
        g->mergeColor[i++] |= b*9;      // write b in both cur and next
      } // end for t
    } // end for b 
  } // end for j
  assert(i==g->mergeLen); 
  // extra check on mergeColor, can be commented out
  #ifndef NDEBUG
  customInt cnt[MAX_NUMBER_OF_BWTS] = {0};
  for(i=0;i<g->mergeLen;i++) cnt[g->mergeColor[i]&0x7]++;
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
static void init_arrays8_largealpha(g_data *g)
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
  init_arrays8(g);
  for(int i=0;i<g->numBwt;i++)
    free(g->bwtOcc[i]);
  free(g->bwtOcc);
  g->bwtOcc=NULL;
} 

// single iteration of the Gap algorithm
// input is head of the irrelevant lists (fin and fout) and an empty liquid block
// return true if the whole sequence has become irrelevant.  

// macros to access colors: one color is in the last 3 bits, the other is in bits 4-6
#define get_mergeColor8(k,round) ((round)? (g->mergeColor[(k)]>>3)&0x07 : (g->mergeColor[(k)]&0x07))
#define set_mergeColor8(k,c,round) (g->mergeColor[(k)] = (round)? \
        ( (g->mergeColor[(k)]&0xF8) | (c)) : ((g->mergeColor[(k)]&0xC7) | ((c)<<3))  )

static bool addCharToPrefix8(solidBlockFile *solidHead, liquidBlock *liquid, uint32_t prefixLength, bool *mergeChanged, const int round, g_data *g) {
  assert(liquid->empty);
  liquid->beginsAt = liquid->endsAt = 0;  
  for(int i=0;i<liquid->occ_size;i++) assert(liquid->occ[i]==0);
  // copy first column to F 
  array_copy(g->F, g->firstColumn, g->sizeOfAlpha); //initialize char positions
  // pointer inside each BWT (k_0 & k_1 in the pseudocode)
  array_clear(g->inCnt,g->numBwt,0);
  if(g->extMem) rewind_bw_files(g); // set file pointers at the beginning of each BWT 
  // id for each character, init with an invalid id
  customInt blockID[g->sizeOfAlpha];
  array_clear(blockID,g->sizeOfAlpha, g->mergeLen);  // mergeLen is an invalid id
  customInt id = 0, k; 
  int m = (round%2==1) ? 1 : 2; // mask for the bitB array (now inside mergeColor)
  uint64_t lcpWritten =0;

  protoBlock cblock = {.mono = false};
  solidBlock *next = readBlock(solidHead); // first block
  solidBlock *last = NULL;
  
  for (k = 0; k < g->mergeLen; ) { 
    assert(next==NULL || k <= next->beginsAt);  // we did not pass next block
    assert(last==NULL || last->nextBlock == next); // last is the immediately preceeding block

    // check if we are entering a block, and if the block is at least 2 iterations old
    bool start_block, last_block_recent=true; // for k=0 a new block starts, so last id properly initialized
    start_block = tba_block_test_set8(g->mergeColor,k,m,&last_block_recent);    
    if(start_block && last_block_recent && g->lcpCompute)
      {writeLcp(k,prefixLength-2,g); lcpWritten++;} // save lcp value found in previous iteration
    if (start_block) {
      // if the block we just left is a singleton we add it to liquid that remains active
      if(!last_block_recent && cblock.mono==true && (cblock.beginsAt==k-1)) {
        cblock.endsAt = k;           // solidifiable singleton block just ended
        add_singleton2liquid(&cblock, liquid);
      }
      // if the block we left is not recent, monochrome, we are only computing BWTs and not using extMem add it  
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
        cblock.beginsAt = k; cblock.mono=true; 
        if(!g->extMem) {
          cblock.color =  get_mergeColor8(k,round);
          cblock.start = &g->bws[cblock.color][g->inCnt[cblock.color]]; // bwt-position of first char in block
        }
      }
      else cblock.mono = false; // not a candidate for solid block, wait next iteration 
      if(last_block_recent) 
         id = k;    // id of the new block 
    }
    // processing a char in a relevant block
    int currentColor = get_mergeColor8(k,round);   // g->mergeColor[k] b in pseudocode
    int currentChar=0;  // c in pseudocode: next char in current BWT 
    if(g->extMem) {
      int e = fread(&currentChar,sizeof(symbol),1,g->bwf[currentColor]);
      if(e!=1) {die(__func__);} g->inCnt[currentColor]++;
    }
    else 
      currentChar =  g->bws[currentColor][g->inCnt[currentColor]++]; // c in pseudocode
    // add currentChar/Color to proto block  
    cblock.lastChar =  currentChar;   // save lastchar, only useful for singleton blocks
    cblock.lastColor =  currentColor; // save lastcolor, only useful for singleton blocks
    if(!g->extMem)
      if(currentColor != cblock.color) cblock.mono = false;       // block is not monochrome
    // write color in new Z array, except 0 chars
    if(currentChar!=0) {
      customInt positionToUpdate = g->F[currentChar]++;
      set_mergeColor8(positionToUpdate,currentColor,round); // g->newMergeColor[positionToUpdate] = currentColor;
      if(g->bwtOnly && !*mergeChanged && get_mergeColor8(positionToUpdate,round)!=currentColor)
        *mergeChanged=true;  // remember there is a difference from the previous iteration 
      // create new block?
      if (blockID[currentChar] != id) {
        if(last_block_recent) tba_mark8_if0(g->mergeColor,positionToUpdate,m);
        blockID[currentChar] = id; // update block id, always!
      }
    }
    k++;
  } // end main loop, check block that just ended  
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
  if(!g->bwtOnly && lcpWritten>0)
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

  return everything_irrelevant;
}


// entry point for the gap bwt/lcp merging procedure with at most 8 input sequences 
// we assume lcpMerge==false so blockBeginsAt is replaced by a bit array 
void gap8(g_data *g, bool lastRound) {
  // check we can really use 4 bits
  assert(g->numBwt<=8 && g->lcpMerge==false);
  if(g->verbose>0) puts("BWT merging with gap8");
  if(g->lcpCompute) {       // we compute LCP values only if we are at the last round 
    assert(!g->bwtOnly && lastRound);
    open_unsortedLCP_files(g);
    if(g->verbose>0) puts("Computing LCP values");
  }
  else assert(g->bwtOnly); 

  // init local global vars
  check_g_data(g);    
  if(g->extMem) open_bw_files(g);
  // allocate an array containing Z (merge) Znew and B 
  // together they use 8 bits 1+1+3+3 
  assert(g->blockBeginsAt==NULL);
  alloc_merge_array(g);

  // allocate other useful arrays
  g->inCnt = malloc(g->numBwt*sizeof(customInt)); 
  g->firstColumn = malloc(g->sizeOfAlpha*sizeof(customInt));
  g->F = malloc(g->sizeOfAlpha*sizeof(customInt));
  if(!g->inCnt || !g->firstColumn || !g->F) die(__func__);
  // init the above arrays
  if(g->smallAlpha) init_arrays8(g);
  else init_arrays8_largealpha(g);
  // now we are ready to give mmap advise
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
  uint32_t prefixLength = 1;
  int lcpSize = POS_SIZE + BSIZE;   // number of bytes for each pos,lcp pair, see writeLcp()
  int round=0;
  bool merge_completed;
  do {
    prefixLength+= 1;
    if(g->lcpCompute && prefixLength-2>MAX_LCP_SIZE) {fprintf(stderr,"LCP too large: %u\n", prefixLength-2);die(__func__);}
    bool mergeChanged = false; // the Z vector has changed in this iteration (used when g->bwtOnly)
    ibList->fout = gap_tmpfile(g->outPath);
    merge_completed=addCharToPrefix8(ibList,liquid,prefixLength,&mergeChanged,round,g);
    if (g->verbose>1 && lastRound) {
      #if MALLOC_COUNT_FLAG
        printf("Lcp: %u. Memory peak/current: %.2lf/%.2lf bytes/symbol. ibList: %ju\n", 
           prefixLength-1, (double)malloc_count_peak()/g->mergeLen,
           (double)malloc_count_current()/g->mergeLen, (uintmax_t) ftello(ibList->fout));
      // also EOF are written to unsortedLcp file so percentages are not accurate     
      #endif
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

  #if MALLOC_COUNT_FLAG
    if (g->verbose>0) {
      if(lastRound)
        printf("Merge8 completed (%d bwts). Mem: %zu peak, %zu current, %.2lf/%.2lf bytes/symbol\n", g->numBwt, malloc_count_peak(),
             malloc_count_current(), (double)malloc_count_peak()/g->mergeLen,
             (double)malloc_count_current()/g->mergeLen);
      else if(g->verbose>1)
        printf("Merge8 completed (%d bwts). Mem: %zu peak, %zu current\n", g->numBwt, malloc_count_peak(),
             malloc_count_current());
    }
  #else
    if (g->verbose>0) {
        printf("Merge8 completed (%d bwts).\n", g->numBwt);
    }
  #endif
  liquid_free(liquid);
  ibHead_free(ibList);  
  if(g->extMem) close_bw_files(g);

  // computation complete, do the merging. 
  // The following call writes the merged BWT back to g->bws[0] (if extMem to the BWTs file)
  mergeBWT8(g,lastRound);  
  if(g->lcpCompute) {
    assert(lastRound);
    // close lcp file (and merge them?)
    close_unsortedLCP_files(g);
    if(g->verbose>0) printf("Remind to run lcpmerge to obtain the final LCP array\n"); 
  }

  free(g->F); // last five arrays deallocated
  free(g->firstColumn); 
  free(g->inCnt);
  free_merge_array(g);
}

#undef get_mergeColor8
#undef set_mergeColor8
