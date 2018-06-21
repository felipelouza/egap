// mergegap with numBwt <= 16; merge and newMerge share in the same array g->mergeColor array
// supports bwtOnly and lcpMerge NOT lcpCompute and extMem
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
// init Z, newZ B, and first Column array using g->bwtOcc[i][j]
static void init_arrays16(g_data *g)
{
  assert(g->numBwt<=16);
  customInt i=0; // position inside Z newZ and B  

  for(int j=0;j<g->sizeOfAlpha;j++) {
    if(g->bwtOnly) tba_or_m(g->bitB,i,1);
    else g->blockBeginsAt[i]=1;  // start of symbol j, correct lcp is 0
    g->firstColumn[j] = i;  // symbol j starts at position i    
    for(int b=0;b<g->numBwt;b++) {
      for(customInt t=0;t<g->bwtOcc[b][j];t++) {
        if(j==0) { // zero chars are all different, Z are newZ do not change 
          if(g->bwtOnly) tba_or_m(g->bitB,i,1);
          else g->blockBeginsAt[i]=1;  
        }
        g->mergeColor[i++] = b*17;      // write b in both cur and next
      } // end for t
    } // end for b 
  } // end for j
  assert(i==g->mergeLen);
  // extra check on mergeColor, can be commented out
  #ifndef NDEBUG
  customInt cnt[MAX_NUMBER_OF_BWTS] = {0};
  for(i=0;i<g->mergeLen;i++) cnt[g->mergeColor[i]&0xF]++;
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
static void init_arrays16_largealpha(g_data *g)
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
  init_arrays16(g);
  for(int i=0;i<g->numBwt;i++)
    free(g->bwtOcc[i]);
  free(g->bwtOcc);
  g->bwtOcc=NULL;
}

// single iteration of the Gap algorithm
// input is head of the irrelevant list
#define get_mergeColor(k,round) ((round)? (g->mergeColor[(k)]>>4)&0x0F : (g->mergeColor[(k)]&0x0F))
#define set_mergeColor(k,c,round) (g->mergeColor[(k)] = (round)? \
        ( (g->mergeColor[(k)]&0xF0) | (c)) : ((g->mergeColor[(k)]&0x0F) | ((c)<<4))  )


static bool addCharToPrefix16(solidBlockFile *solidHead, liquidBlock *liquid, customInt prefixLength, bool *mergeChanged, const int round, g_data *g) {
  assert(prefixLength <= MAX_LCP_SIZE);
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
  int m = (round%2==1) ? 1 : 2;  // mask for the bitB array

  protoBlock cblock = {.mono = false};
  solidBlock *next = readBlock(solidHead); // first block
  solidBlock *last = NULL;
  
  for (k = 0; k < g->mergeLen; ) { 
    assert(next==NULL || k <= next->beginsAt);  // we did not pass next block
    assert(last==NULL || last->nextBlock == next); // last is the immediately preceeding block

    // check if we are entering a block, and if the block is at least 2 iterations old
    bool start_block, last_block_recent=true; // for k=0 a new block starts, so last id properly initialized
    if(!g->bwtOnly) {
      start_block = (g->blockBeginsAt[k]>0) && (g->blockBeginsAt[k] < prefixLength);
      if(start_block) last_block_recent = g->blockBeginsAt[k]>=prefixLength-1;
    }
    else {
      start_block = tba_block_test_set(g->bitB,k,m,&last_block_recent);    
    }
    if (start_block) {
      // check if the block we just left is not recent and monochrome 
      if(cblock.mono==true && !last_block_recent) {
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
        *mergeChanged=true;  // remember there is a difference from the previous ieration 
      // create new block?
      if (blockID[currentChar] != id) {
        if(g->bwtOnly) { // no lcp just mark B array 
          if(last_block_recent) tba_mark_if0(g->bitB,positionToUpdate,m);
        } else    // update lcp
          if(last_block_recent && g->blockBeginsAt[positionToUpdate]==0) // only 0 values in B are overwritten
            g->blockBeginsAt[positionToUpdate] = prefixLength;
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


// entry point for the gap bwt/lcp merging procedure with at most 16 input sequences 
// if lcpMerge==false blockBeginsAt is replaced by a bit array 
void gap16(g_data *g, bool lastRound) {

  // check we can really use 4 bits
  assert(g->numBwt<=16);
  if(g->numBwt<=8 && g->bwtOnly && !g->algorithm) return gap8(g,lastRound);
  if(g->verbose>0) puts("BWT merging with gap16");
  assert(!g->lcpCompute);
  assert(!g->extMem);
  // init local global vars
  check_g_data(g);
  // allocate and clear bit/int array B 
  if(!g->lcpMerge) g->bitB = tba_alloc(g->mergeLen, g->mmapB);
  else alloc0_B_array(g); // alloc and clear blockBeginsAt array
    
  // allocate Z (merge) Znew 
  alloc_merge_array(g);
  // allocate other useful arrays
  g->inCnt = malloc(g->numBwt*sizeof(customInt)); 
  g->firstColumn = malloc(g->sizeOfAlpha*sizeof(customInt));
  g->F = malloc(g->sizeOfAlpha*sizeof(customInt));
  if(!g->inCnt || !g->firstColumn || !g->F) die(__func__);
    
  // init the above arrays
  if(g->smallAlpha) init_arrays16(g);
  else init_arrays16_largealpha(g);  
  customInt prefixLength = 1;
  
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
      if(prefixLength>MAX_LCP_SIZE && !g->bwtOnly) {fprintf(stderr,"LCP too large\n");die(__func__);}
      bool mergeChanged = false; // the Z vector has changed in this iteration (used when g->bwtOnly)
      ibList->fout = gap_tmpfile(g->outPath);
      merge_completed=addCharToPrefix16(ibList,liquid,prefixLength,&mergeChanged,round,g);
      if (g->verbose>1 && lastRound) {
        printf("Lcp: "CUSTOM_FORMAT". Memory: %zu peak, %zu current, %.4lf/%.4lf bytes/symbol\n", 
             prefixLength-1, malloc_count_peak(),
             malloc_count_current(), (double)malloc_count_peak()/g->mergeLen,
             (double)malloc_count_current()/g->mergeLen);
      }
      if(g->bwtOnly && !mergeChanged) {
        if(g->verbose>1) puts("Gap bwt-only early termination");
        fclose(ibList->fout);
        break;
      }
      round  = 1 - round; // change round parity
      if(ibList->fin!=NULL) fclose(ibList->fin);
      rewind(ibList->fout);
      ibList->fin = ibList->fout;
    } while(!merge_completed);  // end main loop
    if(ibList->fin!=NULL) fclose(ibList->fin);
  }
  if (g->verbose>0) {
    if(lastRound)
      printf("Merge16 completed (%d bwts). Mem: %zu peak, %zu current, %.2lf/%.2lf bytes/symbol\n", g->numBwt, malloc_count_peak(),
           malloc_count_current(), (double)malloc_count_peak()/g->mergeLen,
           (double)malloc_count_current()/g->mergeLen);
    else if(g->verbose>1)
      printf("Merge16 completed (%d bwts). Mem: %zu peak, %zu current\n", g->numBwt, malloc_count_peak(),
           malloc_count_current());
  }
  liquid_free(liquid);
  ibHead_free(ibList);  

  // normalize mergeColor
  if(round!=0) for(customInt i=0;i<g->mergeLen;i++) g->mergeColor[i] >>=4;
  else for(customInt i=0;i<g->mergeLen;i++) g->mergeColor[i] &= 0x0F;
  // computation complete, do the merging. The following call writes the
  // (possibly remapped) merged BWT back to g->bws[0]; and if lcpMerge==true the merged LCP to g->lcps[0]     
  mergeBWTandLCP(g,lastRound);

  // free B array 
  if(!g->lcpMerge) tba_free(g->bitB, g->mergeLen, g->mmapB);
  else free_B_array(g); 
  free(g->F); // last five arrays deallocated
  free(g->firstColumn); 
  free(g->inCnt);
  free_merge_array(g);// contains merge and newMerge
}

#undef get_mergeColor
#undef set_mergeColor
