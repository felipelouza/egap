// external memory version only
// using merge/mergecolor to store up to 128 colors + 1 bit 

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
 * Note, that since we are computing only the  WT instead of B we use a 
 * virtual array bitB containing 2 bits per entry with the following meaning:
 *   never set->00,  recently set->01 or 10,  set at least 2 iterations before->11
 * during this initialization we write 01 for recently set entries,
 * therefore in the first iteration the mask for access to bitB should be 10 (eg 2)  
 * */     
// init Z, newZ B, and first Column array using g->bwtOcc[i][j]
static void init_arrays128ext(g_data *g)
{
  assert(g->numBwt>0 && g->numBwt<128);
  customInt i=0; // position inside Z newZ and B  
  // open merge files 
  FILE *fnewmerge=NULL, *fmerge=NULL;
  fnewmerge = fopen(g->newmerge_fname,"wb");
  if(!fnewmerge) die("mergegap:init_arrays:fnewmerge open");
  fmerge = fopen(g->merge_fname,"wb");
  if(!fmerge) die("mergegap:init_arrays:fmerge open");
  // scan bwtOcc
  for(int j=0;j<g->sizeOfAlpha;j++) {
    bool bit=1; // printf("1 at %ld\n",i);  // start of symbol j
    g->firstColumn[j] = i;  // symbol j starts at position i
    for(int b=0;b<g->numBwt;b++) {
      for(customInt t=0;t<g->bwtOcc[b][j];t++) { 
        if(j==0) { // zero chars are all different, Z, Znew never change 
          bit=1;
          fwrite_color(b,fnewmerge); //for 0-chars write b also to newmerge aka newZ 
        }
        // all colors are written to merge aka Z
        fwrite_color(b|(bit<<7),fmerge);
        i++; bit=0;
      } // end for t
    } // end for b 
  } // end for j
  assert(i==g->mergeLen); 
  // extra check on mergeColor, can be commented out
  assert(ftell(fmerge)==g->mergeLen*sizeof(palette));
  if(fclose(fmerge)!=0) die("init_arrays:fmerge close"); 
  if(fclose(fnewmerge)!=0) die("mergegap:init_arrays:fnewmerge close"); 
}


// single iteration of the Gap algorithm
// input is head of the irrelevant lists (fin and fout) and an empty liquid block
// return true if the whole sequence has become irrelevant.  

static bool addCharToPrefix128ext(solidBlockFile *solidHead, liquidBlock *liquid, uint32_t prefixLength, bitfile *b, g_data *g) {
  assert(liquid->empty);
  liquid->beginsAt = liquid->endsAt = 0;  
  for(int i=0;i<liquid->occ_size;i++) assert(liquid->occ[i]==0);
  // copy first column to F 
  array_copy(g->F, g->firstColumn, g->sizeOfAlpha); //initialize char positions
  // pointer inside each BWT (k_0 & k_1 in the pseudocode)
  array_clear(g->inCnt,g->numBwt,0);
  rewind_bw_files(g);  // set file pointers at the beginning of each BWT 
  bitfile_rewind(b);    
  open_merge_files(g); // open merge file for reading and newmerge files for writing
  // id for each character, init with an invalid id
  customInt blockID[g->sizeOfAlpha];
  array_clear(blockID,g->sizeOfAlpha, g->mergeLen);  // mergeLen is an invalid id
  customInt id = 0, k; 
  uint64_t lcpWritten =0;

  protoBlock cblock = {.beginsAt = 0};
  solidBlock *next = readBlock(solidHead); // first block
  solidBlock *last = NULL;
  
  for (k = 0; k < g->mergeLen; ) { 
    assert(next==NULL || k <= next->beginsAt);  // we did not pass next block
    assert(last==NULL || last->nextBlock == next); // last is the immediately preceeding block
    // read newblock & color
    int currentColor = fread_color(g->fmergeColor);// read color and new block bit 
    bool new_block = ((currentColor & 0x80)!=0);   // extract new block bit 
    currentColor &= 0x7F;                          // delete new block bit from color
    // read the old block bit: it is set if the block is at least 2 iterations old
    assert(bitfile_tell(b)==k);
    bool old_block = bitfile_read_or_write(b,new_block);
    if(new_block && !old_block && g->lcpCompute)
      {writeLcp(k,prefixLength-2,g); lcpWritten++;} // save lcp value found in previous iteration
    if (old_block) {
      // if the block we just left is a singleton we add it to liquid that remains active
      if(cblock.beginsAt==k-1) {
        cblock.endsAt = k;           // solidifiable singleton block just ended
        add_singleton2liquid(&cblock, liquid);
      }
      else { // proto block cannot be added, close current liquid
        if(!liquid->empty)
          last = finalize_liquid(last,liquid,next,solidHead); // this is the only point where a new solid block is created
        assert(liquid->empty);
        liquid->beginsAt=liquid->endsAt=k; // start empty liquid block
      }
      assert(liquid->endsAt==k);
      // block ending at k considered, now look forward 
      if(next!=NULL && next->beginsAt==k) { // entering an irrelevant block
        skip128ext(next, b, g);             // skip block
        k = next->endsAt;                   // update k (note after this there is a continue)
        // merge liquid with next block and possibly previous 
        if(last==NULL || last->endsAt!=liquid->beginsAt) {
            if(!liquid->empty) merge_liquid(liquid,next,solidHead); // simple merge 
            if(last!=NULL) writeBlock(last,solidHead);  // save current last  
            last = next;   // advance last           
        } 
        else  //three way merge: next is freed, last does not change
          merge_sls(last,liquid,next,solidHead); // only point where a solid block can be destroyed 
        assert(liquid->empty);
        liquid->beginsAt=liquid->endsAt=k; // start empty liquid block
        next = readBlock(solidHead);       // next has become last, update next (was: next = last->nextBlock; )
        last->nextBlock = next;
        assert(k==last->endsAt);
        continue;                   // resume from the end of the block 
      }  // end skip solid block
      // we are entering a relevant block: make it a candidate for solidification  
      cblock.beginsAt = k; 
    }   // end if(new_block || old_block)
    // processing a char in a relevant block
    assert(ftell(g->fmergeColor)==(k+1)*sizeof(palette));   // we already have currentColor but we check the file pointer is at the right position  
    assert(ftell(g->bwf[currentColor])==(g->inCnt[currentColor]+(g->bws[currentColor]-g->bws[0])+g->symb_offset)*sizeof(symbol)); 
    int currentChar=0;
    int e = fread(&currentChar,sizeof(symbol),1,g->bwf[currentColor]);
    if(e!=1) die("mergegap:addCharToPrefix:bwt[color]Read"); 
    g->inCnt[currentColor]++;
    // add currentChar/Color to proto block 
    cblock.lastChar =  currentChar;   // save lastchar, only useful for singleton blocks
    cblock.lastColor =  currentColor; // save lastcolor, only useful for singleton blocks

    // write color and newblock bit in newMerge array, except 0 chars
    if(new_block) id = k;
    if(currentChar!=0) {
      int bit = (blockID[currentChar] != id) ? 1 : 0;
      cwriter_put(&g->fnewMergeColor[currentChar],currentColor|(bit<<7));
      g->F[currentChar]++; 
      assert(cwriter_tell(&g->fnewMergeColor[currentChar])==g->F[currentChar]*sizeof(palette));
      blockID[currentChar] = id; // update block id, always!
    }
    k++;
  } // end main loop
  assert(bitfile_tell(b)==g->mergeLen);
  bitfile_flush(b);    // save to file last bits 
  assert(k==g->mergeLen);
  assert(next==NULL); 
  if(cblock.beginsAt==k-1) {
    cblock.endsAt = k;           // solidifiable singleton block just ended
    add_singleton2liquid(&cblock, liquid);
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
  assert(ftell(g->fmergeColor)==g->mergeLen*sizeof(palette));
  assert(g->F[g->sizeOfAlpha-1]==g->mergeLen);
  assert(cwriter_tell(&g->fnewMergeColor[g->sizeOfAlpha-1])==g->mergeLen*sizeof(palette));
  close_merge_files(g);
  char *tmp=g->merge_fname; g->merge_fname = g->newmerge_fname; g->newmerge_fname = tmp;
  return everything_irrelevant;
}


/* The main loop on the following function is based on the variable prefix_length
 * 
 * prefix length is the position in the context we are considering:
 * prefix_lenght==1 is the first char in the prefix, and we consider it
 * by counting the characters in the input BWT. The first actual iteration 
 * is done with prefix length==2 which means we are looking at the second 
 * character in the prefix. If two suffixes are in different blocks at
 * the begining of iteration with prefix_length=i this means their contexts differ
 * in one of the positions 1,2,...i-1 and if they are put for the first time 
 * in different block during iteration with prefix_length=i, theier context 
 * differ for the first time in position i, hence their LCP is i-1 */




// entry point for the gap bwt/lcp merging procedure with at most 128 input sequences 
// we assume lcpMerge==false so blockBeginsAt is replaced by a bit array 
void gap128ext(g_data *g, bool lastRound) {
  // we are specialized for external memory
  assert(g->extMem); 
  // check we can really use 7 bits
  assert(g->numBwt<=128 && g->lcpMerge==false);
  if(g->verbose>0) puts("BWT merging with gap128ext");
  if(g->lcpCompute) {       // we compute LCP values only if we are at the last round 
    assert(!g->bwtOnly && lastRound);
    open_unsortedLCP_files(g);
    if(g->verbose>0) puts("Computing LCP values");
  }
  else assert(g->bwtOnly); 

  check_g_data(g);
  // open files for BWTs (only reading)
  open_bw_files(g);
  // create 0 initialized bitfile (reading and writing using bitfile_* functions)
  bitfile b;
  bitfile_create(&b,g->mergeLen,g->outPath,g->dbOrder);
  // allocate Z (merge) Znew (reading only Z, writing only newZ) 
  alloc_merge_arrays(g);
  
  // allocate other useful arrays
  g->inCnt = malloc(g->numBwt*sizeof(customInt)); 
  g->firstColumn = malloc(g->sizeOfAlpha*sizeof(customInt));
  g->F = malloc(g->sizeOfAlpha*sizeof(customInt));
  if(!g->inCnt || !g->firstColumn || !g->F) die(__func__);
  // init the above arrays
  init_arrays128ext(g);

  // init liquid block (containing list of allocated mem)
  liquidBlock *liquid = liquid_new(g); 
  // init list (on disk) of irrelevant blocks, initially empty 
  solidBlockFile *ibList = ibHead_new(g);
  uint64_t maxSolid = 0;   // maximum space used by a pair of solid block files 

  // main loop
  uint32_t prefixLength = 1;                      
  int lcpSize = POS_SIZE + BSIZE;   // number of bytes for each pos,lcp pair, see writeLcp()
  bool merge_completed;
  do {
    prefixLength+= 1; 
    if(g->lcpCompute && prefixLength-2>MAX_LCP_SIZE) {fprintf(stderr,"LCP too large: %u\n", prefixLength-2);die(__func__);}
    ibList->fout = gap_tmpfile(g->outPath);
    merge_completed=addCharToPrefix128ext(ibList,liquid,prefixLength,&b,g);
    if (g->verbose>1 && lastRound) {
      // at this iteration we are discovering suffixes with LCP=prefixLength-1
      printf("Lcp: %u. Memory peak/current: %.2lf/%.2lf bytes/symbol. ibList: %ju\n", 
           prefixLength-1, (double)malloc_count_peak()/g->mergeLen,
           (double)malloc_count_current()/g->mergeLen, (uintmax_t) ftello(ibList->fout));
      // also EOF are written to unsortedLcp file so percentages are not accurate     
      if(g->unsortedLcp) printf("   unsorted lcp values: %ju (%.2lf%%)\n",  
      (uintmax_t) ftello(g->unsortedLcp)/lcpSize, (double) 100*ftello(g->unsortedLcp)/(lcpSize*g->mergeLen));
    }
    // compute current space usage for the gap files 
    uint64_t totSolid = ftell(ibList->fout);
    if(ibList->fin!=NULL) {totSolid += ftell(ibList->fin); fclose(ibList->fin);}
    if(totSolid>maxSolid) maxSolid = totSolid;
    // update solid block files:
    rewind(ibList->fout);
    ibList->fin = ibList->fout;
  } while(!merge_completed && (prefixLength!=g->dbOrder));  // end main loop
  if(ibList->fin!=NULL) fclose(ibList->fin);

  if (g->verbose>0) {
    if(lastRound) {
      printf("Merge128ext completed (%d bwts). Mem: %zu peak, %zu current, %.2lf/%.2lf bytes/symbol\n", g->numBwt, malloc_count_peak(),
           malloc_count_current(), (double)malloc_count_peak()/g->mergeLen,
           (double)malloc_count_current()/g->mergeLen);
      printf("Peak solid block disk space: %lu, %.2lf bytes/symbol\n", maxSolid, (double)maxSolid/g->mergeLen);
    }
    else if(g->verbose>1)
      printf("Merge128ext completed (%d bwts). Mem: %zu peak, %zu current\n", g->numBwt, malloc_count_peak(),
           malloc_count_current());
           
  }
  liquid_free(liquid);
  ibHead_free(ibList);
  bitfile_destroy(&b);
  close_bw_files(g);
  
  // for dbGraph info we need to extract the hi bit from the merge arrray
  if(g->dbOrder>1) 
    extract_bitfile(g->merge_fname, g->mergeLen, g->outPath, g->dbOrder);
  
  // computation complete, do the merging. 
  // The following call writes the merged BWT back to file g->bwfname
  mergeBWT128ext(g,lastRound);
  if(g->lcpCompute) {
    assert(lastRound);
    close_unsortedLCP_files(g);
    if(g->verbose>0) printf("Remind to run lcpmerge to obtain the final LCP array\n"); 
  }
  // if dbOrder>1 rename BWT file, since it is a partially ordered BWT
  if(g->dbOrder>1) {
    char tmp1[Filename_size];
    snprintf(tmp1,Filename_size,"%s.%d.%s",g->outPath,g->dbOrder, BWT_EXT);
    if(rename(g->bwfname,tmp1)!=0)
      die("Cannot rename BWT file");
  }

  free(g->F); // last five arrays deallocated
  free(g->firstColumn); 
  free(g->inCnt);
  free_merge_arrays(g);// merge and newMerge
}
