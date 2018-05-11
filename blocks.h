/* ******************************************************
 * Space economical implementation of irrelevant blocks
 * 
 * blocks can be of three kinds: proto, liquid or solid
 * 
 * solid blocks are those containing a range that should be skippend in 
 * succesive iterations and for this reason they store the number
 * of occ per symbol and per documents

 * Solid blocks have an occ array of smallSolidInt until their dimension is <= SMALLSOLID_LIMIT
 * When blocks become bigger, occ arrays are inflated to an array of customInt.
 * 
 * At any given moment there is at most one liquid block: it contains 
 * adjacent monochrome blocks. It becomes solid if its size reaches
 * solid_limit. If it touches a solid block it is fused to it.
 * 
 * The proto block consists of the last seen characters, it can become part
 * of a liquid block if it turns out to be monochrome. 
 *
 */

// irrelevant solid block to be skipped at each occurrence
typedef struct block {
  customInt beginsAt;
  customInt endsAt;
  union {
    customInt *occ;    // [0..numberOfBWTs-1] [numberOfBWTs .. ]
    smallSolidInt *smallOcc;
  };
  struct block *nextBlock;
} solidBlock;


// a liquid Block is a merge of monochrome blocks that is made solid if 
// it reaches size solid_limit. If it touches a solid block  
// it is merged to it 
typedef struct liqblock {
  customInt beginsAt;
  customInt endsAt;   // only used in assertions 
  int numBwt;
  int sizeOfAlpha;
  int occ_size;
  customInt *occ;    // [0..numberOfBWTs-1] [numberOfBWTs .. ]   
  int solid_limit;   // must be reached to create a solid block    
  bool empty;
} liquidBlock;


// protoBlock: will become a block if turns out to be monochrome
typedef struct {
  customInt beginsAt;
  customInt endsAt;
  symbol *start;    // occ of alpha symbol (this is a problem in external memory)
  bool mono;
  int color;        // color of the block (to recognize monotone blocks)
  int lastChar;     // last char of block, useful only for singleton blocks 
  int lastColor;    // last color of block, useful only for singleton blocks 
} protoBlock;


typedef struct {
  FILE *fin;
  FILE *fout;
  int occ_size;
  solidBlock *solidList;
  customInt *occList;
  smallSolidInt *smallOccList;
} solidBlockFile;


static void writeBlock(solidBlock *s, solidBlockFile *sf);

// --- debug functions

// scan the list of irrelevant blocks and make sure there are no mergeable
// blocks. If verbose report the number of blocks by type 
// NOTE: this is essentially a development/debug function
// to be changed for use with solidBlockFile
customInt checkBlocks(solidBlock *b, customInt pfx, g_data *g, bool lastRound)
{
#if 0 
  customInt tot=0, mergeable=0, small_solid=0;

  while(b!=NULL) {
    tot++;
    if(b->nextBlock!=NULL && b->endsAt == b->nextBlock->beginsAt) 
      mergeable++; 
    if(b->endsAt-b->beginsAt<= SMALLSOLID_LIMIT) small_solid++;
    b = b->nextBlock;  // check what comes after
  }
  if(g->verbose>2 || (g->verbose>1 && lastRound)) 
    printf("Lcp: "CUSTOM_FORMAT"  Blocks: "CUSTOM_FORMAT" total, "CUSTOM_FORMAT" small_solid\n", pfx-1, tot,small_solid);
  if(mergeable!=0)
    printf("Lcp: "CUSTOM_FORMAT"  Fatal error: "CUSTOM_FORMAT" mergeable blocks found\n",pfx-1,mergeable);
  return mergeable;
#else
  (void) b, (void) pfx, (void) g, (void) lastRound; 
  return 0;   
#endif
}

// debug only
void showIrrelevants(solidBlock *b)
{
  #if 0
  printf("Irrelevants: ");
  while(b!=NULL) {
    assert(b->occ!=NULL);
    if (b->endsAt - b->beginsAt <= SMALLSOLID_LIMIT) 
      printf("("CUSTOM_FORMAT" "CUSTOM_FORMAT") ",b->beginsAt,b->endsAt);
    else 
      printf("["CUSTOM_FORMAT" "CUSTOM_FORMAT"] ",b->beginsAt,b->endsAt);
    b = b->nextBlock;
  }
  puts("END");
#else
  (void) b; 
#endif  
}

// debug only: show occ's
void showBlock(solidBlock *b, g_data *g, int num)
{
  printf("Beg "CUSTOM_FORMAT" End "CUSTOM_FORMAT" ",b->beginsAt, b->endsAt);
  if(b->endsAt - b->beginsAt > SMALLSOLID_LIMIT) {
    printf("occ: ");
    for(int i=0;i<num;i++) printf(CUSTOM_FORMAT" ",b->occ[i]);
  } else {
    printf("smallOcc: ");
    for(int i=0;i<num;i++) printf(SMALLSOLID_FORMAT" ",b->smallOcc[i]);
  }
  puts("SolEND");
}


void showLiquid(liquidBlock *b)
{
  printf("Beg "CUSTOM_FORMAT" End "CUSTOM_FORMAT" ",b->beginsAt, b->endsAt);
  customInt totA=0, totB=0;
  for(int i=0;i<b->numBwt;i++) totB += b->occ[i];
  for(int i=0;i<b->sizeOfAlpha;i++) totA += b->occ[b->numBwt+i];
  printf("Chars: "CUSTOM_FORMAT" Bwts: "CUSTOM_FORMAT" Empty: %d ",totA,totB,b->empty);  
  puts("LiqEND");
}



// --- alloc/dealloc functions

// create an empty liquid block, called once for iteration 
liquidBlock *liquid_new(g_data *g){
  liquidBlock *newB = malloc(sizeof(*newB));
  if(newB==NULL) die("Out of mem in liquid_new()");
  newB->beginsAt = 0; 
  newB->endsAt = 0; 
  newB->empty = true; 
  newB->occ_size = g->numBwt+g->sizeOfAlpha;
  newB->numBwt = g->numBwt;
  newB->sizeOfAlpha = g->sizeOfAlpha;
  newB->solid_limit = max(g->solid_limit,newB->occ_size);
  // solid limit cannot be smaller than mergeLen  
  newB->solid_limit = min(newB->solid_limit,g->mergeLen);  
  newB->occ = calloc(newB->occ_size,sizeof(customInt));
  if(newB->occ==NULL) die("Out of mem in liquid_new()");
  return newB;
}

// free liquid block, called once for iteration 
void liquid_free(liquidBlock *b)
{
  assert(b!=NULL && b->occ!=NULL);
  free(b->occ);
  free(b);
}

// init ibHead containing file pointers for solid blocks and local heap
solidBlockFile *ibHead_new(g_data *g) {
  solidBlockFile *ibList = malloc(sizeof(*ibList));
  if(ibList==NULL) die("Out of mem in inHead_new");
  ibList->fin = ibList->fout = NULL;
  ibList->occ_size = g->sizeOfAlpha + g->numBwt;
  ibList->solidList = NULL;
  ibList->occList = NULL;
  ibList->smallOccList = NULL;
  return ibList;
}

// free the local heap: called once for iteration 
void ibHead_free(solidBlockFile *b)
{
  // free list of solid blocks
  assert(b!=NULL);
  {
    solidBlock *next=b->solidList; int tot=0;
    while(next!=NULL) {
      solidBlock *tmp = next->nextBlock;
      free(next); next = tmp; tot++;
    }
    assert(tot<=3); // blocks are on disk: at most three in mem simultaneously 
  }
  // free list of occList
  {customInt *next=b->occList;
  while(next!=NULL) {
    customInt **n = (customInt **) next; // pretend next is an array of pointers customInt *
    customInt *tmp = n[0];
    free(next); next = tmp;
  }}
  // free list of smallOccList
  {smallSolidInt *next= b->smallOccList;
  while(next!=NULL) {
    smallSolidInt **n = (smallSolidInt **) next; // pretend next is an array of pointers customInt *
    smallSolidInt *tmp = n[0];
    free(next); next = tmp;
  }}
  free(b);
}  
  


// ----- get solid block from s->solidList or malloc
solidBlock *get_block(solidBlockFile *s)
{
  solidBlock *newB =  s->solidList;
  if(newB!=NULL) {
    s->solidList = newB->nextBlock;
  }
  else {
    newB = malloc(sizeof(*newB));
    if(newB==NULL) die("Out of mem in get_block()");
  }
  return newB;
}

// get an occ list: if possible from s->occList otherwise via malloc 
customInt *get_occ(solidBlockFile *s) {
  customInt *occ = s->occList;  // try occlist
  if(occ == NULL) {
    occ = malloc(s->occ_size* sizeof(customInt));
    if(occ==NULL) die("Out of mem in get_occ()");
  }
  else { // update list occList
    customInt **n = (customInt **) occ;
    s->occList = n[0];
  }
  return occ;  
}

// get a  smallOcc list: if possible from s->occList otherwise via malloc 
smallSolidInt *get_smallocc(solidBlockFile *s) {
  smallSolidInt *occ = s->smallOccList;  // try smallOcclist
  if(occ == NULL) {
    occ = malloc(s->occ_size* sizeof(smallSolidInt));
    if(occ==NULL) die("Out of mem in get_smallocc()");
  }
  else { // update list occList
    smallSolidInt **n = (smallSolidInt **) occ;
    s->smallOccList = n[0];
  }
  return occ;  
}


// add to s->occList the array occ
void free_occ(customInt *occ, solidBlockFile *s) {
    customInt **n = (customInt **) occ; // pretend occ is an array of pointers
    n[0] = s->occList;
    s->occList = occ;
}

// add to s->occList the array occ
void free_smallocc(smallSolidInt *occ, solidBlockFile *s) {
    smallSolidInt **n = (smallSolidInt **) occ; // pretend occ is an array of pointers
    n[0] = s->smallOccList;
    s->smallOccList = occ;
}



// free the memory used by  solid block 
// saving its space to x->solidList;
// called after a three way merge merge_sls when a block is destroyed
// or after a block has been written to disk by writeBlock 
void block_free(solidBlock *b, solidBlockFile *x)
{
  assert(b!=NULL && b->occ != NULL);
  // free occ's
  if(b->endsAt -b->beginsAt > SMALLSOLID_LIMIT)
    free_occ(b->occ,x);
  else // smallOcc  
    free_smallocc(b->smallOcc,x);
  // free block itself  
  b->nextBlock = x->solidList;
  x->solidList = b;
}

/* * substitutes the smallOcc array with a customInt array
 * \param b  pointer to block to be processed
 * \param s: size occ array 
 * */
static void block_inflate(solidBlock *b, solidBlockFile *s) {
  assert(b->endsAt - b->beginsAt <= SMALLSOLID_LIMIT);
  customInt *bigOcc = get_occ(s);
  for (int i = 0; i < s->occ_size; ++i) {
    bigOcc[i] = (customInt) b->smallOcc[i];
  }
  free_smallocc(b->smallOcc,s);
  b->occ = bigOcc;
}


/* * creates a new block to store data from the liquid block s
 *  \param s pointer to liquid block
 * 
 * Note: this is the only function where solid blocks are created
 * */
static solidBlock *liquidToSolid(liquidBlock *s, solidBlockFile *sf) {
  assert(s->endsAt - s->beginsAt >= s->solid_limit);
  
  // ----- get solid block from s->solidList or malloc
  solidBlock *newB = get_block(sf);
  // --- init block size
  newB->beginsAt = s->beginsAt;
  newB->endsAt = s->endsAt;
  // --- get and init occ/smallOcc array
  if (newB->endsAt - newB->beginsAt <= SMALLSOLID_LIMIT) {
    // small solid block  try get mem from s->occList
    newB->smallOcc = get_smallocc(sf);
    for (int i = 0; i < s->occ_size; ++i)  // copy occ from s
      newB->smallOcc[i] = (smallSolidInt) s->occ[i];
  } 
  else { // large solid block  try get mem from s->occList
    newB->occ = get_occ(sf);
    for (int i = 0; i < s->occ_size; ++i) // copy occ from s
      newB->occ[i] = s->occ[i];
  }
  // clear liquid block 
  s->empty=true;
  memset(s->occ, 0, (s->occ_size) * sizeof(customInt));
  return newB;
}



// --------- single block functions


/* * skip an irrelevant block
 * \param this  pointer to block to be skipped
 * \param g:    all general data
 * */
static void skip(solidBlock *this, g_data *g)
{
  if(g->extMem && g->fmergeColor!=NULL) {// this is because merge8 supports extMem but not extMem colors 
    fseek(g->fmergeColor,(this->endsAt-this->beginsAt)*sizeof(palette),SEEK_CUR);
    assert(ftell(g->fmergeColor)==this->endsAt*sizeof(palette));
  }
  if(this->endsAt - this->beginsAt <= SMALLSOLID_LIMIT){
    // skip in bwt array 
    for (int col = 0; col < g->numBwt; col++) 
      if(this->smallOcc[col]>0) {
        g->inCnt[col] += this->smallOcc[col];    //skipping this
        if(g->extMem) fseek(g->bwf[col], this->smallOcc[col]*sizeof(symbol), SEEK_CUR);
      }
    // skip in new_merge array  
    for (int c = 0; c < g->sizeOfAlpha; ++c)
      if(this->smallOcc[g->numBwt+c]>0) {
        g->F[c]+= this->smallOcc[g->numBwt+c];
        if(g->extMem && c>0 && g->fnewMergeColor!=NULL) {
          cwriter_skip(&g->fnewMergeColor[c],this->smallOcc[g->numBwt+c]);
          assert(cwriter_tell(&g->fnewMergeColor[c])==g->F[c]*sizeof(palette));
        }
      }
  } 
  else { 
    for (int col = 0; col < g->numBwt; col++) 
      if(this->occ[col]>0) {
        g->inCnt[col] += this->occ[col];        //skipping this
        if(g->extMem) fseek(g->bwf[col], this->occ[col]*sizeof(symbol), SEEK_CUR);
      }
    for (int c = 0; c < g->sizeOfAlpha; ++c) 
      if(this->occ[g->numBwt+c]>0) {
        g->F[c] += this->occ[g->numBwt+c];
        if(g->extMem && c>0 && g->fnewMergeColor!=NULL) {
          cwriter_skip(&g->fnewMergeColor[c],this->occ[g->numBwt+c]);
          assert(cwriter_tell(&g->fnewMergeColor[c])==g->F[c]*sizeof(palette));
        }
      }
  }
}

static void skip128ext(solidBlock *this, bitfile *b, g_data *g)
{
  // we have already read the color at position this->beginsAt (we needed the bit value) 
  fseek(g->fmergeColor,((this->endsAt-this->beginsAt)-1)*sizeof(palette),SEEK_CUR);
  assert(ftell(g->fmergeColor)==this->endsAt*sizeof(palette));
  
  if(this->endsAt - this->beginsAt <= SMALLSOLID_LIMIT){
    // skip in bwt array 
    for (int col = 0; col < g->numBwt; col++) 
      if(this->smallOcc[col]>0) {
        g->inCnt[col] += this->smallOcc[col];    //skipping this
        fseek(g->bwf[col], this->smallOcc[col]*sizeof(symbol), SEEK_CUR);
      }
    // skip in new_merge array  
    for (int c = 0; c < g->sizeOfAlpha; ++c)
      if(this->smallOcc[g->numBwt+c]>0) {
        g->F[c]+= this->smallOcc[g->numBwt+c];
        if(c>0) {
          cwriter_skip(&g->fnewMergeColor[c],this->smallOcc[g->numBwt+c]);
          assert(cwriter_tell(&g->fnewMergeColor[c])==g->F[c]*sizeof(palette));
        }
      }
  } 
  else { 
    for (int col = 0; col < g->numBwt; col++) 
      if(this->occ[col]>0) {
        g->inCnt[col] += this->occ[col];        //skipping this
        fseek(g->bwf[col], this->occ[col]*sizeof(symbol), SEEK_CUR);
      }
    for (int c = 0; c < g->sizeOfAlpha; ++c) 
      if(this->occ[g->numBwt+c]>0) {
        g->F[c] += this->occ[g->numBwt+c];
        if(c>0) {
          cwriter_skip(&g->fnewMergeColor[c],this->occ[g->numBwt+c]);
          assert(cwriter_tell(&g->fnewMergeColor[c])==g->F[c]*sizeof(palette));
        }
      }
  }
  // skip in Bit array
  long oldt = bitfile_tell(b);
  bitfile_skip(b,(this->endsAt-this->beginsAt)-1);
  if(bitfile_tell(b)!=this->endsAt) 
    printf("oldt: %ld, tell: %ld, start %ld, end: %ld\n",oldt,bitfile_tell(b),this->beginsAt,this->endsAt);
}



// add proto block to the current liquid block
static void add_proto2liquid(protoBlock *proto,liquidBlock *liquid)
{
  assert(liquid->endsAt==proto->beginsAt);
  assert(proto->endsAt>proto->beginsAt+1);
  liquid->endsAt=proto->endsAt;
  // update color occurrences
  assert(proto->color==proto->lastColor); 
  liquid->occ[proto->color] += (proto->endsAt-proto->beginsAt);
  // update alpha occurences 
  for(int i=0;i<proto->endsAt-proto->beginsAt;i++)
    liquid->occ[liquid->numBwt + (proto->start)[i]] +=1;
  liquid->empty = false;
}


// add singleton proto block to the current liquid block
static void add_singleton2liquid(protoBlock *proto,liquidBlock *liquid)
{
  assert(liquid->endsAt==proto->beginsAt);
  assert(proto->endsAt==proto->beginsAt+1);
  liquid->endsAt=proto->endsAt;
  // update color occurrences 
  liquid->occ[proto->lastColor] += 1;
  // update alpha occurences 
  liquid->occ[liquid->numBwt + proto->lastChar] +=1;
  liquid->empty = false;
}


// ---- merge functions 

// merge liquid block b to solid s
// b is cleared after the merging 
static void merge_liquid(liquidBlock *b, solidBlock *s, solidBlockFile *sf) {
  assert(!b->empty);
  // --- merge occs. 
  if((s->endsAt - s->beginsAt) + (b->endsAt - b->beginsAt) <= SMALLSOLID_LIMIT){
    for (int c = 0; c < b->occ_size; ++c) // merged block is small 
      s->smallOcc[c] += (smallSolidInt) b->occ[c];
  } 
  else {
    if((s->endsAt - s->beginsAt)<= SMALLSOLID_LIMIT)
      block_inflate(s,sf);       // enlarge if required 
    for (int c = 0; c < b->occ_size; ++c) 
      s->occ[c] += b->occ[c];
  }
  // --- update extremes
  if(b->endsAt==s->beginsAt)
    s->beginsAt = b->beginsAt;  // liquid before solid
  else {
    assert(s->endsAt==b->beginsAt); // solid before liquid
    s->endsAt = b->endsAt;
  }
  // --- clear b
  b->empty = true;
  memset(b->occ, 0, (b->occ_size) * sizeof(customInt));
}


// check if liquid block b can be merged to last or made a solid by itself
// here is where solid_limit is used  
// in any case b is emptied
solidBlock *finalize_liquid(solidBlock *last, liquidBlock *b, solidBlock *next, solidBlockFile *sf)
{
  assert(!b->empty);
  if(last!=NULL && last->endsAt==b->beginsAt) {
    // merge b into last and empty it 
    merge_liquid(b,last,sf);  // last and next do not change
  }
  else if(b->endsAt-b->beginsAt>= b->solid_limit) {
    // create new solid and empty b which becomes last block 
    solidBlock *s = liquidToSolid(b,sf);  
    s->nextBlock = next;
    if(last!=NULL) writeBlock(last,sf); // it was: last->nextBlock = s;
    last = s;
  }
  else { // nothing to do, just empty b
    b->empty = true;
    memset(b->occ, 0, (b->occ_size) * sizeof(customInt));
  }
  assert(b->empty);
  return last;
}


// merge solid s to liquid b to solid z
// b is cleared after the merging. z is deallocated 
// NOTE: this is the only point where a solid block is destroyed 
static void merge_sls(solidBlock *s, liquidBlock *b, solidBlock *z, solidBlockFile *sf) {
  // showBlock(s,NULL,0); showBlock(z,NULL,0); showLiquid(b);
  assert(s && z && b && !b->empty);
  assert(s->endsAt==b->beginsAt && b->endsAt==z->beginsAt);
  customInt ssize = s->endsAt - s->beginsAt;
  customInt zsize = z->endsAt - z->beginsAt;
  bool zreduced = false;
  // --- merge occs 
  if(ssize + zsize + (b->endsAt - b->beginsAt) <= SMALLSOLID_LIMIT){
    for (int c = 0; c < b->occ_size; ++c) // all blocks are small 
      s->smallOcc[c] += z->smallOcc[c] + (smallSolidInt) b->occ[c];
  } 
  else if(ssize>SMALLSOLID_LIMIT && zsize>SMALLSOLID_LIMIT) {
      for (int c = 0; c < b->occ_size; ++c) // all blocks are large  
        s->occ[c] += z->occ[c] + b->occ[c];
  }
  else { // mixed small/large blocks
    if(zsize>SMALLSOLID_LIMIT) {// s must be  small, swap with z
      customInt *tmp = z->occ; z->smallOcc = s->smallOcc; s->occ = tmp; zreduced=true;
    }
    else {  // z is small, s can be small or large 
      if(ssize<=SMALLSOLID_LIMIT)  // if s is small enlarge it
        block_inflate(s,sf);
    }
    // now s is large and z is small
    for (int c = 0; c < b->occ_size; ++c) 
        s->occ[c] += z->smallOcc[c] + b->occ[c];
  }
  // ----  update endpoints and nextBlock
  s->endsAt = z->endsAt;
  s->nextBlock = z->nextBlock;
  // --- clear b and free z
  b->empty = true;
  memset(b->occ, 0, (b->occ_size) * sizeof(customInt));
  if(zreduced) z->endsAt = z->beginsAt; // a large z got a smallocc from s, make it appear small 
  block_free(z,sf);
}


// functions handling the two-bit array used in lieu of B
// for bwt-only gap algorithm. 

// alloc and set to zero a two-bit array
uint64_t *tba_alloc(customInt n, bool mmapB)
{
  uint64_t *a;
  customInt j = (n+31)/32;
  if(mmapB) {
    a = mmap(NULL,j*sizeof(*a),PROT_READ|PROT_WRITE,MAP_PRIVATE|MAP_ANONYMOUS,-1,0);
    if(a== MAP_FAILED) die(__func__);
  }
  else {
    a = calloc(j,sizeof(*a));
    if(a==NULL) die(__func__);
  }
  return a;
}

void tba_free(uint64_t *a, customInt n, bool mmap)
{
  if(mmap) {
    customInt j = (n+31)/32;
    int e = munmap(a,j*sizeof(*a));
    if(e!=0) die(__func__);
  }
  else free(a);
}
// test a[i]==0 (not used)
static inline bool tba_is_zero(uint64_t *a,customInt i)
{
  customInt q = i/32;
  customInt r = i%32;
  uint64_t mask = (uint64_t) 3 << (2*r);
  return (a[q] & mask)==0;
}
// test a[i]==11 (not used)
static inline bool tba_is_11(uint64_t *a,customInt i)
{
  customInt q = i/32;
  customInt r = i%32;
  return ((a[q] >> (2*r)) & 3) ==3; // shift and check last 2 bits 
} 
// set a[i] = a[i] | m
static inline void tba_or_m(uint64_t *a, customInt i, uint64_t m)
{
  assert(m==1 || m==2);
  customInt q = i/32;
  customInt r = i%32;
  a[q] = a[q] | (m<< (2*r));
}
static inline void tba_mark_if0(uint64_t *a, customInt i, uint64_t m)
{
  assert(m==1 || m==2);
  customInt q = i/32;
  customInt r = i%32;
  if( ((a[q]>> (2*r)) & 3) == 0)
    a[q] = a[q] | (m<< (2*r));
}
// return true if we are at the begining of a block 
// that is if a[i]==11 || a[i]== ~m   (note: m is 01 or 10)
// if a[i]==11 set recent to false, 
// if a[i]==~m set recent to true and set a[i] to 11 (for next iterations)
// otherwise return false and do not change recent 
static inline bool tba_block_test_set(uint64_t *a,customInt i, int m, bool *recent)
{
  assert(m==1 || m==2);
  customInt q = i/32;
  customInt r = i%32;
  int b2 = (a[q] >> (2*r)) & 3;       // shift and check last 2 bits
  if(b2==3-m) {
    *recent = true;
    a[q] = a[q] | ( (uint64_t) 3 << (2*r)); // make the block not recent 
    return true;
  }
  else if(b2==3) {
    *recent = false;
    return true;
  }
  return false; // do not change recent 
} 
// report value
int tba_get(uint64_t *a,customInt i)
{
  customInt q = i/32;
  customInt r = i%32;
  int b2 = (a[q] >> (2*r)) & 3;       // shift and check last 2 bits
  return b2;
}


// ---- same tba functions for the case the B array stored in the two most signficat bits of the merge array 

// set a[i] = a[i] | m
static inline void tba_or_m8(uint8_t *a, customInt i, uint64_t m)
{
  assert(m==1 || m==2);
  a[i] = a[i] | (m<< 6);
}
static inline void tba_mark8_if0(uint8_t *a, customInt i, uint64_t m)
{
  assert(m==1 || m==2);
  if( (a[i] & 0xC0) == 0)   // if 2 msb are 0
    a[i] = a[i] | (m<< 6);  // write m on them 
}
// return true if we are at the begining of a block 
// that is if a[i]==11 || a[i]== ~m   (note: m is 01 or 10)
// if a[i]==11 set recent to false, 
// if a[i]==~m set recent to true and set a[i] to 11 (for next iterations)
// otherwise return false and do not change recent  
static inline bool tba_block_test_set8(uint8_t *a,customInt i, int m, bool *recent)
{
  assert(m==1 || m==2);
  int b2 = (a[i] >> 6) & 3;       // shift and check 2 msb 
  if(b2==3-m) {
    *recent = true;
    a[i] = a[i] | ( 3 << 6); // make the block not recent 
    return true;
  }
  else if(b2==3) {
    *recent = false;
    return true;
  }
  return false; // do not change recent 
} 


// ---- same tba functions for the case the B array id stored in the two most signficant bits 
// of the merge array (each entry now 16 bits)

// set a[i] = a[i] | m
static inline void tba_or_m16(uint16_t *a, customInt i, uint64_t m)
{
  assert(m==1 || m==2);
  a[i] = a[i] | (m<< 14);
}
// set a[i] = a[i] | m if the two most significant bit are 0
static inline void tba_mark16_if0(uint16_t *a, customInt i, uint64_t m)
{
  assert(m==1 || m==2);
  if( (a[i] & 0xC000) == 0)   // if 2 msb are 0
    a[i] = a[i] | (m<< 14);  // write m on them 
}
// return true if we are at the begining of a block 
// that is if a[i]==11 || a[i]== ~m   (note: m is 01 or 10)
// if a[i]==11 set recent to false, 
// if a[i]==~m set recent to true and set a[i] to 11 (for next iterations)
// otherwise return false and do not change recent  
static inline bool tba_block_test_set16(uint16_t *a,customInt i, int m, bool *recent)
{
  assert(m==1 || m==2);
  int b2 = (a[i] >> 14) & 3;       // shift and check 2 msb 
  if(b2==3-m) {
    *recent = true;
    a[i] = a[i] | ( 3 << 14); // make the block not recent 
    return true;
  }
  else if(b2==3) {
    *recent = false;
    return true;
  }
  return false; // do not change recent 
} 





// ===== handle blocks on file ===== 

// read a solid block from file 
// the liquid block s is used only to access the local heap and for occ_size
static solidBlock *readBlock(solidBlockFile *sf) 
{
  if(sf->fin==NULL) return NULL;
  // read size of block
  customInt buffer[SIZE_OF_ALPHABET+MAX_NUMBER_OF_BWTS];
  int e = fread(buffer,sizeof(customInt),2,sf->fin);
  if(e==0) return NULL; // no more blocks
  if(e!=2) die("tmp file read error in readBlock (1)");
  // ----- get solid block from s->solidList or malloc
  solidBlock *newB = get_block(sf);  
  newB->beginsAt = buffer[0];
  newB->endsAt = buffer[1];
  // printf("Reading [%ld %ld]\n",newB->beginsAt,newB->endsAt); //!!!!!!
  // --- get and init occ/smallOcc array
  if (newB->endsAt - newB->beginsAt <= SMALLSOLID_LIMIT) {
    // small solid block
    newB->smallOcc = get_smallocc(sf);
    e = fread(newB->smallOcc,sizeof(smallSolidInt),sf->occ_size,sf->fin);
    if(e!=sf->occ_size) die("tmp file read error in readBlock (2)");
  }
  else { // large solid block 
    newB->occ = get_occ(sf);
    memset(newB->occ,0,sf->occ_size*sizeof(customInt));// zero all values
    uint8_t *byteBuffer = (uint8_t *) buffer; // transform to an uint8_t buffer
    e = fread(byteBuffer,5,sf->occ_size,sf->fin); // read 5 int for customInt
    if(e!=sf->occ_size) die("tmp file read error in readBlock (3)");
    // fill newB->occ using five bytes for entry
    for(int j=0;j<sf->occ_size;j++) {
      for(int i=0;i<5;i++) 
        newB->occ[j] |= ((customInt) byteBuffer[j*5+i]) << (8*i);  
    }
  }
  newB->nextBlock=NULL;
  return newB;
}    

// write solid block s to file and deallocate it
static void writeBlock(solidBlock *s, solidBlockFile *sf)
{
  assert(s!=NULL);
  assert(sf!=NULL && sf->fout!=NULL);
  customInt buffer[SIZE_OF_ALPHABET+MAX_NUMBER_OF_BWTS];
  // printf("Writing [%ld %ld]\n",s->beginsAt,s->endsAt); //!!!!!!
  buffer[0] = s->beginsAt;
  buffer[1] = s->endsAt;
  int e = fwrite(buffer,sizeof(customInt),2,sf->fout);
  if(e!=2) die("tmp file write in writeBlock (1)");
  if (s->endsAt - s->beginsAt <= SMALLSOLID_LIMIT) {
    // small solid block
    e = fwrite(s->smallOcc,sizeof(smallSolidInt),sf->occ_size,sf->fout);
    if(e!=sf->occ_size) die("tmp file write error in writeBlock (2)");
  }  
  else { // large solid block 
    uint8_t *byteBuffer = (uint8_t *) buffer;
    // fill newB->occ using five bytes for entry
    for(int j=0;j<sf->occ_size;j++) {
      for(int i=0;i<5;i++)
        byteBuffer[j*5+i] = (s->occ[j]>>(8*i)) & 0xFF;  
    }
    e= fwrite(byteBuffer,5,sf->occ_size,sf->fout); // read 5 int for customInt
    if(e!=sf->occ_size) die("tmp file write error in writeBlock (3)");
  }
  block_free(s,sf);
}
