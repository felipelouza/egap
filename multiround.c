#include "util.h"
#include "io.h"
#include "alphabet.h"
#include "mergegap.h"
#include "mergehm.h"
#include "threads.h"

#define min(a,b) ((a)<(b) ? (a) : (b))


// ------ merge of multiple BWTs possibly in parallel
// the number of BWTs to be merged is in input->numBWT
// in each round at most group_size BWTs can be merged this number
// is limited by the size of the elements used to store IDs for BWTs (ZSIZE)
// The merging is done in rounds; at each round the total number of active BWTs decrease
// by a factor group_size. Each round must be completed before we can start the next one
// A round may consist of several parallel merges, and in this cases they can be computed
// in parallel using multiple threads. 
// A producer/consumers system and the relative threads are setup for each
// invocation of multiround, but all rounds share the the same system and threads
// Since consumer threads do not termiate at the end of a round, the main thread 
// uses a condition variable to wait for all thread computations to be completed 
void multiround(bool hm, int group_size,char *path, g_data *input, int num_threads)
{
  customInt offset, tot_symb=0;
  // threads related
  pc_system merge;
  pthread_t t[num_threads];
  g_data buffer[Threads_buf_size];
  int e;
  
  // init threads
  if(num_threads>0) {
    pc_system_init(&merge,Threads_buf_size);
    merge.buffer = (void *) buffer;
    merge.hm = hm;
    for(int i=0;i<num_threads;i++) {
      e = pthread_create(&t[i], NULL, merger, &merge);
      if(e) die("multiround create");
    }
    if(input->verbose>0) printf("%d merger threads created\n",num_threads);
  }
  
  // multiround merging  
  for(int r=0; ;r++) {     // do as many rounds as necessary
    
    if(input->numBwt<= 2*group_size-1 ) 
       break;  // we have reached the last round (or almost)
       
    // main loop to be executed by parallel threads
    // offset is within the number of bwt
    if(num_threads>0) {
      e = pthread_mutex_lock(&merge.remutex);
      if(e) die("remaining lock 1");
      merge.remaining = (input->numBwt+group_size-1)/group_size; // parallel merges in this round
      e= pthread_mutex_unlock(&merge.remutex);
      if(e) die("remaining unlock 1");
    }
    tot_symb = 0; 
    for(offset=0; offset < input->numBwt; offset += group_size) {
      g_data g=*input; // copy current status
      if(g.verbose>2) printf("Round %d, NumBwt offset "CUSTOM_FORMAT"\n",r,offset);
      g.bws = input->bws+offset;
      g.bwtLen = input->bwtLen+offset;
      if(input->smallAlpha) g.bwtOcc = input->bwtOcc+offset;
      g.numBwt = min(group_size, input->numBwt - offset);
      g.symb_offset = tot_symb;
      g.mergeLen=0;
      for(int i=0;i<g.numBwt;i++) g.mergeLen += g.bwtLen[i];
      if(input->lcpCompute) {
        assert(!input->bwtOnly && !input->lcpMerge);
        g.bwtOnly = true;      // since this is not the last round we can only
        g.lcpCompute = false;  // compute the BWT and ignore LCP  
      }
      // prepare access to lcp values
      if(g.lcpMerge) g.lcps = input->lcps+offset;
      // execute merge, possibly using threads 
      if(num_threads>0) {
        e = sem_wait(&merge.free_slots); if(e) die("multiround producer wait");
        check_g_data(&g);
        buffer[merge.pindex++ % merge.buf_size]=g; // copy g_data to buffer
        e = sem_post(&merge.ready);  if(e) die("multiround producer post");
      }
      else { // no helper threads
        if (hm) holtMcMillan(&g, false);
        else gap(&g, false);
        // update input via g
        for(int i=1;i<g.numBwt;i++) {
          g.bwtLen[0] += g.bwtLen[i];
          if(g.smallAlpha)
            for(int j=0;j<g.sizeOfAlpha;j++)
              g.bwtOcc[0][j] += g.bwtOcc[i][j];
        }
        assert(g.bwtLen[0]==g.mergeLen);
      }
      // update tot_symb seen so far
      tot_symb += g.mergeLen;
    } // end of inner parallelizable loop
    assert(input->mergeLen==tot_symb);
    // --- wait till all merges for this round have been completed 
    if(num_threads>0) {
      e = pthread_mutex_lock(&merge.remutex); if(e) die("remaining lock 2");
      while(merge.remaining>0) {
        e = pthread_cond_wait(&merge.recond,&merge.remutex);
        if(e) die("remaining wait");
      }
      e = pthread_mutex_unlock(&merge.remutex); if(e) die("remaining unlock 2");
      assert(merge.pindex==merge.cindex);
      if(input->verbose>1) puts("Threads syncronized, round completed");
    }
    // merged done for all segments. merged bwt/lcp are in g.bws[0][] and g.lcps[0] (if g.lcpMerge==true)

    // update input so that it consists of fewer, larger segments. move Len and Occ real entries to the left
    for(offset=0; offset < input->numBwt; offset += group_size) {
      input->bwtLen[offset/group_size] = input->bwtLen[offset];
      input->bws[offset/group_size] = input->bws[offset];
      if(input->lcpMerge) input->lcps[offset/group_size] = input->lcps[offset];;
      // replace with a swap ??? 
      if(input->smallAlpha)
        for(int j=0;j<input->sizeOfAlpha;j++)
          input->bwtOcc[offset/group_size][j] = input->bwtOcc[offset][j];
    }
    input->numBwt = (input->numBwt+group_size-1)/group_size;
    check_g_data(input);
    if(input->verbose>0)  printf("Round %d complete: %d merged bwt/lcp saved\n",r,input->numBwt);
  } // end of all rounds except last
  if(num_threads>0) {  // last round does not use threads, kill consumers
    g_data g; // dummy task
    g.numBwt = 0;
    g.verbose = input->verbose;
    for(int i=0;i<num_threads;i++) {
        e = sem_wait(&merge.free_slots); if(e) die("cleaning producer wait");
        buffer[merge.pindex++ % merge.buf_size]=g; // copy g_data to buffer
        e = sem_post(&merge.ready);  if(e) die("cleaning producer post");
    }
    for(int i=0;i<num_threads;i++) {
      e = pthread_join(t[i], NULL);
      if(e) die("multiround join");
    }
    pc_system_destroy(&merge);
    if(input->verbose>1) printf("%d merge threads destroyed",num_threads);
  }
  // here maybe there is an additional round to go from n<2g, to g;
  // there is some code duplication, but it is an important special case 
  if(input->numBwt > group_size) { 
    offset = group_size-1; // first BWT to consider 
    tot_symb = 0;          // total number of symbols to skip
    for(int i=0;i<offset;i++)
      tot_symb += input->bwtLen[i];
    // prepare merge task  
    g_data g=*input; // copy current status
    if(g.verbose>2) printf("Special round, NumBwt offset "CUSTOM_FORMAT"\n",offset);
    g.bws = input->bws+offset;
    g.bwtLen = input->bwtLen+offset;
    if(input->smallAlpha) g.bwtOcc = input->bwtOcc+offset;
    g.numBwt = input->numBwt - offset;
    assert(g.numBwt <= group_size);
    g.symb_offset = tot_symb;
    g.mergeLen=0;
    for(int i=0;i<g.numBwt;i++) g.mergeLen += g.bwtLen[i];
    if(input->lcpCompute) {
      assert(!input->bwtOnly && !input->lcpMerge);
      g.bwtOnly = true;      // since this is not the last round we can only
      g.lcpCompute = false;  // compute the BWT and ignore LCP  
    }
    if(g.lcpMerge) g.lcps = input->lcps+offset;
    // execute merge 
    if (hm) holtMcMillan(&g, false);
    else gap(&g, false);
    // update input via g
    for(int i=1;i<g.numBwt;i++) {
      g.bwtLen[0] += g.bwtLen[i];
      if(g.smallAlpha)
        for(int j=0;j<g.sizeOfAlpha;j++)
          g.bwtOcc[0][j] += g.bwtOcc[i][j];
    }
    assert(g.bwtLen[0]==g.mergeLen);
    // update input so that it consists of fewer, larger segments. move Len and Occ real entries to the left
    input->numBwt = (input->numBwt-g.numBwt+1);
    assert(input->numBwt==group_size);
    check_g_data(input);
  }
  assert(input->numBwt <= group_size);
  // execute last round (only point where gap/hm with lastRound==true is called) 
  if (hm) holtMcMillan(input, true);
  else gap(input, true);
}
