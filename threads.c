#include "malloc_count/malloc_count.h"
#include "threads.h"
#include "util.h"
#include "mergegap.h"
#include "mergehm.h"


#define min(a,b) ((a)<(b) ? (a) : (b))


// init a produce/consumers system
void pc_system_init(pc_system *pc, int buf_size)
{
  pc->buf_size = buf_size; 
  pc->pindex=0;
  pc->cindex=0;
  int e = sem_init(&(pc->free_slots),0,buf_size);
  if(e) die(__func__);
  e = sem_init(&(pc->ready),0,0);
  if(e) die(__func__);
  e = pthread_mutex_init(&(pc->cmutex),NULL);
  if(e) die(__func__);
  // re condition variable and related mutex
  e = pthread_mutex_init(&(pc->remutex),NULL);
  if(e) die(__func__);
  e = pthread_cond_init(&(pc->recond),NULL);
  if(e) die(__func__);
}

// destroy a producer/consumers system
void pc_system_destroy(pc_system *pc)
{
  assert(pc->pindex==pc->cindex);  
  int e = pthread_cond_destroy(&(pc->recond));
  if(e) die(__func__);
  e = pthread_mutex_destroy(&(pc->remutex));
  if(e) die(__func__);
  e = pthread_mutex_destroy(&(pc->cmutex));
  if(e) die(__func__);
  e = sem_destroy(&(pc->ready));
  if(e) die(__func__);
  e = sem_destroy(&(pc->free_slots));
  if(e) die(__func__);
}

// execute a single call to HM or Gap
void *merger(void *v)
{
  pc_system *pc = (pc_system *) v;
  g_data g, *buffer= (g_data *)pc->buffer;
  
  int tot=0;
  while(1) {
    int e=sem_wait(&pc->ready); // wait there is something to do
    if(e) die("consumer wait");
    e = pthread_mutex_lock(&pc->cmutex); // get esclusive access to queue
    if(e) die("consumer lock");
    g = buffer[(pc->cindex)++ % pc->buf_size];
    e = pthread_mutex_unlock(&pc->cmutex); // fine accesso a zona comune
    if(e) die("consumer unlock");
    e = sem_post(&pc->free_slots);
    if(e) die("consumer post");
    if(g.numBwt==0) break;
    else {
      if(g.verbose>2) printf("Working on range ["CUSTOM_FORMAT","CUSTOM_FORMAT")\n", g.symb_offset,g.symb_offset+g.mergeLen-1);
      if(pc->hm)  holtMcMillan(&g, false);
      else gap(&g, false);
      // merge statistics 
      for(int i=1;i<g.numBwt;i++) {
        g.bwtLen[0] += g.bwtLen[i];
        if(g.smallAlpha)
          for(int j=0;j<g.sizeOfAlpha;j++)
            g.bwtOcc[0][j] += g.bwtOcc[i][j];
      }
      assert(g.bwtLen[0]==g.mergeLen);
      tot++;
      // update # remaining segments when we reach 0 round is completed
      e = pthread_mutex_lock(&pc->remutex); 
      if(e) die("remaining locking in thread");
      if(--pc->remaining==0) {
        e = pthread_cond_signal(&pc->recond);
        if(e) die("remaining signaling in thread");
      }
      e = pthread_mutex_unlock(&pc->remutex); 
      if(e) die("remaining unlocking in thread");
    }    
  }
  if(g.verbose>1)  printf("Thread terminated. %d merges processed\n",tot);
  pthread_exit(NULL); 
}
