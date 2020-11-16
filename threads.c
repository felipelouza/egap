#if MALLOC_COUNT_FLAG
  #include "malloc_count/malloc_count.h"
#endif
#include "threads.h"
#include "util.h"
#include "mergegap.h"
#include "mergehm.h"


#define min(a,b) ((a)<(b) ? (a) : (b))

// force named semaphores for MacOS
#if defined(__APPLE__) && defined(__MACH__) 
#define USE_NAMED_SEMAPHORES
#endif

#ifdef USE_NAMED_SEMAPHORES /* OS X */
  #pragma message "Compiling using named semaphores for MacOS systems"
#endif

// ---- semaphores

// creation and destruction of semaphores using named for MacOS and unnamed for Linux
// for named semaphores 
//   if s==NULL a new semaphore is created using the static num to ensure distinct names 
//   if s!=NULL the last created semaphore is destroyed (unlinked)
// for unnamed semaphores
//   if s==NULL a new sem is allocated initialized and a pointer to it returned  
//   if s!=NULL the passed semaphore is destroyed 
// also check the error message from sem_unlink or sem_destroy functions 
// return a pointer to the newly created semaphore or NULL if ti was s!=NULL
static int Thread_error_wait=3;
static sem_t *xsem_create_destroy(sem_t *s, unsigned int value, int line, const char *file)
{
#ifdef USE_NAMED_SEMAPHORES
  // -----------------------------------------------
  // MacOS: close/open a named semaphore
  // static variables containing pid and # created semaphores 
  char tmp[NAME_MAX+1];
  static int num=0;
  static intmax_t pid;
  // save pid once
  if(num==0) pid = (intmax_t) getpid();

  if(s!=NULL) {
    // check there is a named semaphore still around 
    if(num==0) {
      perror("Error deleting non-existing named semaphore");
      fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
      sleep(Thread_error_wait);  // give some extra time to other threads 
      exit(1);
    }
    // delete last created named semaphore
    int e = snprintf(tmp,NAME_MAX,"%jd.%d",pid,--num);
    if(e<0 || e>NAME_MAX) {
      perror("Error creating semaphore name");
      fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
      sleep(Thread_error_wait);  // give some extra time to other threads 
      exit(1);
    }
    e = sem_unlink(tmp);
    if(e!=0) {
      perror("Error deleting named semaphore");
      fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
      sleep(Thread_error_wait);  // give some extra time to other threads 
      exit(1);
    }
    return NULL;
  }
  // create a new named semphore; create unique name 
  int e = snprintf(tmp,NAME_MAX,"%jd.%d",pid,num++);
  if(e<0 || e>NAME_MAX) {
    perror("Error creating semaphore name");
    fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
    sleep(Thread_error_wait);  // give some extra time to other threads 
    exit(1);
  }
  // open and init semaphore 
  s = sem_open(tmp,O_CREAT| O_EXCL, S_IRUSR | S_IWUSR ,value);
  if(s==SEM_FAILED) {
    perror("Error opening named semaphore"); 
    fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
    sleep(Thread_error_wait);  // give some extra time to other threads 
    exit(1);
  }
  return s;
#else
  // ---------------------------------------------- 
  // linux destroy/init an unnamed semaphore
  if(s!=NULL) {
    if(sem_destroy(s) !=0) {
      perror("Error destroying unnamed semaphore");
      fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
      sleep(Thread_error_wait);  // give some extra time to other threads 
      exit(1);
    }
    free(s);
    return NULL;
  }
  // allocate init/ and return an unnamed sem_t
  s = malloc(sizeof(sem_t));
  if(s==NULL) {
    perror("malloc error");
    fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
    sleep(Thread_error_wait);  // give some extra time to other threads 
    exit(1);
  }
  // init with value and no sharing 
  if(sem_init(s,0,value)!=0) {
    perror("sem_init error");
    fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
    sleep(Thread_error_wait);  // give some extra time to other threads 
    exit(1);
  }
  return s;
#endif
}



// init a producer/consumers system
void pc_system_init(pc_system *pc, int buf_size)
{
  pc->buf_size = buf_size; 
  pc->pindex=0;
  pc->cindex=0;
  pc->free_slots = xsem_create_destroy(NULL,buf_size,__LINE__,__FILE__);
  pc->ready = xsem_create_destroy(NULL,0,__LINE__,__FILE__);
  int e = pthread_mutex_init(&(pc->cmutex),NULL);
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
  xsem_create_destroy(pc->ready,0,__LINE__,__FILE__);
  xsem_create_destroy(pc->free_slots,0,__LINE__,__FILE__);
}

// execute a single call to HM or Gap
void *merger(void *v)
{
  pc_system *pc = (pc_system *) v;
  g_data g, *buffer= (g_data *)pc->buffer;
  
  int tot=0;
  while(1) {
    int e=sem_wait(pc->ready); // wait there is something to do
    if(e) die("consumer wait");
    e = pthread_mutex_lock(&pc->cmutex); // get esclusive access to queue
    if(e) die("consumer lock");
    g = buffer[(pc->cindex)++ % pc->buf_size];
    e = pthread_mutex_unlock(&pc->cmutex); // fine accesso a zona comune
    if(e) die("consumer unlock");
    e = sem_post(pc->free_slots);
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
