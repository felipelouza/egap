// Sample implementation of the GAP algorithm 
// Copyright 2017 Lavina Egidi, Giovanni Manzini
//
// External memory version:
//  * solid blocks are kept in two external files:
//    one for current blocks, the other for next iteration blocks
//  * no squeezing: 1) it was not very effective in ram, 2) it increases
//    alphabet size and therefore the number of access points to mmaped memory
//    3) squeezing was still not working for multi-bwt inputs       
#include "util.h"
#include "alphabet.h"
#include "gap.h"
#include "malloc_count/malloc_count.h"

static void writeOutputBWT(g_data *g, char *path, bool hm);

// from multiround.c
void multiround(bool hm, int group_size,char *path, g_data *input, int);


void usage(char *name, g_data *g){
  printf("\nUsage: %s [options] PATH\n\n",name);
  puts("Merges BWTs (and optionally LCPs) using the HM or Gap (default) algorithm");
  printf("each LCP value is represented using %zu byte(s).\n\n",sizeof(lcpInt));

  printf("The input bwt's and lcp's must be in PATH.bwt and PATH.%zu.lcp\n", sizeof(lcpInt));
  puts("and the bwt's lengths stored as 8 byte uint's in PATH."LEN_EXT"\n"); 

  puts("Available options:");
  puts("\t-h    this help message");
  puts("\t-o    basename for output files (def. input basename)");
  puts("\t-r    merge lcp values (overwrite input LCPs)");
  puts("\t-l    compute lcp values");
  puts("\t-d D  create document array using D bytes per entry, ext: ."DA_EXT);
  puts("\t-S S  create suffix array using S bytes per entry, ext: ."SA_EXT);
  puts("\t-x    compute lcp without external mergesort");
  puts("\t-a    assume alphabet is small");   
  printf("\t-g G  max # BWTs merged simultaneously (def %llu)\n", MAX_NUMBER_OF_BWTS);   
  puts("\t-A a  preferred gap algorithm to use (see doc or leave it alone)");   
  puts("\t-m    use H&M algorithm");
  printf("\t-s S  minimum solid block size (def %d)\n",g->solid_limit);
  puts("\t-p P  use P parallel threads for merging (def 0)");
  puts("\t-E    run in external memory");
  puts("\t-T    mmap input BWT arrays (overwrite input BWTs)");
  puts("\t-Z    mmap merge arrays");
  puts("\t-B    mmap B array");
  puts("\t-D k  compute order-k deBruijn graph info (def No)");
  puts("\t-v    verbose output (more v's for more verbose)\n");
}

int main(int argc, char *argv[]) {
  extern int optind, opterr, optopt;
  extern char *optarg;
  int c, group_size;
  char *path;
  g_data g;

  // default option values
  g.verbose=0; g.solid_limit = 256;
  group_size=MAX_NUMBER_OF_BWTS;      
  g.mwXMerge = g.bwtOnly = true;
  bool hm = false;
  g.unsortedLcp = NULL;
  g.outPath = NULL;
  g.algorithm = 0;
  g.extMem = g.smallAlpha=g.mmapZ=g.mmapBWT=g.mmapB= g.lcpMerge = g.lcpCompute = false;
  g.outputDA = 0;
  g.outputSA = 0;
  g.dbOrder = 0;           // order for deBruijn graph 
  int num_threads = 0;
  while ((c=getopt(argc, argv, "vhalrxmd:p:g:A:s:o:EZTBD:S:")) != -1) {
    switch (c) 
      {
      case 'v':
        g.verbose++; break;         // verbose output (can be repeated)
      case 'o':
        g.outPath=optarg; break;    // basename for output files 
      case 'l':
        g.lcpCompute = true;        // compute lcp array from scratch 
        g.bwtOnly=false; break;  
      case 'r':
        g.lcpMerge = true;          // merge lcp values 
        g.bwtOnly=false; break;  
      case 'x':
        g.mwXMerge=false; break;      // do not use external multiway merge sort for computing lcp values
      case 'd':
        g.outputDA = atoi(optarg); break;  // output Document Array (for last iteration only) 
      case 'S':
        g.outputSA = atoi(optarg); break;  // output Suffix Array (for last iteration only) 
      case 'm':
        hm=true; break;                 // use hm algorithm 
      case 'a':
        g.smallAlpha=true; break;       // assume alphabet is small and use bwtOcc[][]
      case 'g':
        group_size = atoi(optarg);      // size for multigroup algorithm 
        break;      
      case 'A':
        g.algorithm = atoi(optarg);     // preferred algorithm 
        break;      
      case 'D':
        g.dbOrder = atoi(optarg);       // if > 1 compute order-k db graph 
        break;      
      case 's':
        g.solid_limit = atoi(optarg);   // smallest solid block 
        break;                
      case 'p':
        num_threads = atoi(optarg);     // number of consumer threads 
        break;       
      case 'E':
        g.extMem=true; break;           // use external memory (see mergegap.c)
      case 'Z':
        g.mmapZ=true; break;      // mmap merge and newmerge
      case 'B':
        g.mmapB=true; break;      // mmap B array
      case 'T':
        g.mmapBWT=true; break;    // mmap input BWTs
      case 'h':                   // usage instruction 
      case '?':
        usage(argv[0],&g);
        exit(EXIT_FAILURE);
      }
  }
  if(g.dbOrder<0 || g.dbOrder==1) { // illegal order
    printf("dbOrder must be at least 2\n");
    exit(EXIT_FAILURE);
  }
  if(g.dbOrder>1) { // valid order !=0
    if(!g.extMem) {
      printf("Option -D forces option -E\n");
      g.extMem = true;
    }
    if(g.algorithm!=128) {
      printf("Option -D forces option -A 128\n");
      g.algorithm = 128;
    }
    if(g.lcpMerge) {
      printf("Option -D incompatible with merging lcp values\n");
      exit(EXIT_FAILURE);
    }
  }
  if(num_threads <0) {
    printf("Invalid number of threads, must be non negative\n");
    exit(EXIT_FAILURE);
  }
  if(group_size<2 || group_size>MAX_NUMBER_OF_BWTS) {
    printf("Invalid group size. Must be in range [2,%llu]\n",MAX_NUMBER_OF_BWTS);
    exit(EXIT_FAILURE);
  } 
  if(hm && g.lcpCompute) {// not sure this is required, maybe changing something inside holtMcMillan() will suffice 
    printf("You cannot compute lcp values with H&M (only merge)\n");
    exit(EXIT_FAILURE);
  }
  if(g.lcpMerge && g.lcpCompute) {
    printf("You can *either* merge *or* compute lcp values\n");
    exit(EXIT_FAILURE);
  }
  if(!g.mwXMerge && !g.lcpCompute) {
    printf("Option -x can only be used with -l\n");
    exit(EXIT_FAILURE);
  }
  if(g.extMem) {
    if(hm) {
      printf("You cannot run H&M in external memory\n");
      exit(EXIT_FAILURE);
    }
    if(g.mmapBWT) { // extMem reads from the input bwt files from disk, it does not mmap them
      printf("Option -E incompatible with -T\n");
      exit(EXIT_FAILURE);
    }
    if(!g.smallAlpha) {
      printf("Option -E forces option -a\n");
      g.smallAlpha = true;
    }
  }
  if(optind+1==argc) {
    path=argv[optind];
  }
  else  {
    usage(argv[0],&g);
    exit(EXIT_FAILURE);
  }

  if(g.verbose>0) {
    puts("Command line:");
    for(int i=0;i<argc;i++)
      printf(" %s",argv[i]);
    puts("");
  }

  if(g.verbose>1) {
    puts(                    "***************************************************");
    if(g.bwtOnly)       puts("* Merge BWTs");
    else if(g.lcpMerge) puts("* Merge BWTs and LCPs");
    else                puts("* Merge BWTs and compute LCPs");
    if(hm)              puts("* Using multiround HM algorithm");
    else                puts("* Using multiround Gap algorithm");
    printf(                  "* Max #BWT %llu, Max LCP: %llu, Max mergeSize: %llu\n",
                     MAX_NUMBER_OF_BWTS, MAX_LCP_SIZE, MAX_OUTPUT_SIZE);  
    puts(                    "* Using Timo Bingmann's malloc_count to measure heap usage");
    puts(            "***************************************************");
  }

  // input check completed start measuring time 
  clock_t start = clock();
  time_t start_wc = time(NULL);
  double elapsed, elapsed_wc;

  // read BWT values. init g.numBWT, g.bwtLen[]; if !extMem also g.bws[]
  // also  remap BWTs init g->sizeOfAlpha and if g->smallAlpha is true init also g>bwtOcc[i]
  g.bwtOcc=NULL; g.sizeOfAlpha = 0;      //these two will be initialized later
  bool something_to_do = readBWTsingle(path, &g);
  
  if(g.numBwt> group_size) { // multiround computation required
    if(g.outputDA>0 || g.outputSA>0) {
      puts("Sorry, multiround DA or SA computation not supported");
      exit(EXIT_FAILURE);
    }
  }
  if(!something_to_do) {
    puts("Single BWT in input and no LCP/dBG computation. I have nothing to do!"); 
    //g.mergeLen = g.sizeOfAlpha = 1;
  }
  else {
    // init additional g fields 
    g.lcpinPath = path;
    if(g.outPath==NULL) g.outPath = path;  
    g.symb_offset = 0;
    g.blockBeginsAt = NULL; g.lcps = NULL; // prevent unintended use
    // if requested allocate space for and read lcp Values (B array)
    if(g.lcpMerge)
       initLCPmem(&g); // init g.lcps. lcp values are stored in mmapped memory 
          
    elapsed = (clock()-start)/(double)(CLOCKS_PER_SEC);
    elapsed_wc = difftime(time(NULL),start_wc);
    // printf("Partial elapsed time: %.4lf secs\n", elapsed_wc);
    // printf("Partial computing time: %.4lf secs\n", elapsed);    
    printf("Merge starting (%d bwts). Mem: %zu peak, %zu current\n", g.numBwt, malloc_count_peak(),
             malloc_count_current());
      
    // do the merging
    multiround(hm,group_size,path,&g, num_threads);
    
    // computation done the result is in g->bws[0], and g->lcps[0] (if lcpMerge) 
    // if lcpCompute the result is already in the outputfile or in the pair files 
    // and a message is printed to run mergelcp 
  
    // write content of g.bws[0] to output file (if necessary) and/or free/munmap it 
    if(g.mmapBWT) {
      int e = munmap(g.bws[0],g.mergeLen*sizeof(symbol));
      if(e) die("main (unmap bws)");
    }
    else if(!g.extMem) { 
      writeOutputBWT(&g,path,hm); // if we are working in external memory output overwrites input
      free(g.bws[0]); // deallocate all bwt's (they are contiguous)
    }
    // free other bwt related stuff
    free(g.bws);
    if(g.smallAlpha) {free(g.bwtOcc[0]); free(g.bwtOcc);}
    free(g.bwtLen);
    // free lcp related stuff
    if(g.lcpMerge) {
      int e = munmap(g.lcps[0],g.mergeLen*sizeof(lcpInt));
      if(e) die("main (unmap lcps)");
      free(g.lcps);      
    }
  }
  
  // report running times  
  elapsed = (clock()-start)/(double)(CLOCKS_PER_SEC);
  elapsed_wc = difftime(time(NULL),start_wc);
  printf("##\n");
  printf("Elapsed time: %.4lf secs\n", elapsed_wc);
  // fprintf(stderr, "%.6lf\n", elapsed_wc);
  printf("Computing time: %.4lf secs\n", elapsed);
  printf("##\n");
  printf("Size of merged BWT: "CUSTOM_FORMAT, g.mergeLen);
  printf(" Alphabet size: %d\n", g.sizeOfAlpha);
  printf("##\n");
  printf("Microseconds per symbol: %.4lf\n", elapsed_wc*1000000.0/g.mergeLen);
  //  fprintf(stderr, "%.6lf\n", elapsed_wc*1000000.0/g.mergeLen);
  printf("Peak memory allocation: %zu bytes, %.4lf bytes/symbol\n",
           malloc_count_peak(), (double)malloc_count_peak()/g.mergeLen);
  printf("##\n");

  return 0; // Done!!!
}


/**
 * create output bwt file and write the content of g->bwt[0] (already remapped 
 * and unsqueezed) to it.
 * Only called if BWT was not mmapped becaues in that case the output overwrites the input 
 * */
static void writeOutputBWT(g_data *g, char *path, bool hm)
{
  char filename[Filename_size];

  // write to output files
  if(g->verbose>1) puts("Writing output files");
  assert(!g->mmapBWT && !g->extMem);
  // open bwt output file 
  snprintf(filename,Filename_size,"%s.%s",path,hm?HM_BWT_EXT:BWT_EXT);
  FILE *f = fopen(filename,"wb");
  if(f==NULL) die(__func__);
  size_t w = fwrite(g->bws[0], sizeof(symbol), g->mergeLen, f);
  if(w!=g->mergeLen) die("Error writing final BWT");
  if(fclose(f)!=0) die(__func__);
}

