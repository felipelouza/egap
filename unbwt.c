/* *********************************************************************
   Inversion in RAM of a bwt for a collection of sequences
   It is assumed that the multibwt uses 0 as the EOS symbol
   RAM usage: 4n bytes if n<2**32, 5n bytes for 2**32 <= n < 2**40.  
   
   Copyright (C) 2018  Giovanni Manzini  (giovanni.manzini@uniupo.it)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

   ********************************************************************* */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <time.h>
#include <limits.h>
#include <inttypes.h>


#define _BW_ALPHA_SIZE 256
typedef uint32_t sa_int;


static  int32_t get_char_from_first_column64(uint64_t i, uint64_t *first_col);

void die(const char *s) {
  fprintf(stderr,"%s\n",s);
  exit(1);
}

void perror_die(const char *s) {
  perror(s);
  exit(2);
}


#define BSIZE 1000000
void invert_bwt(char *bnam, char *tnam)
{
  FILE *b, *t;
  uint64_t occ[_BW_ALPHA_SIZE];      // emulation of first column of bwt matrix 
  uint64_t start_sa_range[_BW_ALPHA_SIZE+1];
  uint64_t i,tot;
  uint32_t *rankprev, *bwt;
  uint8_t *high8 = NULL;

  b = fopen(bnam,"r");
  if(b==NULL) perror_die("Error opening infile");
  t = fopen(tnam,"w");
  if(t==NULL) perror_die("Error opening outfile");

  // get BWT size
  long bsize = fseek(b,0,SEEK_END);
  if(bsize<0) perror_die("fseek error");
  bsize = ftell(b);
  if(bsize<0) perror_die("ftell error");
  rewind(b);
  
  if(bsize>0xFFFFFFFFFF) 
    die("BWT too large for internal memory inversion");// maximum limit is 2^40 -1
  if(bsize>0xFFFFFFFF) {
    fprintf(stderr,"Large file (>4GB): using 5n bytes ram\n");
    // we need an extra array to store the highest 8 bits 
    high8 = malloc(bsize*sizeof(*high8));
    if(high8==NULL) die("Out of mem");
  }
  bwt = rankprev = (uint32_t *) malloc(bsize*sizeof(*rankprev));
  if(rankprev==NULL) die("Out of mem");

  // clear occ
  for(i=0;i<_BW_ALPHA_SIZE;i++) occ[i]=0;  
  // read bwt from file and compute occ
  for(i=0;i<bsize;i++) {
    int e = fread(bwt+i,1,1,b);
    if(e!=1) perror_die("Error reading BWT");
    occ[bwt[i]]++; // increment count
  }
  if(fclose(b)!=0) perror_die("Error closing BWT file");
  
  // compute start sa_range  
  for(i=0,tot=0;i<_BW_ALPHA_SIZE;i++) {
    start_sa_range[i] = tot; tot+= occ[i];
  }
  start_sa_range[_BW_ALPHA_SIZE] = tot;  
  assert(tot==bsize);

  fprintf(stderr,"Computing rankprev\n");
  // bwt -> rankprev inplace
  for(i=0;i<bsize;i++) {
      // if(i%10000==0) fprintf(stderr,"%ld outof %ld\n",i,bsize);
      uint64_t rp = (start_sa_range[bwt[i]]++);
      assert(rp<bsize);
      rankprev[i] = (uint32_t) rp;
      if(high8) high8[i]=(uint8_t)(rp>>32);
  } 
  assert(start_sa_range[_BW_ALPHA_SIZE-1]==start_sa_range[_BW_ALPHA_SIZE]);
  fprintf(stderr,"Done\n");  
  // recover the original start_sa_range values
  tot=0;
  for(i=0;i<_BW_ALPHA_SIZE;i++) {
    uint64_t temp = start_sa_range[i];
    start_sa_range[i]=tot;
    tot=temp;
  }
  assert(tot==start_sa_range[_BW_ALPHA_SIZE]);

  // start recovering of sequences
  // init buffer
  long s_size=1024; long s_step = 1024;
  char *s = malloc(s_size*sizeof(*s));
  if(!s) die("out of mem");
  long s_pos = 0;
  for(int i=0;i<occ[0];i++) {
    if(i%10000==0) fprintf(stderr,"Recovering sequecence %d\n",i);
    s_pos=0;
    long old_rank = i;
    while(true) {
      // move to first column
      uint64_t rank = rankprev[old_rank];
      if(high8!=NULL) rank += ((uint64_t) high8[old_rank]) << 32;
      int32_t c = get_char_from_first_column64(rank, start_sa_range);
      //fprintf(stderr,"%ld --> %d\n",rank,c); 
      if(c==0) break; //end of sequence
      if(s_pos>=s_size) {
        assert(s_pos==s_size);
        s_size += s_step;
        s = realloc(s,s_size*sizeof(*s));
      }
      assert(s_pos<s_size);  
      s[s_pos++] = c;
      old_rank = rank;
    }
    if(s_pos==0)
      die("Empty string in BWT file");
    // ---- reverse string ----
    for(long k=0;k<s_pos/2;k++) {
     char c = s[k]; s[k] = s[s_pos-1-k]; s[s_pos-1-k]=c;
    }
    // write string
    // fprintf(t,">\n");
    int e = fwrite(s,1,s_pos,t);
    if(e!=s_pos) perror_die("Error writing ouput sequence");
    fprintf(t,"\n");
  }
  fprintf(stderr,"Recovered %ld sequences\n",occ[0]);
  if(fclose(t)!=0) perror_die("Error closing output file");      
  free(rankprev);
  if(high8) free(high8);
}


// given an index returns the corresponding char in the bwt matrix
// doing a binary search in the first column representation
static int32_t get_char_from_first_column64(uint64_t index, uint64_t *first_col)
{
  // binary search
  int32_t med,lo=0,hi=_BW_ALPHA_SIZE-1;
  // invariant: first_col[lo]<= index < first_col[hi+1]
  assert(first_col[lo]<=index && index<first_col[hi+1]);
  while(lo<hi) {
    med = (lo+hi+1)/2;
    if(index < first_col[med])
      hi = med-1;
    else
      lo = med;
  }
  assert(lo==hi);  
  return lo;
  #if 0  
  int32_t i;
  // sequential search
  for(i=0;i<_BW_ALPHA_SIZE;i++)
    if(first_col[i]<=index && index < first_col[i+1]) {
      assert(lo==i);
      return i;
    }
  _bwtext_fatal_error(__func__);
  return 0; 
  #endif 
}


int main(int argc, char *argv[])
{
  char *bnam=NULL, *tnam=NULL; 

  /* ------------- read options from command line ----------- */
  if(argc==3) {
    bnam=argv[1]; tnam = argv[2];
  }
  else {
    fprintf(stderr,"Usage:  %s infile outfile\n\n", argv[0]);
    fprintf(stderr,"Takes as input a multi-bwt that uses 0x0 as the EOS symbol and recovers\n");
    fprintf(stderr,"the original sequences writing them to <outfile> separated by a newline\n\n");
    exit(1);
  }
  
  // do the inversion
  invert_bwt(bnam,tnam);
  return 0;
}
