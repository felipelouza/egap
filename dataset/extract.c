// vim: noai:ts=2:sw=2

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <limits.h>


#define min(a,b) ((a) < (b) ? (a) : (b))

#define MB 1048576

/**********************************************************************/

/* Returns the file extension
 */
const char *get_filename_ext(const char *filename) {
        const char *dot = strrchr(filename, '.');                                                         
        if(!dot || dot == filename) return "";                                                            
        return dot + 1;                                                                                   
}

/**********************************************************************/

/* Returns the file without extension
 */
const char *get_filename_without_ext(const char *filename) {                                              
        char *dot = strrchr(filename, '.');
        if(!dot || dot == filename) return filename;

        *dot = '\0';
        return filename;
}

/*******************************************************************/

// read line by line
int load_multiple_txt(FILE* f_in, FILE* f_out, size_t size, size_t len) {

        size_t l=0, n=0, sum=0, k=0, max=0;
        char *buf=NULL;
        while((n = getline(&buf, &l, f_in))!=-1){

                if(len){
                        if(n<len){//discard
                                free(buf);
                                l= 0; buf = NULL;
                                continue;
                       }
                        else{ //cut
                                buf[len] = '\n';
                                buf[len-1] = '\0';
                                n=len;
                        }
                }

                k++;
                if(n>max) max=n;
                sum += n;
                if(sum > size){
                        k--;
                        printf("%lu\t%.2lf\t%lu\n", k, (double)((sum-n)/k), max);
                        free(buf);
                        return 0;
                }
                fprintf(f_out, "%s", buf);
                free(buf);
                l= 0; buf = NULL;
        }

        printf("%lu bytes\n", sum);
        printf("%lu\t%.2lf\t%lu\n", k, (double)(sum/k), max);
        free(buf);

return 1;
}

/**********************************************************************/

// read sequences separeted by '@' line
int load_multiple_fastq(FILE* f_in, FILE* f_out, size_t size, size_t len){

        size_t l=0, n=0, sum=0, k=0, max=0;
        char *buf1 = NULL, *buf2, *buf3, *buf4;

        while(getline(&buf1, &l, f_in)!=-1){ // @'s line

                l= 0; buf2 = NULL;
                n = getline(&buf2, &l, f_in); // sequence

                if(len){
                        if(n<len){//discard

                                l= 0; buf3 = NULL;
                                getline(&buf3, &l, f_in); // +'s line
                                l= 0; buf4 = NULL;
                                getline(&buf4, &l, f_in); // @'s line

                                free(buf1); free(buf2); free(buf3);     free(buf4);
                                l= 0; buf1 = NULL;

                                continue;
                        }
                        else{ //cut
                                buf2[len-1] = '\n';
                                buf2[len] = '\0';
                                n=len;
                        }
                }

                k++;

                if(n>max) max=n;
                l= 0; buf3 = NULL;
                getline(&buf3, &l, f_in); // +'s line
                l= 0; buf4 = NULL;
                getline(&buf4, &l, f_in); // @'s line
    sum += n;
                if(sum > size){
                        printf("%lu bytes\n", sum-n);
                        free(buf1); free(buf2); free(buf3);     free(buf4);
                        k--;
                        if(k) printf("%lu\t%.2lf\t%lu\n", k, (double)((sum-n)/k), max);
                        return 0;
                }

                fprintf(f_out, "%s", buf1);
                fprintf(f_out, "%s", buf2);
                fprintf(f_out, "%s", buf3);
                fprintf(f_out, "%s", buf4);

                free(buf1); free(buf2); free(buf3);     free(buf4);
                l= 0; buf1 = NULL;
        }

        printf("%lu bytes\n", sum);
        if(k) printf("%lu\t%.2lf\t%lu\n", k, (double)(sum/k), max);
        free(buf1);

return 1;
}

/**********************************************************************/
// read sequences separeted by '>' line
int load_multiple_fasta(FILE* f_in, FILE* f_out, size_t size, size_t len){
      size_t l=0, n=0, sum1=0, sum2=0, k=0, max=0, curr=0;

        char *c_buffer = NULL;

        char *buf1 = NULL, *buf2=NULL;
        getline(&buf1, &l, f_in);// first >'s line

        l= 0;
        size_t p=0;
        int nalloc = MB;
        c_buffer = malloc(nalloc*sizeof(char));

//f_out = stdout;
        while((n=getline(&buf2, &l, f_in))!=-1){ // >'s line

        //      sum2+=n-1;
                curr+=n-1;

                strcpy(&c_buffer[p], buf2);
                p+=strlen(buf2)-1;

                free(buf2);
                l= 0; buf2=NULL;
                while((n=getline(&buf2, &l, f_in))!=-1){

                        if(buf2[0] == '>'){

                                if(len) curr=min(len,curr);

                                if(sum2+curr<=size){

                                        if(len){
                                                if(curr<len){//discard
                                                        k--;
                                                }
                                                else{

                                                        c_buffer[len-1] = '\n';
                                                        c_buffer[len] = '\0';
                                                        curr=len;
                                                        sum2+=curr;
                                                        fprintf(f_out, "%s", buf1); //header
                                                        fprintf(f_out, "%s", c_buffer); //sequence
                                                        if(curr>max) max=curr;

                                                }
                                        }
                                        else{
                                                sum2+=curr;
                                                fprintf(f_out, "%s", buf1); //header
                                               fprintf(f_out, "%s", c_buffer); //sequence
                                                if(curr>max) max=curr;
                                        }

                                        sum1=sum2;
                                        free(buf1);
                                        curr=0;
                                        buf1 = (char*) malloc(nalloc*sizeof(char));
                                }
                                else{
                                        printf("%lu bytes\n", sum1);
                                        if(k) printf("%lu\t%.2lf\t%lu\n", k, (double)((sum1)/k), max);
                                        free(c_buffer);
                                        free(buf1);
                                        free(buf2);
                                        return 0;
                                }

                                k++;

                                strcpy(buf1, buf2);
                                p=0; nalloc=MB;
                                free(c_buffer);
                                c_buffer = malloc(nalloc*sizeof(char));

                                free(buf2);
                                buf2=NULL; l=0;

                                break;
                        }

//                      sum2+=n-1;
                        curr+=n-1;

                        if(p+l>nalloc){
                                nalloc += l+MB;
                                c_buffer = realloc(c_buffer, sizeof(char) * nalloc);
                        }

                        strcpy(&c_buffer[p], buf2);
                        p+=strlen(buf2)-1;

                        free(buf2);
                        l= 0; buf2 = NULL;
                }

                free(buf2);
                l= 0; buf2 = NULL;
        }

  if(len) curr=min(len,curr);
        if(sum2<=size){
        k++;
                if(len){
                        if(curr<len){//discard
                                k--;
                        }
                        else{

                                c_buffer[len-1] = '\n';
                                c_buffer[len] = '\0';
                                curr=len;
                                sum2+=curr;
                                fprintf(f_out, "%s", buf1); //header
                                fprintf(f_out, "%s", c_buffer); //sequence
                                if(curr>max) max=curr;

                        }
                }
                else{
                        sum2+=curr;
                        fprintf(f_out, "%s", buf1); //header
                        fprintf(f_out, "%s", c_buffer); //sequence
                        if(curr>max) max=curr;
                }
                sum1=sum2;
        }

        free(c_buffer);
        free(buf1);
        free(buf2);

        printf("%lu bytes\n", sum1);
        if(k) printf("%lu\t%.2lf\t%lu\n", k, (double)((sum1)/k), max);

return 1;
}

/*******************************************************************/

int file_load_multiple(char* c_in, size_t size, size_t len) {

/* .ext
 * .txt   - strings per line
 * .fasta - strings separated by '>' line
 * .fastq - strings separated by four lines
 */

        FILE* f_in = fopen(c_in, "rb");
        if(!f_in) return 0;

  printf("INPUT:\t%s\n", c_in);

        const char *type = get_filename_ext(c_in);
  char c_out[500];
        const char *name = get_filename_without_ext(c_in);
  //sprintf(c_out, "%s.%zu.%s", name, size, type);
  if(len)
          sprintf(c_out, "%s.%lu.%s", name, len, type);
        else
          sprintf(c_out, "%s.X.%s", name, type);

        FILE* f_out = fopen(c_out, "wb");
        if(!f_out) return 0;

  printf("OUTPUT:\t%s\n", c_out);

        if(strcmp(type,"txt") == 0)
                load_multiple_txt(f_in, f_out, size, len);
        else if(strcmp(type,"fasta") == 0)
                load_multiple_fasta(f_in, f_out, size, len);
        else if(strcmp(type,"fastq") == 0)
                load_multiple_fastq(f_in, f_out, size, len);
        else{
                printf("Error: file not recognized (.txt, .fasta, .fastq)\n");
                return 0;
        }

//printf("k\tavg(len)\tmax(len)\n");

        fclose(f_in);
        fclose(f_out);

return 1;
}

/*******************************************************************/

void usage(char *name){
  printf("\n\tUsage: %s [options] FILE N len \n\n",name);
  puts("Output:\tFILE.X contains strings of total size < N.\n");
  puts("If len>0, it discards strings < len and cut strings > len .\n");
  puts("Available options:");
  puts("\t-h\tthis help message");
  //puts("\t-m RAM  available memory (in KB) to be used by gSACA-k algorithm");
  //puts("\t-g D\tLCP in gap format D bytes per entry (ext: .D.lcp)");
  exit(EXIT_FAILURE);
}

/**********************************************************************/

int main(int argc, char **argv) {

  extern char *optarg;
  extern int optind, opterr, optopt;

        int c=0, time=0, verbose=0;

  char *c_file=NULL;
  size_t size=1024, len=100;

  while ((c=getopt(argc, argv, "vth")) != -1) {
    switch (c)
      {
//      case 'k':
//        heap_size=atoi(optarg); break;  // compute LCP and output in Gap format
      case 'v':
        verbose++; break;
      case 't':
                                time++; break;
      case 'h':
        usage(argv[0]); break;       // show usage and stop
      //case 'm':
        //RAM=(size_t)atoi(optarg)*KB; break;
      case '?':
        exit(EXIT_FAILURE);
      }
  }
  free(optarg);

  if(optind+3==argc) {
    c_file=argv[optind++];
    //size= (size_t) atoi(argv[optind++]);
    sscanf(argv[optind++], "%zu", &size);
    sscanf(argv[optind++], "%zu", &len);
  }
  else{
                usage(argv[0]);
        }

        if(size==0) size = SSIZE_MAX;
        file_load_multiple(c_file, size, len);

}

