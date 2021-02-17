#include "file.h"

/* Returns the file extension
 */
const char *get_filename_ext(const char *filename) {                                                  
    const char *dot = strrchr(filename, '.');                                                         
    if(!dot || dot == filename) return "";                                                            
    return dot + 1;                                                                                   
}   

/*******************************************************************/

/* Changes to a working directory, where everything will be read
 * from and written to
 */ 
int file_chdir(char* dir){
  
  char* oldwd = getcwd(NULL,0);
  if (!oldwd) die(__func__);
  if (chdir(dir) == -1) die(__func__);

  free(oldwd);
return 0;
}

/*******************************************************************/

//Open a file and returns a pointer
FILE* file_open(char *c_file, const char * c_mode){
  
  FILE* f_in;
  
  f_in = fopen(c_file, c_mode);
  if (!f_in) {perror(__func__); exit(EXIT_FAILURE);}
  fseek(f_in, 0, SEEK_SET);
  
return f_in;
}


int file_close(FILE* f_in){
  
  fclose(f_in);
  if (!f_in) {perror (__func__); exit(EXIT_FAILURE);}
  
return 0;
}

/*******************************************************************/

size_t file_size(FILE* f_in){

    fseek(f_in, 0, SEEK_END);
    size_t length = ftell(f_in);
    rewind(f_in);
  
return length;
}

/*******************************************************************/

uint_t file_read(FILE* f_in){

  uint_t value;

  int e = fread(&value, sizeof(uint_t), 1, f_in);
  if(e!=1) {perror("read error"); exit(1);}

return value;
}

/*******************************************************************/

int file_write(FILE* f_out, uint_t value){

  fwrite(&value, sizeof(uint_t), 1, f_out);

return 0;
}


/*******************************************************************/

int file_write_array(FILE* f_out, int_t *A, int_t n, ssize_t size){

  int_t i;
  for(i=0; i<n; i++) fwrite(A+i, size, 1, f_out);

return 1;
}

/*******************************************************************/

char* file_load(FILE* f_in) {
  
  size_t len = 0;
  ssize_t size = 0;
  char *c_aux = NULL;
  
  /*Read one line*/
  size = getline(&c_aux, &len, f_in);
  if (size == -1) {perror(__func__); exit(EXIT_FAILURE);}

  /*Copy line to c_buffer*/
  char *c_buffer = (char*) malloc((size+1)*sizeof(char));
  strncpy(c_buffer, c_aux, size);

  c_buffer[size-1] = END_MARKER;
  c_buffer[size] = '\0';

  free(c_aux);

return c_buffer;
}

/*******************************************************************/

// read line by line
int_t* count_multiple_txt(FILE* f_in, int_t *k, uint_t chunk_size, size_t *n, int_t *chunks, ssize_t** pos) {

  if((*k)==0) (*k)=I_MAX;
  
  size_t sum=0;

  int nalloc = 128;
  int_t *K = (int_t*) malloc(nalloc*sizeof(int_t));
  (*pos) = (ssize_t*) malloc((nalloc+1)*sizeof(ssize_t));

	ssize_t tell=0;
	(*pos)[0]=tell;

  int i;
  int_t d=0;
  for(i=0; i<(*k); i++){

    size_t len = 0;
    char* c_buffer = NULL;
      
    ssize_t size = getline(&c_buffer, &len, f_in);
    if (size == -1){//EOF
      *k = i;
      free(c_buffer);
      break;   
      //return 0;
    }
    free(c_buffer);

    d++;//number of strings
    sum += size;//sum of lengths

    if(sum > chunk_size){
      //printf("K[%d] = %d\tsum = %lu\n", (*chunks), max(1,d-1), sum);
      
      K[(*chunks)++]=max(1,d-1);
      (*pos)[(*chunks)]=tell;
			if(d==1) d=0;
      sum = size; d=1;

      if((*chunks) > nalloc){
        nalloc += 128;
        K = realloc(K, sizeof(int_t) * nalloc);
        (*pos) = realloc((*pos), sizeof(ssize_t) * (nalloc+1));
      }
    }

    (*n) += (size_t) size;
		tell = ftell(f_in);
  }

  {//last chunk
    //printf("**K[%d] = %d\tsum = %lu\n", (*chunks), d, sum);
    K[(*chunks)++]=d;
  }

return K;
}

int_t* count_multiple_fastq(FILE* f_in, int_t *k, uint_t chunk_size, size_t *n, int_t *chunks, ssize_t **pos) {

  if((*k)==0) (*k)=I_MAX;
  size_t sum=0;

  int nalloc = 128;
  int_t *K = (int_t*) malloc(nalloc*sizeof(int_t));
  (*pos) = (ssize_t*) malloc((nalloc+1)*sizeof(ssize_t));
    

  int i;
  int_t d=0;
  size_t len = 0;
  char *buf = NULL;
	ssize_t tell=0;
	(*pos)[0]=tell;

  for(i=0; i<(*k); i++){

    buf = NULL; len = 0;  
    ssize_t size = getline(&buf, &len, f_in); // @'s line
    //printf("Line of size: %d\n", strlen(buf));
    free(buf);buf = NULL;
    //if (size <= 1) return 0;
    if (size == -1){//EOF
      *k = i;
      break;   
      //return 0;
    }

    char* c_buffer = NULL; len = 0; 
    size = getline(&c_buffer, &len, f_in); // read line
    //printf("Line of size: %d\n", strlen(c_buffer[i]));
    free(c_buffer);
    
    d++;//number of strings
    sum += size;//sum of lengths

    if(sum > chunk_size){
      //printf("K[%d] = %d\tsum = %lu\n", (*chunks), d-1, sum-size);
      
      K[(*chunks)++]=d-1;
      (*pos)[(*chunks)]=tell;
      sum = size; d=1;

      if((*chunks) > nalloc){
        nalloc += 128;
        K = realloc(K, sizeof(int_t) * nalloc);
        (*pos) = realloc((*pos), sizeof(ssize_t) * (nalloc+1));
      }
    }

    (*n) += (size_t)size;

    buf = NULL; len = 0;  
    size = getline(&buf, &len, f_in); // +'s line
    // printf("Line of size: %d\n", strlen(buf));
    free(buf);
    buf = NULL; len = 0;  
    size = getline(&buf, &len, f_in); // @'s line
    //printf("Line of size: %d\n", strlen(buf));
    free(buf);buf = NULL;
    //printf("Line %d OK\n", i);
		tell = ftell(f_in);
  }

  {//last chunk
    //printf("K[%d] = %d\tsum = %lu\n", (*chunks), d, sum);
    K[(*chunks)++]=d;
  }

return K;
}

int_t* count_multiple_fasta(FILE* f_in, int_t *k, uint_t chunk_size, size_t *n, int_t *chunks, ssize_t** pos) {

  if((*k)==0) (*k)=I_MAX;
  size_t sum=0;
  ssize_t size=0;

  int nalloc = 128;
  int_t *K = (int_t*) malloc(nalloc*sizeof(int_t));
  (*pos) = (ssize_t*) malloc((nalloc+1)*sizeof(ssize_t));

	ssize_t tell_cat=0, tell=0;
	(*pos)[0]=tell;

/**/
  char *buf = NULL;
  size_t len = 0;

  size = getline(&buf, &len, f_in);// first sequence
  free(buf);

  int_t d=0;
  int count=0;
  int i;
  for(i=0; i<(*k); i++){

    char *buf = NULL;

    if(i!=count) return 0;

    size_t p=0;
    buf = NULL; len = 0;  
    while((size = getline(&buf, &len, f_in))!=-1){

      if(buf[0] == '>'){
        count++;
        break;
      }
			tell_cat = ftell(f_in);
      p+=strlen(buf)-1;
    }

    d++;//number of strings

    if (size == -1){//EOF
      (*k) = i+1;
      free(buf);
      sum+=p+1;
      (*n) += p+1;
      break;   
      //return 0;
    }

    free(buf);
    size = p+1;

    sum += size;//sum of lengths

    if(sum > chunk_size){

      //printf("K[%d] = %d\tsum = %lu\n", (*chunks), d-1, sum-size);

      K[(*chunks)++]=max(1,d-1);
      (*pos)[(*chunks)]=tell;
      sum = size; d=1;

      if((*chunks) > nalloc){
        nalloc += 128;
        K = realloc(K, sizeof(int_t) * nalloc);
        (*pos) = realloc((*pos), sizeof(ssize_t) * (nalloc+1));
      }
    }

    (*n) += p+1;
		tell=tell_cat;
  }

  {//last chunk
    //printf("**K[%d] = %d\tsum = %u\n", (*chunks), d, sum);
    K[(*chunks)++]=d;
  }

return K;
}

/*******************************************************************/

// read line by line
char** load_multiple_txt(FILE* f_in, int_t k, size_t *n) {

  char **c_buffer = (char**) malloc(k*sizeof(char*));

  int_t i;
  for(i=0; i<k; i++){
    size_t len = 0;
    c_buffer[i] = NULL;
      
    ssize_t size = getline(&c_buffer[i], &len, f_in);
    if (size == -1){
      printf("K = %" PRIdN" \n", i);
      return 0;
    }
    c_buffer[i][size-1] = 0;

    (*n) += (size_t)size;
  }


return c_buffer;
}

// read sequences separeted by '@' line
char** load_multiple_fastq(FILE* f_in, int_t k, size_t *n){

  char **c_buffer = (char**) malloc(k*sizeof(char*));

  int_t i;
  for(i=0; i<k; i++){
    size_t len = 0;
    char *buf = NULL;

    ssize_t size = getline(&buf, &len, f_in); // @'s line
    //printf("Line of size: %d\n", strlen(buf));
    free(buf);buf = NULL;
    if (size <= 1) return 0;

    c_buffer[i] = NULL; len=0;
    size = getline(&c_buffer[i], &len, f_in); // read line
    //printf("Line of size: %d\n", strlen(c_buffer[i]));
    c_buffer[i][size-1] = 0;
    (*n) += (size_t)size;

    buf = NULL; len = 0;  
    size = getline(&buf, &len, f_in); // +'s line
    // printf("Line of size: %d\n", strlen(buf));
    free(buf);
    buf = NULL; len = 0;  
    size = getline(&buf, &len, f_in); // @'s line
    //printf("Line of size: %d\n", strlen(buf));
    free(buf);buf = NULL;
    //printf("Line %d OK\n", i);
    // free(buf);
  }
  return c_buffer;
}

/*******************************************************************/

// read sequences separeted by '@' line
char** load_multiple_qs(FILE* f_in, int_t k){

  char **c_buffer = (char**) malloc(k*sizeof(char*));

  int_t i;
  for(i=0; i<k; i++){
    size_t len = 0;
    char *buf = NULL;

    ssize_t size = getline(&buf, &len, f_in); // @'s line
    //printf("Line of size: %d\n", strlen(buf));
    free(buf);buf = NULL;
    if (size <= 1) return 0;

    buf = NULL; len = 0;  
    size = getline(&buf, &len, f_in); // read line
    // printf("Line of size: %d\n", strlen(buf));
    free(buf);
    buf = NULL; len = 0;  
    size = getline(&buf, &len, f_in); // +'s line
    //printf("Line of size: %d\n", strlen(buf));
    free(buf);buf = NULL;
    //printf("Line %d OK\n", i);
    // free(buf);

    c_buffer[i] = NULL; len=0;
    size = getline(&c_buffer[i], &len, f_in); // read line
    //printf("Line of size: %d\n", strlen(c_buffer[i]));
    c_buffer[i][size-1] = 0;
  }
  return c_buffer;
}

char** file_load_multiple_qs_chunks(char* c_file, int_t k, FILE *f_in) {

/* .ext
 * .fastq - strings separated by four lines
 */

  const char *type = get_filename_ext(c_file);
  char **c_buffer = NULL; // = (char**) malloc(k*sizeof(char*));

  if(strcmp(type,"fastq") == 0 || strcmp(type,"fq") == 0)
    c_buffer = load_multiple_qs(f_in, k);
  else{
    printf("Error: file not recognized (.fastq or .fq)\n");
    exit (EXIT_FAILURE);
  }

return c_buffer;
}

/*******************************************************************/

// read sequences separeted by '>' line
char** load_multiple_fasta(FILE* f_in, int_t k, size_t *n){

  char **c_buffer = (char**) malloc(k*sizeof(char*));

  char *buf = NULL;
  size_t len = 0;
  ssize_t size=0;

  size = getline(&buf, &len, f_in);// first sequence
  free(buf);

  size_t tell=0;

  int count=0;
  int_t i;
  for(i=0; i<k; i++){


    if(i!=count) return 0;

    len = 0;
    size_t nalloc = 128;
    c_buffer[i] = malloc(nalloc*sizeof(char));

    size_t p=0;
    buf = NULL; len = 0;  
    while((size = getline(&buf, &len, f_in))!=-1){

      if(buf[0] == '>'){
        count++;
        break;
      }

      if(p+len>nalloc){
        nalloc += len+128;
        c_buffer[i] = realloc(c_buffer[i], sizeof(char) * nalloc);
      }

      strcpy(&c_buffer[i][p], buf);
      p+=strlen(buf)-1;
      
      tell = ftell(f_in);
      free(buf); len=0; buf=NULL;
    }
/*
    if (size == -1){//EOF
      printf("K = %" PRIdN "\n", i);
      free(buf);
      return 0;
    }
*/

    free(buf);
    c_buffer[i][p] = 0;
    (*n) += p+1;
  }

fseek(f_in , tell, SEEK_SET);
return c_buffer;
}

/*******************************************************************/

int_t* file_count_multiple(char* c_file, int_t *k, uint_t chunk_size, int_t *chunks, size_t *n, FILE *f_in, ssize_t **pos) {

/* .ext
 * .txt   - strings per line
 * .fa .fasta - strings separated by '>' line
 * .fastq - strings separated by four lines
 */

  const char *type = get_filename_ext(c_file);
  int_t *K = NULL; 
  *pos = NULL; 

  if(strcmp(type,"txt") == 0)
    K = count_multiple_txt(f_in, k, chunk_size, n, chunks, pos);

  else if(strcmp(type,"fasta") == 0 || strcmp(type,"fa")==0 )
    K = count_multiple_fasta(f_in, k, chunk_size, n, chunks, pos);

  else if(strcmp(type,"fastq") == 0 || strcmp(type,"fq")==0 )
    K = count_multiple_fastq(f_in, k, chunk_size, n, chunks, pos);

  else{
    printf("Error: file not recognized (.txt, .fa, .fasta, .fastq, .fq)\n");
    exit(EXIT_FAILURE);
  }

  //sets the pointer to the beginning of the file
  rewind(f_in);

return K;
}

/*******************************************************************/

char** file_load_multiple_chunks(char* c_file, int_t k, size_t *n, FILE *f_in) {

/* .ext
 * .txt   - strings per line
 * .fa .fasta - strings separated by '>' line
 * .fastq - strings separated by four lines
 */

  const char *type = get_filename_ext(c_file);
  char **c_buffer = NULL; // = (char**) malloc(k*sizeof(char*));

  if(strcmp(type,"txt") == 0)
    c_buffer = load_multiple_txt(f_in, k, n);

  else if(strcmp(type,"fasta") == 0 || strcmp(type,"fa")==0 )
    c_buffer = load_multiple_fasta(f_in, k, n);

  else if(strcmp(type,"fastq") == 0 || strcmp(type,"fq")==0 )
    c_buffer = load_multiple_fastq(f_in, k, n);

  else{
    printf("Error: file not recognized (.txt, .fa, .fasta, .fastq, .fq)\n");
    exit(EXIT_FAILURE);
  }

return c_buffer;
}
/*******************************************************************/

char** file_load_multiple(char* c_file, int k, size_t *n) {

/* .ext
 * .txt   - strings per line
 * .fa .fasta - strings separated by '>' line
 * .fastq - strings separated by four lines
 */

  FILE* f_in = file_open(c_file, "rb");
  if(!f_in) return 0;

  const char *type = get_filename_ext(c_file);
  char **c_buffer = NULL; // = (char**) malloc(k*sizeof(char*));

  if(strcmp(type,"txt") == 0)
    c_buffer = load_multiple_txt(f_in, k, n);

  else if(strcmp(type,"fasta") == 0 || strcmp(type,"fa")==0 )
    c_buffer = load_multiple_fasta(f_in, k, n);

  else if(strcmp(type,"fastq") == 0 || strcmp(type,"fq")==0 )
    c_buffer = load_multiple_fastq(f_in, k, n);
  else{
    printf("Error: file not recognized (.txt, .fa, .fasta, .fastq, .fq)\n");
    exit(EXIT_FAILURE);
  }

  fclose(f_in);

return c_buffer;
}

// read sequences using the length files 
unsigned char** file_load_concat(char* c_file, char *len_file, int k, int_t *n) {

  unsigned char **c_buffer = (unsigned char**) malloc(k*sizeof(unsigned char*));
  if(c_buffer==NULL) die("out of mem");

  FILE *f = fopen(c_file,"rb");
  if(f==NULL) {perror("sequence file missing"); fprintf(stderr, "<%s> ", c_file); die("file_load_concat");}
  FILE *fl = fopen(len_file,"rb");
  if(fl==NULL) {perror("sequence file missing"); fprintf(stderr, "<%s> ", len_file); die("file_load_concat");}
  
  int i;
  for(i=0;i<k;i++) {
    uint32_t size;
    size_t r = fread(&size, 4, 1, fl);// assume big endian!
    if(r!=1) die("len file too short");
    c_buffer[i] = malloc((size+1)*sizeof(char));
    if(!c_buffer[i]) die("out of mem");
    r = fread(c_buffer[i], 1, size, f);
    if(r!=size) die("concat file too short");
    c_buffer[i][size] = 0;
    *n += size+1;
  }
  fclose(fl);
  fclose(f);
  return c_buffer;
}

/*******************************************************************/
