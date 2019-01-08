#ifndef FILE_H
#define FILE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "utils.h"


typedef struct _file_pair{
	int_t		docs;
	ssize_t pos;
} file_pair;

/*******************************************************************/
int file_chdir(char* dir);

FILE* file_open(char *c_file, const char * c_mode);
int file_close(FILE* f_in);

size_t file_size(FILE* f_in);

int file_write(FILE* f_out, uint_t array);
uint_t file_read(FILE* f_in);

int file_write_array(FILE* f_out, int_t *A, int_t n);

char* file_load(FILE* f_in) ;
char** file_load_multiple(char* c_file, int k, size_t* n) ;

char** file_load_multiple_chunks(char* c_file, int_t k, size_t* n, FILE *f_in);
int_t* file_count_multiple(char* c_file, int_t *k, uint_t chunk_size, int_t* chunks, size_t* n, FILE *f_in, ssize_t **pos);

unsigned char** file_load_concat(char* c_file, char *len_file, int k, int_t *n);
/*******************************************************************/



#endif
