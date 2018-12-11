#include "utils.h"

/**********************************************************************/

void time_start(time_t *t_time, clock_t *c_clock){

	*t_time = time(NULL);
	*c_clock =  clock();
}

double time_stop(time_t t_time, clock_t c_clock){

	double aux1 = (clock() - c_clock) / (double)(CLOCKS_PER_SEC);
	double aux2 = difftime (time(NULL),t_time);
	
	printf("CLOCK = %lf TIME = %lf\n", aux1, aux2);
	
	return aux1;
}


/**********************************************************************/
void die(const char* where) {

	printf("Error at %s: %s.\n",where,errno ? strerror(errno) : "errno not set");
	exit(errno?errno:1);
}
 
void dies(const char* where, char* format, ...) {

	int_t err = errno;
	va_list val;

	printf("Error at %s: %s.\n",where,errno ? strerror(errno) : "errno not set");
	va_start(val,format);
	vprintf(format,val);
	va_end(val);
	printf("\n");

	exit(err ? err : 1);
}
/**********************************************************************/
int_t print_char(char* A, int_t n){

	int_t i;
	for(i=0; i<n; i++)	
		printf("%" PRIdN ") %c\n", i, A[i]);

	printf("\n");

return 0;
}
/**********************************************************************/
int_t print_int(int_t* A, int_t n){

	int_t i;
	for(i=0; i<n; i++)	
		printf("%" PRIdN ") %" PRIdN "\n", i, A[i]);

	printf("\n");

return 0;
}
/**********************************************************************/
int_t min_range(int_t* A, int_t l, int_t r){


	if(r>l)return 0;

	printf("[l, r] = [%" PRIdN ", %" PRIdN "]\n", l, r);

	int_t min = INT_MAX;
	int_t i;
	for(i=l; i<=l; i++)
		min = (A[i]<min?A[i]:min);

	printf("min = %" PRIdN "\n", min);

return min;
}
/*******************************************************************/
int_t* cat_int(unsigned char** R, int k, int_t *n){

	(*n)++; //add 0 at the end

	int_t i, j;
	int_t l=0;
	int_t *str_int = (int_t*) malloc((*n)*sizeof(int_t));

	for(i=0; i<k; i++){
		int_t m = strlen((char*)R[i]);
		for(j=0; j<m; j++){
			//removes symbols > 255
			if(R[i][j]+1<256) str_int[l++] = R[i][j]+(k+1);
			else (*n)--;
		}
//		for(j=0; j<m; j++)
//			str_int[l++] = R[i][j]+(k+1);
		str_int[l++] = i+1; //add $_i as separator
	}
	
	str_int[l++]=0;
        if(*n>l){
		str_int = (int_t*) realloc(str_int, (l)*sizeof(int_t));
		printf("N = %" PRIdN "\n", l);
	}
	*n = l;

return str_int;
}
/*******************************************************************/
#if REVERSE_SCHEME==1
unsigned char* cat_char_rev(unsigned char** R, int k, size_t *n){

	(*n)++; //add 0 at the end

	int_t i, j;
	int_t l=0;
	unsigned char *str = (unsigned char*) malloc((*n)*sizeof(unsigned char));

	for(i=0; i<k; i++){
		int_t m = strlen((char*)R[i]);
		//removes empty strings
		if(m==0){
			(*n)--;
			continue;
		}
		for(j=m-1; j>=0; j--){
			//removes symbols > 255
			if(R[i][j]+1<256) str[l++] = R[i][j]+1;
			else (*n)--;
		}
		str[l++] = 1; //add 1 as separator
	}

	str[l++]=0;
  if(*n>l){
		str = (unsigned char*) realloc(str, (l)*sizeof(unsigned char));
		printf("N = %" PRIdN "\n", l);
	}
	*n = l;

return str;
}
/*******************************************************************/
#elif	REVERSE_SCHEME==2
unsigned char* cat_char_rev(unsigned char** R, int k, size_t *n){

	(*n)++; //add 0 at the end

	int_t i, j;
	int_t l=0;
	unsigned char *str = (unsigned char*) malloc((*n)*sizeof(unsigned char));

	for(i=k-1; i>=0; i--){

		int_t m = strlen((char*)R[i]);
		//removes empty strings
		if(m==0){
			(*n)--;
			continue;
		}
		for(j=m-1; j>=0; j--){
			//removes symbols > 255
			if(R[i][j]+1<256) str[l++] = R[i][j]+1;
			else (*n)--;
		}
		str[l++] = 1; //add 1 as separator
	}

	str[l++]=0;
  if(*n>l){
		str = (unsigned char*) realloc(str, (l)*sizeof(unsigned char));
		printf("N = %" PRIdN "\n", l);
	}
	*n = l;
return str;
}
#endif
/*******************************************************************/
unsigned char* cat_char(unsigned char** R, int k, size_t *n){

	(*n)++; //add 0 at the end

	int_t i, j;
	int_t l=0;
	unsigned char *str = (unsigned char*) malloc((*n)*sizeof(unsigned char));

	for(i=0; i<k; i++){
		int_t m = strlen((char*)R[i]);
		//removes empty strings
		if(m==0){
			(*n)--;
			continue;
		}
		for(j=0; j<m; j++){
			//removes symbols > 255
			if(R[i][j]+1<256) str[l++] = R[i][j]+1;
			else (*n)--;
		}
		str[l++] = 1; //add 1 as separator
	}

	str[l++]=0;
  if(*n>l){
		str = (unsigned char*) realloc(str, (l)*sizeof(unsigned char));
		printf("N = %" PRIdN "\n", l);
	}
	*n = l;

return str;
}

/**********************************************************************/

static void swap2(void *x, void *y, size_t l) {
   char *a = x, *b = y, c;
   while(l--) {
      c = *a;
      *a++ = *b;
      *b++ = c;
   }
}

static void sort(char *array, size_t size, int (*cmp)(void*,void*), int begin, int end) {
   if (end > begin) {
      void *pivot = array + begin;
      int l = begin + size;
      int r = end;
      while(l < r) {
         if (cmp(array+l,pivot) <= 0) {
            l += size;
         } else if ( cmp(array+r, pivot) > 0 )  {
            r -= size;
         } else if ( l < r ) {
            swap2(array+l, array+r, size);
         }
      }
      l -= size;
      swap2(array+begin, array+l, size);
      sort(array, size, cmp, begin, l);
      sort(array, size, cmp, r, end);
   }
}

void qsort2(void *array, size_t nitems, size_t size, int (*cmp)(void*,void*)) {
   sort(array, size, cmp, 0, nitems*size);
}

/**********************************************************************/

unsigned char map(unsigned char c){

  if(c==0) return 0; //$ and # symbols 

  switch(c){
    case 'A': 
    case 'a': c = 1<<5;
              break;

    case 'C':
    case 'c': c = 2<<5;
              break;

    case 'G': 
    case 'g': c = 3<<5;
              break;

    case 'T': 
    case 't': c = 4<<5;
              break;

    //including N
    default:  c = 5<<5;
              break;
  }

return c;
}

unsigned char unmap(unsigned char c){

  c = c>>5;

  switch(c){
    case 1: c = 'a';
            break;

    case 2: c = 'c';
            break;

    case 3: c = 'g';
            break;

    case 4: c = 't';
            break;

    //including N
    default:  c = 'n';
              break;
  }

return c;
}

unsigned char rle(unsigned char c, unsigned char run){

  run--; //run=0 means run=1
  c = map(c);

return c|run;
}

/**********************************************************************/

int bwt(int_t sa, unsigned char* str){

    int c;

    if(sa==0 || str[sa-1]==1)
      c = 0;
    else
      c = str[sa-1]-1;

return c;
}

/**********************************************************************/

