/******************************************************************************
 * malloc_count0.c
 *
 * dummy replacement for functions 
 *   malloc_count_current()
 *   malloc_count_peak()
 *****************************************************************************/
#include <stdlib.h>

/* user function to return the currently allocated amount of memory */
extern size_t malloc_count_current(void)
{
    return 0;
}

/* user function to return the peak allocation */
extern size_t malloc_count_peak(void)
{
    return 0;
}

