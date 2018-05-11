#ifndef THREADS_H_INCLUDED
#define THREADS_H_INCLUDED

#include "config.h"
void pc_system_init(pc_system *pc, int);
void pc_system_destroy(pc_system *pc);
void *merger(void *v);
void *squeezer(void *v);
#endif
