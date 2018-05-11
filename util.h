#ifndef UTIL_H_INCLUDED
#define UTIL_H_INCLUDED

#include "config.h"
#include "alphabet.h"

void alloc0_B_array(g_data *g);
void free_B_array(g_data *g);
void alloc_merge_arrays(g_data *g);
void free_merge_arrays(g_data *g);
void alloc_merge_array(g_data *g);
void free_merge_array(g_data *g);
void alloc_merge_array16(g_data *g);
void free_merge_array16(g_data *g);
void alloc_array32(g_data *g);
void free_array32(g_data *g);

bool readBWTsingle(char *path, g_data *);
void initLCPmem(g_data *g);

void open_unsortedLCP_files(g_data *g);
void close_unsortedLCP_files(g_data *g);
void writeLcp(customInt k, uint32_t lcp, g_data *g);
void writeLcp_EOF(uint64_t size, g_data *g);

void array_clear(customInt* array, customInt size, customInt value);
void array_copy(customInt *dest, customInt *src, customInt size);

void testPrintMergeColor(customInt, palette *array, customInt size);
void testPrintI(char *name, customInt param, customInt *array, customInt size);
void testPrintLcp(char *name, customInt param, lcpInt *array, customInt size);

void writeBWTandLCP(FILE *bwtout, FILE *lcpout, FILE *lcpin[], g_data *g);

void mergeBWTandLCP(g_data *g, bool lastRound);
void mergeBWTandLCP256(g_data *g, bool lastRound);
void mergeBWT128(g_data *g, bool lastRound);
void mergeBWT128ext(g_data *g, bool lastRound);
void mergeBWT8(g_data *g, bool lastRound);

void check_g_data(g_data *g);
void die(const char* where);


#endif
