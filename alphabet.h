#ifndef ALPHABET_H_INCLUDED
#define ALPHABET_H_INCLUDED

#include "config.h"
int floorlog(int n);
int intpow(int b, int e);
void init_freq(symbol *b, customInt n, customInt *freq);
void init_freq_no0(symbol *b, customInt n, customInt *freq);
int init_alpha_maps(g_data *g); // symbol *b, customInt n);
int alpha_reduce(int c);
int alpha_enlarge(int c);
// void remap_string(symbol *b, customInt n, customInt *freq);

int squeezable(int n);
int intpow(int b, int e);
void ktuple_init_freq(symbol *b, customInt n, customInt *freq);
void ktuple_init(int a, int k);
bool ktuple_has0(int s);
bool ktuple_tail0(int s);
int ktuple_decode0(int s);
int ktuple_headlen(int s);
int ktuple_last(int s);
void ktuple_info(FILE *f);
void ktuple_print(FILE *f, int s);
int ktuple_lcp(int x, int y);
// void ktuple_remap_string(symbol *bwt, customInt n, customInt *freq, int a, int k);
#endif
