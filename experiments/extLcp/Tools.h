/**
 ** BCR is part of:
 ** BEETL: Burrows-Wheeler Extended Tool Library
 ** Documentation in: doc/BEETL.md
 **
 ** Copyright (c) 2011-2014 Illumina, Inc. **
 ** BEETL software package is
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** Citations: 
 ** 
 ** Markus J. Bauer, Anthony J. Cox and Giovanna Rosone
 ** Lightweight BWT Construction for Very Large String Collections.
 ** Proceedings of CPM 2011, pp.219-231
 ** 
 ** Markus J. Bauer, Anthony J. Cox, Giovanna Rosone and Marinella Sciortino
 ** Lightweight LCP Construction for Next-Generation Sequencing Datasets. 
 ** Proceedings of WABI 2012, pp 326-337, 2012
 
 ** Markus J. Bauer, Anthony J. Cox, Giovanna Rosone 
 ** Lightweight algorithms for constructing and inverting the BWT of string collections. 
 ** Theoretical Computer Science 483: 134-148 (2013)
 **  
 ** Anthony J. Cox, Fabio Garofalo, Giovanna Rosone, Marinella Sciortino
 ** Lightweight LCP construction for very large collections of strings. 
 ** Journal of Discrete Algorithms (2016) 37: 17-33
 **
 ** By Giovanna Rosone
 **
 **/
 
 /*
 * Collection of basic tools and defines
 */

#ifndef _TOOLS_H_
#define _TOOLS_H_

#include <iostream>
#include <fstream>

//#include <ctime>

#define TERMINATE_CHAR '\0'
#define TERMINATE_CHAR_LEN '$'

#define SIZE_ALPHA 256   

/*
#ifndef uchar
#define uchar unsigned char
#endif
#ifndef uint
#define uint unsigned int
#endif
#ifndef ulong
#define ulong unsigned long
#endif
*/

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long ulong;

#define dataTypedimAlpha uchar  //size of the alphabet (in biologic case 6 ($,A,C,G,N,T))
#define dataTypelenSeq uint   //length of the sequences (in biologic case 100)
#define dataTypeNumSeq 0    //number of sequences in the input file. If =0 -> uint
#define dataTypeNumChar 1   //numer of characters in the input file (length of the BWT) If =0 -> uint

#if dataTypeNumSeq == 1
#   define dataTypeNSeq ulong
#else
#   define dataTypeNSeq uint
#endif

#if dataTypeNumChar == 1
#   define dataTypeNChar ulong
#else
#   define dataTypeNChar uint
#endif

#define printEGSA 0

#define verboseEncode 0
#define verboseDecode 0
#define deletePartialBWT 0
#define deletePartialSA 0

#define BackByVector 1
#define BUILD_LCP 1
#define BUILD_SA 0

#define BCR_FROMCYC 0     //if BCR_FROMCYC=0 then BCR builds the cyc files before, otherwise BCR does not build the cyc files.
#define BCR_SET 1       //if BCR_SET=1 then BCR computes the EBWT (set) (one can have string of different length in input, so BCR uses the symbol TERMINATE_CHAR_LEN) else BCR_SET=0 the BWT (single sequence)
#define BCR_INPUT_IN_MEMORY 0 //if BCR_INPUT_IN_MEMORY==1 loads the input file in a string and compute the BWT of the string,  else it reads from file and computes the BWT of the reverse string.

#define BUILD_EXT_MEMORY  1 //if BUILD_EXT_MEMORY==1, BCR uses files, if BUILD_EXT_MEMORY==0, BCR uses strings

//Output of EGSA of Felipe Louze
#define BUILD_BCR_FROM_EGSA 0   //if BUILD_BCR_FROM_EGSA == 1, BCR takes in input the output of EGSA and add the symbols of new symbols. It produces the BCR output, i.e. file .pairSA, .bwt, .lcp, etc. (called in BCRexternalBWT)
//if BUILD_BCR_FROM_EGSA == 1, BCR needs of the file outputFileName.occ that contains the occurrences of each symbols in the BWT already computed.

#define OUTPUT_FORMAT_EGSA 0 //if OUTPUT_FORMAT_EGSA == 1, the output format of BCR is as the output of EGSA. BUILD_LCP and BUILD_SA and BUILD_EXT_MEMORY must be 1

#define OUTPUT_linear_SuffixArray 0

#if ((BUILD_BCR_FROM_EGSA == 1) || (OUTPUT_FORMAT_EGSA == 1))

  typedef struct{
    dataTypelenSeq  suff;
    dataTypelenSeq  lcp;
    uchar   bwt;
    dataTypeNSeq  text;
  } t_GSA;
  
  
/*  struct __attribute__((__packed__)) t_GSA {
    dataTypeNSeq  text;
    dataTypelenSeq  suff;
    dataTypelenSeq  lcp;

    uchar   bwt;
  }; */
  

#endif


//struct __attribute__((__packed__)) ElementType {
struct ElementType {
  dataTypelenSeq sa;          //It is the position in the sequence, so it goes from 0 a length read
  dataTypeNSeq numSeq;    //It is the number of the sequence.
};

class Tools
{
private:
    static double startTime;
public:
    static void StartTimer();
    static double GetTime();
    static uchar * GetRandomString(unsigned, unsigned, unsigned &);
  static uchar * GetFileContents(char *, ulong =0);
    static unsigned FloorLog2(ulong);
    static unsigned CeilLog2(ulong);
    static unsigned* MakeTable();
    static unsigned FastFloorLog2(unsigned);
};

#endif
