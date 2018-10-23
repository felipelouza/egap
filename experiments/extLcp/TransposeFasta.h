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
 
 #ifndef TRANPOSEFASTA_INCLUDED
#define TRANPOSEFASTA_INCLUDED


#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
//#include <stdlib.h>

#include "Tools.h"
#include <string.h>    


#define BUFFERSIZE 1024

//#define CYCLENUM 100

//#define uchar unsigned char

using std::string;
using std::vector;


class TransposeFasta
{
public:
    TransposeFasta();
    ~TransposeFasta();

  bool convert (const string& input, char const * fileOutput, const string& output );       //it reads a file in fasta format and builds cyc files
  bool convertFromCycFile(const string& input, char const * fileOutput);        //it reads the file input for find the length and the number of sequences. Does not build cyc files.
  bool convert1Sequence(char const * filename1) ;

  dataTypelenSeq lengthRead;    //Lenght of each text
  dataTypeNChar lengthTexts;   //Total length of all texts without $-symbols

  dataTypeNSeq nSeq;   //number total of texts in filename1
  dataTypeNChar freq[256];  //contains the distribution of the symbols. It is useful only for testing. It depends on the #characters

  #if BCR_INPUT_IN_MEMORY==1
    //std::wstring strInput;
    uchar * strInput;
  #endif

    #if BUILD_BCR_FROM_EGSA == 1
      dataTypeNSeq nAddedTextEGSA_transp;
  #endif


private:

  bool findLengthNseq( const string& input, const string& fileOutput);


//    FILE* outputFiles_[CYCLENUM];
//    uchar buf_[CYCLENUM][BUFFERSIZE];

  ulong readln(char* s, int n, FILE* iop);

};



#endif
