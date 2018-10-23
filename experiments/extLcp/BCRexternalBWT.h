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

#ifndef _BCRexternalBWT_H_
#define _BCRexternalBWT_H_

#include <map>
#include "BWTCollection.h"
#include <fstream>
#include <iostream>


class BCRexternalBWT : public SXSI::BWTCollection {
public:
    /**
     * Constructor
     */
        explicit BCRexternalBWT(char*, char*, int);
    ~BCRexternalBWT();

	int buildBCR(char const *, char const *, char const *);

  void storeBWTFilePartial(uchar const *, dataTypelenSeq);
	#if BUILD_LCP == 1
	    void storeBWTandLCP(uchar const *, dataTypelenSeq);
		void storeEntireLCP(const char*);
    #endif
	void storeEntireBWTFilePartial(const char*);
	#if OUTPUT_FORMAT_EGSA == 1
		virtual int storeEGSAoutput(const char*);
		virtual int storeEGSAoutputFromEntireFiles (string input);
	#endif
	#if BUILD_SA == 1
		void storeEntirePairSA(const char*);
		void storeEntireSAfromPairSA(const char*);
	#endif
  dataTypeNChar rankManySymbolsFilePartial(FILE &, dataTypeNChar *, dataTypeNChar, uchar *);
	#if BUILD_EXT_MEMORY==0
		void storeEntireBWTIntMem(const char*);
		void  storeBWTIntMem(uchar const *, dataTypelenSeq) ;
		dataTypeNChar rankManySymbolsIntMem(dataTypedimAlpha , dataTypeNChar *,  dataTypeNChar, dataTypeNChar , uchar *);
	#endif

private:
	#if BUILD_BCR_FROM_EGSA == 1
		dataTypeNChar readEGSA(char const *);
	#endif
	void printSegments();
	void printOutput(char *);
	void InsertNsymbols(uchar const *, dataTypelenSeq);
	void InsertFirstsymbols(uchar *); //Added/Modified/Removed 2016-02-23


	int createFilePartialBWT();
	FILE * openWriteFilePartialBWT_0();
	
	dataTypeNChar writeFilePartial(uchar * , FILE * ) ;
	FILE * openFilePartialIn( dataTypedimAlpha );
	FILE * openFilePartialOut(dataTypedimAlpha );
	int closeFilePartial(FILE * InFile);
	int renameFilePartial(dataTypedimAlpha );
	dataTypeNChar readOnFilePartial(uchar *, dataTypeNChar , FILE * );
	dataTypeNChar writeOnFilePartial(uchar *, dataTypeNChar , FILE * );
	dataTypeNChar writeSymbolOnFilePartial(uchar , dataTypeNChar , FILE * );
};

#endif
