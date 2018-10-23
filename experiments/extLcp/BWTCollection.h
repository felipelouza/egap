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

#ifndef _SXSI_BWTCollection_h_
#define _SXSI_BWTCollection_h_

#include "Tools.h" // Defines ulong and uchar.
//#include "Sorting.cpp"
#include "Sorting.h"
//#include <utility> // Defines std::pair.
//#include <fstream>
#include <iostream>
using namespace std;


namespace SXSI
{
    /**
     * General interface for a bwt collection
     *
     * Class is virtual, make objects by calling
     * the static method InitBWTCollection().
     */
    class BWTCollection
    {
    public:

		std::vector <sortElement> vectTriple;  //Is is used both encoding, decoding, searching.
		//ulong seqN;  //contains a number of a sequence
		//ulong posN;  //contains the position of the last inserted symbol of the sequence seqN[i]
		//uchar pileN; //contains the number of the pile of the last inserted symbol of the sequence seqN[i]

		dataTypeNSeq nText;  //number total of texts in filename1
		dataTypeNSeq nExamedTexts;  //number total of texts in filename1
		//dataTypeNSeq middle; // number of sequence in filename1
		//dataTypeNChar middleLength; //text[middleLength] = the first letter of the second database (filename2)
		dataTypelenSeq lengthRead; //number of char in each text + $
		dataTypeNChar lengthTot;   //length of the all texts without $
		dataTypeNChar lengthTot_plus_eof; //length of the BWT
		dataTypeNChar** tableOcc; //contains the number of occurrences of each symbol
		dataTypedimAlpha alpha[SIZE_ALPHA]; //Corresponding between the alphabet, the piles and tableOcc
		dataTypedimAlpha sizeAlpha;  //number of the different symbols in the input texts
		dataTypedimAlpha *alphaInverse;  //Corresponding between alpha[i] and the symbol as char

		#if USE_QS==1
			uchar *newSymbQS;
		#endif

		#if BUILD_BCR_FROM_EGSA == 1
			dataTypeNSeq nAddedTextEGSA;
			//char c_aux[512] = "7seqsVar.fasta.7.gesa";
		#endif

		#if BUILD_EXT_MEMORY==0
			std::vector< std::vector<char> > vectVectBWT;
		#endif

		std::vector <bool> vectInsTexts;

		vector< vector< vector<dataTypeNChar> > > vectorOcc;
		vector <dataTypeNChar> numBlocksInPartialBWT;

		vector<sortElement> FirstVector, LastVector;

    static BWTCollection * InitBWTCollection(char*,char*, int);

		/**
         * Virtual destructor
         */
    virtual ~BWTCollection() { };

		virtual int buildBCR(char const *, char const *, char const *) = 0;

    virtual void storeBWTFilePartial(uchar const *, dataTypelenSeq) = 0;
		#if BUILD_LCP == 1
		    virtual void storeBWTandLCP(uchar const *, dataTypelenSeq) =0;
			virtual void storeEntireLCP(const char*) = 0;
        #endif
		virtual void storeEntireBWTFilePartial(const char*) = 0;
		#if BUILD_SA == 1
			virtual void storeEntirePairSA(const char*) = 0;
			virtual void storeEntireSAfromPairSA(const char*) = 0;
		#endif
		
		#if OUTPUT_FORMAT_EGSA == 1
			virtual int storeEGSAoutput(const char*) = 0;
			virtual int storeEGSAoutputFromEntireFiles (string input)= 0;
		#endif
		
		virtual dataTypeNChar rankManySymbolsFilePartial(FILE &, dataTypeNChar *, dataTypeNChar, uchar *)=0;
		#if BUILD_EXT_MEMORY==0
			virtual void  storeBWTIntMem(uchar const *, dataTypelenSeq) =0;
			virtual void storeEntireBWTIntMem(const char*) = 0;
			virtual dataTypeNChar rankManySymbolsIntMem(dataTypedimAlpha , dataTypeNChar *, dataTypeNChar, dataTypeNChar , uchar *) =0;
		#endif

	private:
		#if BUILD_BCR_FROM_EGSA == 1
			virtual dataTypeNChar readEGSA(char const *)=0;
		#endif
		virtual void printSegments()=0;
		virtual void printOutput(char *)=0;
		virtual void InsertNsymbols(uchar const *, dataTypelenSeq) = 0;
		virtual void InsertFirstsymbols(uchar *) = 0;   //Added/Modified/Removed 2016-02-23


		virtual int createFilePartialBWT() =0;
		virtual FILE * openWriteFilePartialBWT_0() =0;
		virtual dataTypeNChar writeFilePartial(uchar * , FILE *) =0;
		virtual FILE * openFilePartialIn(dataTypedimAlpha) =0;
		virtual FILE * openFilePartialOut(dataTypedimAlpha ) =0;
		virtual int closeFilePartial(FILE * InFile)=0;
		virtual int renameFilePartial(dataTypedimAlpha currentPile)=0;
		virtual dataTypeNChar readOnFilePartial(uchar *, dataTypeNChar , FILE * )=0;
		virtual dataTypeNChar writeOnFilePartial(uchar *, dataTypeNChar , FILE * )=0;
		virtual dataTypeNChar writeSymbolOnFilePartial(uchar , dataTypeNChar , FILE * )=0;

    protected:
        // Protected constructor; call the static function InitBWTCollection().
        BWTCollection() { };

        // No copy constructor or assignment
        BWTCollection(BWTCollection const&);
        BWTCollection& operator = (BWTCollection const&);
    };
}
#endif
