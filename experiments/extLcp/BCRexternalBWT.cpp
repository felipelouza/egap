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

#include "BCRexternalBWT.h"
#include "BWTCollection.h"
#include "TransposeFasta.h"
#include "Tools.h"
#include "Timer.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <assert.h>
#include <sys/stat.h>
#include <math.h>
#include <string.h>    

//#include "Sorting.h"
#define SIZEBUFFER 1024
#define DIMBLOCK 1024

//using std::vector;
using namespace std;
using SXSI::BWTCollection;

////////////////////////////////////////////////////////////////////////////
// Class BCRexternalBWT

/**
 * Constructor inits
 */
BCRexternalBWT::BCRexternalBWT(char *file1, char *fileOutput, int mode)
{
  for (dataTypeNSeq j = 0 ; j < nText; j++) {
    vectTriple[j].posN = j+1;  // position of the suffix (1-based)
    vectTriple[j].seqN = j;   // number of sequence
    vectTriple[j].pileN = 0;    //The first symbols are in $-pile
    #if BUILD_LCP == 1
      vectTriple[j].lcpCurN = 0;
      vectTriple[j].lcpSucN = 0;
    #endif
  }

  const char * fileOut = "cyc.";
  if (mode == 0) {
        time_t start,end;
        double dif;
        time (&start);

    std::cerr << "Start BCR encode\n";
    if (BUILD_LCP==0)
      std::cerr << "Compute the BWT by using the BCR (BCRLCPSA)\n";
    else
      std::cerr << "Compute the BWT and the LCP by using the BCR (BCRLCPSA)\n";
    if (BUILD_SA==1)
      std::cerr << "Compute also the SA by using the BCR (BCRLCPSA)\n";
    
    if (OUTPUT_FORMAT_EGSA==1) {
      if ((BUILD_LCP==1) && (BUILD_SA==1) && (BUILD_EXT_MEMORY==1))
        std::cerr << "The output format of BCR is as the output of EGSA.\n";
      else
        std::cerr << "The output format of BCR is as the output of EGSA. BUILD_LCP and BUILD_SA and BUILD_EXT_MEMORY must be 1 (see Tools.h).\n";
    }
    
    int result = -1;

    #if BCR_SET == 1
       result = buildBCR(file1, fileOutput, fileOut);
    #else
          //For a sequence, we read the file one symbol at time
          //we assume that the input file is the inverse file
      result = buildBCR(file1, fileOutput, fileOut);
    #endif
    assert (result != '1');

    #if BUILD_EXT_MEMORY==1
      #if OUTPUT_FORMAT_EGSA == 0
        //Store the entire BWT from sizeAlpha-files
        storeEntireBWTFilePartial(fileOutput);
      #else
        if ((BUILD_LCP == 1) && (BUILD_SA==1)){
          storeEGSAoutput(fileOutput);        
        }
      #endif
    #else
      storeEntireBWTIntMem(fileOutput);
    #endif

    #if ((OUTPUT_FORMAT_EGSA==0) && (BUILD_LCP == 1))
      storeEntireLCP(fileOutput);
    #endif

    #if ((OUTPUT_FORMAT_EGSA==0) && (BUILD_SA==1))
      storeEntirePairSA(fileOutput);
    
      //Optional
      #if OUTPUT_linear_SuffixArray ==  1
        storeEntireSAfromPairSA(fileOutput);
      #endif
    #endif

    #if verboseEncode==1
      printOutput(fileOutput);
    #endif

    /*
    std::cerr << "Removing the auxiliary input file (cyc files)\n";
    char *filename;
    filename = new char[strlen(fileOut)+sizeof(dataTypelenSeq)*8];

    // delete output files
    for(dataTypelenSeq i=0;i<lengthRead;i++ ) {
      sprintf (filename, "%s%u.txt", fileOut, i);
      if (remove(filename)!=0)
        std::cerr << filename <<" BCRexternalBWT: Error deleting file" << std::endl;
    }
    */
    time (&end);
        dif = difftime (end,start);
        std::cerr << "Start builBCR (including the writing the cyc files) " << start << " seconds\n";
        std::cerr << "End   builBCR (including the writing the cyc files) " << end << " seconds\n";
        std::cerr << "builBCR (including the writing the cyc files) tooks " << dif << " seconds\n";
  }
  else
    std::cerr << "Mode Error \n";

}



//Computes the rank function and returns the number of symbols that it read.
//The rank function computes the number char less than the symbol c from the starting position (startPos) in the BWT to the position pos (endPos) in the BWT.
//Here, we compute the number of occurrences of each symbol from from the starting position (startPos) in the BWT to the position pos (endPos) in the BWT.
//The startPos is the position of the File pointer InFileBWT, the endPos depends on toRead
//In the original definition of the rank, startPos corresponds to the position 1 and endPos corresponds to the previous symbol.
//Here, we work by using \sigma partial BWTs.
//toRead is the number of symbols that I have to read before to find the symbol in B corresponding to the symbol in F.
dataTypeNChar BCRexternalBWT::rankManySymbolsFilePartial(FILE & InFileBWT, dataTypeNChar *counters, dataTypeNChar toRead, uchar *foundSymbol)
{
  dataTypeNChar numchar, cont=0;  //cont is the number of symbols already read!
  uchar *buffer = new uchar[SIZEBUFFER];

  //it reads toRead symbols from the fp file (Partial BWT)
  while (toRead > 0) {            //((numchar!=0) && (toRead > 0)) {
    if (toRead <= SIZEBUFFER) {    //Read toRead characters
      numchar = fread(buffer,sizeof(uchar),toRead,&InFileBWT);
      assert(numchar == toRead); // we should always read/write the same number of characters
      *foundSymbol = buffer[numchar-1];     //The symbol of the sequence k.  It is the symbol in the last position in the partial BWT that we have read.

      if (*foundSymbol == TERMINATE_CHAR) {
        std::cerr << "--> Rank toRead=" << (unsigned int)toRead << " foundSymbol= " << (unsigned int)(*foundSymbol) << std::endl;
      }

    }
    else {   //Read sizebuffer characters
      numchar = fread(buffer,sizeof(uchar),SIZEBUFFER,&InFileBWT);
      assert(numchar == SIZEBUFFER); // we should always read/write the same number of characters
    }

    //For each symbol in the buffer, it updates the number of occurrences into counters
    for (dataTypeNChar r=0; r<numchar; r++)
      counters[alpha[(unsigned int)buffer[r]]]++;    //increment the number of letter symbol into counters


    cont   += numchar;  //number of read symbols
    toRead -= numchar;  //number of remaining symbols to read
    if ((numchar == 0) && (toRead > 0)) {  //it means that we have read 0 character, but there are still toRead characters to read
      std::cerr << "rankManySymbolsFilePartial: read 0 character, but there are still " << toRead << " characters to read  " << std::endl;
      exit (EXIT_FAILURE);
    }
  }
  delete [] buffer;

  return cont;
}

#if BUILD_EXT_MEMORY==0
dataTypeNChar BCRexternalBWT::rankManySymbolsIntMem(dataTypedimAlpha currentPile, dataTypeNChar *counters, dataTypeNChar alreadyRead, dataTypeNChar toRead, uchar *foundSymbol)
{
  dataTypeNChar cont;  //cont is the number of symbols already read!

  //it reads toRead symbols from the BWT string (Partial BWT)
  //For each symbol in the buffer, it updates the number of occurrences into counters
  assert(toRead <= vectVectBWT[currentPile].size());
  for (cont=alreadyRead; cont<alreadyRead+toRead; cont++)   {
    counters[alpha[(unsigned int)vectVectBWT[currentPile][cont]]]++;    //increment the number of letter symbol into counters
  }

  *foundSymbol = vectVectBWT[currentPile][alreadyRead+toRead-1];     //The symbol of the sequence k.  It is the symbol in the last position in the partial BWT that we have read.

  return toRead;
}
#endif



int BCRexternalBWT::buildBCR(char const * file1, char const * fileOutput, char const * fileOut)
{
  #if  BUILD_EXT_MEMORY == 1    //BCR uses the internal memory for the BWT partial
    std::cout << "BCR uses the external memory for the BWT partial" << endl;
  #else
    std::cout << "BCR uses the internal memory for the BWT partial" << endl;
    #if ((BUILD_LCP == 1) || (BUILD_SA == 1) )
        std::cout << "Error! For computing LCP and SA, then BCR_INPUT_IN_MEMORY must be 0." << std::endl;
        exit (EXIT_FAILURE);
    #endif
  #endif

  #if BCR_SET == 0        //Build BCR for 1 sequence
    std::cout << "BCR of 1 sequence" << endl;


    #if BCR_FROMCYC==1
      std::cerr << "Error: BCR_FROMCYC==1 (cyc files in input) is not implemented! (see BCR_FROMCYC in Tools.h)" << endl;
      exit (EXIT_FAILURE);
    #endif

    #if BUILD_BCR_FROM_EGSA==1
      std::cerr << "Error: BUILD_BCR_FROM_EGSA==1 is not implemented! (see BUILD_BCR_FROM_EGSA in Tools.h)" << endl;
      exit (EXIT_FAILURE);
    #endif

    #if BCR_FROMCYC==1
      std::cerr << "Error: BCR_FROMCYC==1 is not implemented! (see BCR_FROMCYC in Tools.h)" << endl;
      exit (EXIT_FAILURE);
    #endif

    #if BCR_INPUT_IN_MEMORY==1    // BCR reads from string
      std::cout << "loads the input file in a string and compute the BWT of the string." << std::endl;
    #else               // BCR reads from file
      std::cout << "reads from file and computes the BWT of the reverse string." << std::endl;
    #endif
  #else
    std::cout << "BCR of multi-sequences" << endl;
    #if BCR_INPUT_IN_MEMORY==1    // BCR reads from string
      std::cout << "Error BCR_SET == 0, so that BCR_INPUT_IN_MEMORY must be 0." << std::endl;
      exit (EXIT_FAILURE);
    #endif

  #endif

  #if BUILD_BCR_ALTERNATE == 0
    std::cout << "Lexicographic order" << endl;
  #else
    std::cout << "Alternate lexicographic order" << endl;
  #endif




  std::cerr << "Build BCR before of transpose\n";
  Timer timer;
  cerr << "TIMER start buildBCR " << timer<< endl;

  time_t startBuildBCR,endBuildBCR, startTranspose, endTranpose;
    double difTranspose, difBuildBCR;
    time (&startTranspose);

  int res=false;
  TransposeFasta trasp;

  #if BCR_SET == 0        //Build BCR for 1 sequence
    //Build the cyc files
    //Read 1 sequence
    bool controllo = trasp.convert1Sequence(file1);

    if (controllo == false) {  //Error in the reading
      std::cerr << "Reading Error \n";
      exit (EXIT_FAILURE);
    }
  #else
    #if BCR_FROMCYC==0
      //Build the cyc files
      cerr << "Builds cyc files and the builds the BCR " << endl;

      res = trasp.convert( file1, fileOutput, fileOut);

      if (res == false) {  //Error in the reading
        std::cerr << "Error in transpose (builds cyc files)! \n";
        exit (EXIT_FAILURE);
      }

    #else
      //cyc files in input
      cerr << "Builds the BCR (no cyc files)" << endl;
      res = trasp.convertFromCycFile(file1, fileOutput);
      if (res == false) {  //Error in the reading
        std::cerr << "Error in transpose (no build cyc files)! \n";
        exit (EXIT_FAILURE);
      }
    #endif



  #endif


  nText = trasp.nSeq;
  lengthRead = trasp.lengthRead;
  lengthTot = trasp.lengthTexts;

  #if (BUILD_BCR_FROM_EGSA == 1)
    nAddedTextEGSA = trasp.nAddedTextEGSA_transp;
    std::cerr << "Number of added sequences from EGSA: " << nAddedTextEGSA << "\n";
  #endif


  sizeAlpha=0;
  for (dataTypedimAlpha i = 0; i < SIZE_ALPHA-1; ++i)
    if (trasp.freq[i] > 0) {
      alpha[i] = sizeAlpha;
      sizeAlpha++;
    }
  if (trasp.freq[SIZE_ALPHA-1] > 0) {
    alpha[SIZE_ALPHA-1] = sizeAlpha;
    sizeAlpha++;
  }
    std::cerr << "We supposed that the symbols in the input file are:\n";
  for (dataTypedimAlpha i = 0; i < SIZE_ALPHA-1; ++i)
    if (trasp.freq[i] > 0)
      std::cerr << (unsigned int)i << " " << i << " " << trasp.freq[i] << " " << (unsigned int)alpha[i] << "\n";
  if (trasp.freq[SIZE_ALPHA-1] > 0)
    std::cerr << (unsigned int)SIZE_ALPHA-1 << " " << SIZE_ALPHA-1 << " " << trasp.freq[SIZE_ALPHA-1] << " " << (unsigned int)alpha[SIZE_ALPHA-1] << "\n";
  
  std::cerr << "dataTypedimAlpha: sizeof(type size of alpha): " << sizeof(dataTypedimAlpha) << "\n";
  std::cerr << "dataTypelenSeq: sizeof(type of seq length): " << sizeof(dataTypelenSeq) << "\n";
  std::cerr << "dataTypeNSeq: sizeof(type of #sequences): " << sizeof(dataTypeNSeq) << "\n";
  
  #if OUTPUT_linear_SuffixArray ==  1
    std::cerr << "dataTypeNChar: sizeof(type of #characters): " << sizeof(dataTypeNChar) << "\n";
    std::cerr << "RAM Buffer for bufferNChar: " << SIZEBUFFER * sizeof(dataTypeNChar) << "\n";
  #endif
  
  std::cerr << "\nElementType: sizeof(type of pairSA element): " << sizeof(ElementType) << "\n";
  std::cerr << "sortElement: sizeof(type of sortElement): " << sizeof(sortElement) << "\n";

  std::cerr << "\nSizeAlpha: " << (unsigned int)sizeAlpha << "\n";
  std::cerr << "Number of sequences: " << nText << "\n";
  std::cerr << "Length of the longest sequence: " << (unsigned int)lengthRead << "\n\n";
  std::cerr << "Total length (without $): " << lengthTot << "\n";

    lengthTot_plus_eof = lengthTot+nText;
  std::cerr << "Total length (with $): " << lengthTot_plus_eof << "\n";

  std::cerr << "\nRAM for input symbols: " << nText << "\n";
  std::cerr << "RAM for BWT-LCP-PairSA: " << nText * sizeof(sortElement) << "\n";
  std::cerr << "RAM for PairSA: " << nText * sizeof(ElementType) << "\n";
  
  
  std::cerr << "Size Buffer: " << SIZEBUFFER << "\n";
  
  std::cerr << "RAM Buffer for PairSA: " << SIZEBUFFER * sizeof(ElementType) << "\n";
  
  std::cerr << "RAM Buffer for bufferLCP: " << SIZEBUFFER * sizeof(dataTypelenSeq) << "\n";

  std::cerr << "\nRAM Buffer for minLCPcur: " << (unsigned int)sizeAlpha * sizeof(dataTypelenSeq) << "\n";
  std::cerr << "RAM Buffer for minLCPcurFound: " << (unsigned int)sizeAlpha * sizeof(bool) << "\n";
  std::cerr << "RAM Buffer for minLCPsucToFind: " << (unsigned int)sizeAlpha * sizeof(bool) << "\n";
  std::cerr << "RAM Buffer for minLCPsuc: " << (unsigned int)sizeAlpha * sizeof(dataTypelenSeq) << "\n";
  std::cerr << "RAM Buffer for minLCPsucText: " << (unsigned int)sizeAlpha * sizeof(dataTypeNSeq) << "\n";

  dataTypeNChar totRAMinMB = nText +        //symbols to insert - step i
    (nText * sizeof(sortElement)) +       //vectTriple
    (SIZEBUFFER  * sizeof(ElementType)) +     //Buffer for GSA 
    (SIZEBUFFER * sizeof(dataTypelenSeq));    //Buffer for LCP


//    (nText * sizeof(ElementType)) +
//    (nText * sizeof(ElementType)) +
    
  #if OUTPUT_linear_SuffixArray ==  1
    totRAMinMB = totRAMinMB + (SIZEBUFFER * sizeof(dataTypeNChar));
  #endif
  
  
  std::cerr << "\nTotal RAM in MB for BWT-LCP-PairSA: " << (long) totRAMinMB / 1024 << "\n";
  
  totRAMinMB = totRAMinMB + 
      ((unsigned int)sizeAlpha * sizeof(dataTypelenSeq)) + 
      ((unsigned int)sizeAlpha * sizeof(bool)) + 
      ((unsigned int)sizeAlpha * sizeof(bool)) +
      ((unsigned int)sizeAlpha * sizeof(dataTypelenSeq)) +
      ((unsigned int)sizeAlpha * sizeof(dataTypeNSeq));
  std::cerr << "Total RAM  in MB for BWT-LCP-PairSA, computing LCP: " << (long) totRAMinMB / 1024 << "\n";

  time (&endTranpose);
    difTranspose = difftime (endTranpose,startTranspose);
  time (&startBuildBCR);

  char *filename;
  filename = new char[strlen(fileOut)+sizeof(dataTypelenSeq)*8];

  #if BUILD_EXT_MEMORY==1
    createFilePartialBWT();
    std::cerr << "\nPartial File name for input: " << fileOut <<" \n\n";
  #else
    vectVectBWT.resize(sizeAlpha);
  #endif



  uchar *newSymb = new uchar[nText];
  vectInsTexts.resize(nText);
  for (dataTypeNSeq j = 0 ; j < nText; j++) {
      vectInsTexts[j] = 0;
    }
  //vectTriple.resize(nText);
  tableOcc = new dataTypeNChar*[sizeAlpha];
  //Counting for each pile: $-pile, A-pile, C-pile, G-pile, N-pile, T-pile
  for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
      tableOcc[j] = new dataTypeNChar[sizeAlpha];
    }
  for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++)
    for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++)
      tableOcc[j][h]=0;

  //lengthTot = 0;  //Counts the number of symbols

  std::cerr << "\nFirst symbols: "<< "j= "<< 0 <<" - symbols in position " << (long) lengthRead-1 << "\n";
  //cout << "Starting iteration " << 0 << ", time now: " << timer.timeNow();
  //cout << "Starting iteration " << 0 << ", usage: " << timer << endl;

  static FILE *InFileInputText;
  dataTypeNChar num;
  #if BCR_SET == 0        //Build BCR for 1 sequence
    sprintf (filename, "%s", file1);
    //std::cerr << "The filename is ''" << filename <<"'' file." << std::endl;
    InFileInputText = fopen(filename, "rb");
    if (InFileInputText==NULL) {
      std::cerr << filename <<" : Error opening " << std::endl;
      exit (EXIT_FAILURE);
    }

    #if BCR_INPUT_IN_MEMORY==1    // BCR reads from string
      newSymb[0] = trasp.strInput[lengthRead-1];      //lengthRead does not contain the $.
    #else               // BCR reads from file
      num = fread(newSymb,sizeof(uchar),nText,InFileInputText);
      assert( num == nText); // we should always read the same number of characters
    #endif
  #else
    sprintf (filename, "%s%u.txt", fileOut, lengthRead-1);
    InFileInputText = fopen(filename, "rb");
    if (InFileInputText==NULL) {
        std::cerr << filename <<" : Error opening " << std::endl;
        exit (EXIT_FAILURE);
    }
    num = fread(newSymb,sizeof(uchar),nText,InFileInputText);
    assert( num == nText); // we should always read the same number of characters
    //lengthTot += num;   //Increment the number of chars
    fclose(InFileInputText);



  #endif



  std::cerr << "\n"<< "j= "<< 0 <<" - symbols in position " << (long) lengthRead - 1<< "\n";
  InsertFirstsymbols(newSymb);

  #if verboseEncode==1
    printSegments();
  #endif


  #if BUILD_BCR_FROM_EGSA == 1  //BCR concatenate the output of EGSA.
    std::cerr << "BUILD_BCR_FROM_EGSA "  << BUILD_BCR_FROM_EGSA <<" \n";
    dataTypeNChar numCharReadInEGSA = readEGSA(file1);
    std::cerr << "From EGSA we have added "  << numCharReadInEGSA <<" symbols.\n";

  #endif


  std::cerr << "\nInserting symbols after the first position." << endl;
  //cout << "Starting iteration, time now: " << timer.timeNow();
  //cout << "Starting iteration, usage: " << timer << endl;

  //maxLengthRead-2
  for (dataTypelenSeq t = lengthRead-2 ; t > 0; t--) {
    #if verboseEncode==1
      std::cerr << "\n"<< "j= "<< (long) lengthRead - t - 1 <<" - symbols in position " << (long) t << "\n";
    #endif

    //if (t % 1000000 == 0) {
    //  std::cerr << "\n"<< "j= "<< (long) lengthRead - t - 1 <<" - symbols in position " << (long) t << "\n";
      //cout << "Starting iteration " << (long) lengthRead - t - 1  << ", time now: " << timer.timeNow();
      //cout << "Starting iteration " << (long) lengthRead - t - 1 << ", usage: " << timer << endl;
    //}

      //To insert the symbol from position m-3 to position 1
      //The last inserted symbol is in position i+1 (or it is newSymb[j]),
      //the next symbol (to insert) is in position i

    #if BCR_SET == 0        //Build BCR for 1 sequence
      newSymb[0]='\0';

      #if BCR_INPUT_IN_MEMORY==1    // BCR reads from string
        newSymb[0] = trasp.strInput[t];     //t does not contain the $.

        if (newSymb[0] == TERMINATE_CHAR) {
          std::cerr << "--> build---\n"<< "j= "<< (long) lengthRead - t - 1 <<" - symbols in position " << (long) t << "\n";
        }
      #else
        num = fread(newSymb,sizeof(uchar),nText,InFileInputText);
        assert( num == nText); // we should always read the same number of characters
      #endif

        //std::cerr << "The read symbol is: " << newSymb[0] << ".\n";
    #else
      sprintf (filename, "%s%u.txt", fileOut, t);
      InFileInputText = fopen(filename, "rb");
      if (InFileInputText==NULL) {
        std::cerr << filename <<" : Error opening " << std::endl;
        exit (EXIT_FAILURE);
      }
      num = fread(newSymb,sizeof(uchar),nText,InFileInputText);
      assert( num == nText); // we should always read the same number of characters
      //lengthTot += num;   //Increment the number of chars
      fclose(InFileInputText);


    #endif

    /*
    if ((newSymb[0] != '\r') && (newSymb[0] != '\n') && (newSymb[0] != '\0')) {
        //cerr << "block sequence " << newSymb << "." << endl;
          InsertNsymbols(newSymb, t);
        }
        else
           cerr << "Ignored symbol, in j= " << lengthRead - t << ", is " <<(unsigned int)newSymb[0] << "." << endl;
    */

    InsertNsymbols(newSymb, t+1);

    #if verboseEncode==1
      printSegments();
    #endif

  }

  //The last inserted symbol is in position 1 (or it is newSymb[j]),
  //the next symbol (to insert) is in position 0
  std::cerr << "\n"<< "j= "<< (long) lengthRead - 1 <<" - symbols in position " << 0 << "\n";
  //cout << "Starting iteration " << (long) lengthRead - 1 << ", time now: " << timer.timeNow();
    //cout << "Starting iteration " << (long) lengthRead - 1 << ", usage: " << timer << endl;

  #if BCR_SET == 0        //Build BCR for 1 sequence

    #if BCR_INPUT_IN_MEMORY==1    // BCR reads from string
      newSymb[0] = trasp.strInput[0];     //t does not contain the $.
    #else
      num = fread(newSymb,sizeof(uchar),nText,InFileInputText);
      assert( num == nText); // we should always read the same number of characters
      //std::cerr << "The read symbol is: " << newSymb[0] << ".\n";
      fclose(InFileInputText);
    #endif
  #else
    sprintf (filename, "%s%u.txt", fileOut, 0);
    InFileInputText = fopen(filename, "rb");
    if (InFileInputText==NULL) {
        std::cerr << filename <<" : Error opening " << std::endl;
        exit (EXIT_FAILURE);
    }
    num = fread(newSymb,sizeof(uchar),nText,InFileInputText);
    assert( num == nText); // we should always read the same number of characters
    fclose(InFileInputText);
    /*
    std::cerr << filename  << std::endl;
    std::cerr << "read number is " << num << "\n";
    for (dataTypeNSeq j = 0 ; j < nText; j++)
      std::cerr << newSymb[j];
    std::cerr << ".\n";
    */
    //lengthTot += num;   //Increment the number of chars

  #endif

  InsertNsymbols(newSymb, 1);
  //The last inserted symbol is in position 0 (or it is newSymb[j]),
  #if verboseEncode==1
    printSegments();
  #endif


    //the next symbol (to insert) is in position m-1, that is, I have to insert the symbols $
    std::cerr << "\n"<< "j= "<< (long) lengthRead <<" - symbols in position " << (long) lengthRead  << ". Inserting $=" << (unsigned int)TERMINATE_CHAR << "=" << TERMINATE_CHAR << " symbol\n\n";
    //cout << "Starting iteration " << (long) lengthRead << ", time now: " << timer.timeNow();
    //cout << "Starting iteration " << (long) lengthRead << ", usage: " << timer << endl;
  for (dataTypeNSeq j = 0 ; j < nText; j++) {
    #if (BCR_SET==1)
      if ((newSymb[j] == TERMINATE_CHAR) || (newSymb[j] == TERMINATE_CHAR_LEN))
        newSymb[j] = TERMINATE_CHAR_LEN;
      else
        newSymb[j] = TERMINATE_CHAR;
    #else
      newSymb[j] = TERMINATE_CHAR;
    #endif



  }
      //std::cerr << "\n";

  InsertNsymbols(newSymb, 0);
  #if verboseEncode==1
    cerr << "\nAll new characters have been inserted, usage: " << timer << endl;
    printSegments();
  #endif


  // to delete those:
  delete [] newSymb;


  #if BCR_SET == 0        //Build BCR for 1 sequence
    #if BCR_INPUT_IN_MEMORY==1    // BCR reads from string
      delete [] trasp.strInput;
    #endif
  #endif


//  vectTriple.~vector<sortElement>();
  /*
  delete [] seqN;
  delete [] pileN;
  delete [] posN;
  */
  std::cerr << std::endl;
  //std::cerr << "The input text is long " << lengthTot << std::endl;
  dataTypeNChar numCharInTable = 0;
  for (dataTypedimAlpha r = 0; r < sizeAlpha; r++) {
    for (dataTypedimAlpha t = 0; t < sizeAlpha; t++) {
      numCharInTable += tableOcc[r][t];
    }
  }
  std::cerr << "In tableOcc, there are " << numCharInTable << " letters" << std::endl;

  for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
    delete [] tableOcc[j];
    tableOcc[j] = NULL;
  }
  delete [] tableOcc;

  delete[] filename;
  delete[] alphaInverse;
  //------------------

  time (&endBuildBCR);
    difBuildBCR = difftime (endBuildBCR,startBuildBCR);

  std::cerr << "\nStart Preprocessing " << startTranspose << " seconds\n";
    std::cerr << "End   Preprocessing " << endTranpose << " seconds\n";
    std::cerr << "Preprocessing tooks " << difTranspose << " seconds\n";

    std::cerr << "Start builBCR " << startBuildBCR << " seconds\n";
    std::cerr << "End   builBCR " << endBuildBCR << " seconds\n";
    std::cerr << "builBCR tooks " << difBuildBCR << " seconds\n\n";

  std::cerr << "Preprocessing + builBCR tooks " << difTranspose + difBuildBCR << " seconds\n";

  //delete trasp;

  return true;
}

void BCRexternalBWT::InsertFirstsymbols(uchar * newSymb)
{

  char *filenameOut = new char[12];
  char *filename = new char[8];

  const char *ext = ".aux";
  sprintf (filename, "bwt_%d",0);
  sprintf (filenameOut,"%s%s",filename,ext);

  #if BUILD_EXT_MEMORY==1
    FILE *OutFileBWT = openWriteFilePartialBWT_0();
  #endif



  sortElement tripla;
  nExamedTexts=0;
  for (dataTypeNSeq j = 0 ; j < nText; j++) {

    #if (BCR_SET==1)  //The input is a set
    if (newSymb[j] != TERMINATE_CHAR_LEN) {     //we don't insert it into partial BWT_0, we can remove it, we shift it
      tableOcc[0][alpha[(unsigned int)newSymb[j]]]++;       //counting the number of occurrences in BWT of the $-pile
      newSymb[nExamedTexts] = newSymb[j];
      if (j != nExamedTexts)
        newSymb[j] = TERMINATE_CHAR_LEN;



      tripla.posN = nExamedTexts+1;  // position of the suffix (1-based)
      tripla.seqN = j;    // number of the sequence
      tripla.pileN = 0;    //The first symbols are in $-pile
      #if BUILD_LCP == 1
        tripla.lcpCurN = 0;   //it should be -1?
        tripla.lcpSucN = 0;   //it should be -1?
      #endif
      vectTriple.push_back(tripla);
      vectInsTexts[j]=1;
      nExamedTexts++;
    }
  #else
    tableOcc[0][alpha[(unsigned int)newSymb[j]]]++;       //counting the number of occurrences in BWT of the $-pile
    newSymb[nExamedTexts] = newSymb[j];


    tripla.posN = nExamedTexts+1;  // position of the suffix (1-based)
    tripla.seqN = j;    // number of the sequence
    tripla.pileN = 0;    //The first symbols are in $-pile
    #if BUILD_LCP == 1
      tripla.lcpCurN = 0;   //it should be -1?
      tripla.lcpSucN = 0;   //it should be -1?
    #endif
    vectTriple.push_back(tripla);
    vectInsTexts[j]=1;
    nExamedTexts++;
  #endif

  }

  std::cerr << std::endl;
  #if verboseEncode==1
    std::cerr << "First step:" << std::endl;
    std::cerr << "U:  ";
    for (dataTypeNSeq j = 0 ; j < nExamedTexts; j++) {
      std::cerr << newSymb[j] << " ";
    }
    std::cerr << std::endl;
    std::cerr << "Q  ";
    for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
      std::cerr << (unsigned int)vectTriple[g].pileN << " ";
    }
    std::cerr << std::endl;
    std::cerr << "P  ";
    for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
      std::cerr << vectTriple[g].posN  << " ";
    }
    std::cerr << std::endl;
    std::cerr << "N  ";
    for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
      std::cerr << vectTriple[g].seqN  << " ";
    }
    std::cerr << std::endl;
    #if BUILD_LCP == 1
      std::cerr << "C  ";             //LCP current
      for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
        std::cerr << (unsigned int)vectTriple[g].lcpCurN  << " ";
      }
      std::cerr << std::endl;
      std::cerr << "S  ";                     //LCP successive
      for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
        std::cerr << (unsigned int)vectTriple[g].lcpSucN  << " ";
      }
      std::cerr << std::endl;
    #endif
  #endif
  //Store newSymb into $-pile BWT

  //dataTypeNChar num = fwrite (newSymb, sizeof(uchar), nText , OutFileBWT);
  //assert( num == nText); // we should always read the same number of characters
  //We don't store TERMINATE_CHAR_LEN symbol (the strings can have different length)
  //dataTypeNChar num = fwrite (newSymb, sizeof(uchar), nExamedTexts , OutFileBWT);
  #if BUILD_EXT_MEMORY==1
    dataTypeNChar num = writeFilePartial(newSymb, OutFileBWT);
    assert( num == nExamedTexts); // we should always read the same number of characters
    fclose(OutFileBWT);
  #else
    for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
      vectVectBWT[0].push_back (newSymb[g]);
    }
  #endif


  //if (num != nText)
  //  std::cerr << "the written characters is not equal to number of the texts" << num << " and "<< nText <<"\n";


  #if BUILD_LCP == 1
    dataTypelenSeq *vectLCP = new dataTypelenSeq[nExamedTexts];
    //vectLCP.resize(nText);
    for (dataTypeNSeq j = 0 ; j < nExamedTexts; j++) {
      vectLCP[j]=0;
    }
    static FILE *OutFileLCP;                  // output and input file LCP;
    sprintf (filename, "lcp_%d",0);
    sprintf (filenameOut,"%s%s",filename,ext);
    OutFileLCP = fopen(filenameOut, "wb");
    if (OutFileLCP==NULL) {
      std::cerr << "LCP file $: Error opening " << std::endl;
      exit (EXIT_FAILURE);
    }
    num = fwrite (vectLCP, sizeof(dataTypelenSeq), nExamedTexts , OutFileLCP);
    assert( num == nExamedTexts); // we should always read the same number of integers
    //vectLCP.clear();
    fclose(OutFileLCP);

    //Creates one file for each letter in the alphabet for LCP. From 1 to sizeAlpha-1
    for (dataTypedimAlpha i = 1; i < sizeAlpha; i++) {
      sprintf (filename, "lcp_%d", i);
      sprintf (filenameOut,"%s%s",filename,ext);
      OutFileLCP = fopen(filenameOut, "wb");
      if (OutFileLCP==NULL) {
        std::cerr << "LCP file " << (unsigned int)i <<" : Error opening " << std::endl;
        exit (EXIT_FAILURE);
      }
    }
    fclose(OutFileLCP);
    delete [] vectLCP;
  #endif

  //Do we want to compute the extended suffix array (position and number of sequence)?
  #if (BUILD_SA == 1)   //To store the SA
    sprintf (filename, "sa_%d",0);
    sprintf (filenameOut,"%s%s",filename,ext);
    FILE *OutFileSA = fopen(filenameOut, "wb");
    if (OutFileSA==NULL) {
      std::cerr << "SA file: Error opening: " << filenameOut << " (SA file $)" << std::endl;
      exit (EXIT_FAILURE);
    }

    ElementType *newEle = new ElementType[nExamedTexts];
    for (dataTypeNSeq j = 0 ; j < nExamedTexts; j++) {
      // lengthRead --> Length of the longest sequence
      newEle[j].sa=lengthRead;      //(posSymb + 1) % (lengthRead + 1);
      newEle[j].numSeq= vectTriple[j].seqN;
    }
    //Store into $-pile SA
    num = fwrite (newEle, sizeof(ElementType), nExamedTexts , OutFileSA);
    if (num != nExamedTexts)
      std::cerr << "The written characters is not equal to number of the texts in SA" << num << " and "<< nExamedTexts <<"\n";
    assert(num == nExamedTexts);

    fclose(OutFileSA);

    delete[] newEle;
    //------------------

    //Creates one file for each letter in the alphabet. From 1 to sizeAlpha-1
    for (dataTypedimAlpha i = 1; i < sizeAlpha; i++) {
      sprintf (filename, "sa_%d", i);
      sprintf (filenameOut,"%s%s",filename,ext);

      OutFileSA = fopen(filenameOut, "wb");
      if (OutFileSA==NULL) {
        std::cerr << "SA file " << (unsigned int)i <<" : Error opening " << std::endl;
        exit (EXIT_FAILURE);
      }
      fclose(OutFileSA);
    }
  #endif

  delete [] filenameOut;
  delete [] filename;


  #if verboseEncode==1
    #if BUILD_LCP == 1
      dataTypeNChar numchar;
      static FILE *InFileLCP;                  // output and input file LCP;
      char *filenameInLCP = new char[12];
      char *filenameLCP = new char[8];
      //const char *ext = ".aux";
      dataTypelenSeq *bufferLCP = new dataTypelenSeq[SIZEBUFFER];
      dataTypedimAlpha mmm=0;
      while (mmm < sizeAlpha) {
        numchar=sprintf (filenameLCP, "lcp_%d", mmm);
        numchar=sprintf (filenameInLCP,"%s%s",filenameLCP,ext);
        //printf("===currentPile= %d\n",mmm);
        InFileLCP = fopen(filenameInLCP, "rb");
        for (dataTypeNChar g = 0 ; g < SIZEBUFFER; g++)
          bufferLCP[g] = 0;
        numchar = fread(bufferLCP,sizeof(dataTypelenSeq),SIZEBUFFER,InFileLCP);
        std::cerr << "L[" << (unsigned int)mmm << "]:\t";
        if (numchar==0)
          std::cerr  << "empty";
        else
          for (dataTypeNSeq g = 0 ; g < numchar; g++)
            std::cerr  << (unsigned int)bufferLCP[g] << " ";
        std::cerr  << "\n";
        fclose(InFileLCP);
        mmm++;
      }
      delete [] bufferLCP;
      delete [] filenameInLCP;
      delete [] filenameLCP;
    #endif
  #endif


}



void BCRexternalBWT::InsertNsymbols(uchar const * newSymb, dataTypelenSeq posSymb)
{
  #if BUILD_EXT_MEMORY==1
    static FILE *InFileBWT;                  // output and input file BWT;
  #endif
  char *filenameIn = new char[12];
  char *filename = new char[8];
  //const char *ext = ".aux";

  dataTypeNChar *counters = new dataTypeNChar[sizeAlpha];  //it counts the number of each symbol into the i-Pile-BWT
  for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
      counters[i]=0;


  dataTypeNChar toRead = 0;
  //Find the positions of the new symbols
  dataTypeNSeq j = 0;
  //Insert the symbols beloging to old strings
  while (j < nExamedTexts) {
    //std::cerr << "+++j= " << j << " vectInsTexts[j]= " << vectInsTexts[j] << " newSymb[j]= " << newSymb[j]  << std::endl;
    //if ((vectInsTexts[j]==0) && (newSymb[j] == TERMINATE_CHAR_LEN))  //nothing to do
    #if (BCR_SET==1)
      if ((vectInsTexts[vectTriple[j].seqN]==1) && (newSymb[vectTriple[j].seqN] != TERMINATE_CHAR_LEN)) { //The new symbol have to be inserted in some y-pile with y>0
    #else
      if (vectInsTexts[vectTriple[j].seqN]==1) { //The new symbol have to be inserted in some y-pile with y>0
    #endif

      dataTypedimAlpha currentPile = vectTriple[j].pileN;

      #if BUILD_EXT_MEMORY==1
        InFileBWT = openFilePartialIn(currentPile );
      #endif


      dataTypeNSeq k=j;
      //For each pile, we have a different counter of characters
      for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
        counters[i]=0;
      dataTypeNChar cont = 0;   //number of the read symbols
      uchar foundSymbol;
      dataTypeNChar numberRead=0;
      #if BUILD_EXT_MEMORY==0
        dataTypeNChar alreadyRead=0;
      #endif  
      while ((k< nExamedTexts) && (vectTriple[k].pileN == currentPile)) {
         #if (BCR_SET==1)
         if (newSymb[vectTriple[k].seqN] != TERMINATE_CHAR_LEN) {
         #endif
          //if (verboseEncode == 1)
          //  std::cerr << "j-1: Q["<<k<<"]=" << (unsigned int)vectTriple[k].pileN << " P["<<k<<"]=" << (dataTypeNChar)vectTriple[k].posN << " N["<<k<<"]=" << (dataTypeNSeq)vectTriple[k].seqN << "\t";

          //std::cerr << "--k= " << k << " pileN[k]= " << pileN[k] << " posN[k]= " << posN[k] <<  " seqN[k]= " << seqN[k] << std::endl;
          //For any character (of differents sequences) in the same pile
          foundSymbol = '\0';
          //cont is the number of symbols already read!
          toRead = vectTriple[k].posN - cont;

          #if BUILD_EXT_MEMORY==1
            numberRead = rankManySymbolsFilePartial(*InFileBWT, counters, toRead, &foundSymbol);
          #else
            numberRead = rankManySymbolsIntMem(currentPile, counters, alreadyRead, toRead, &foundSymbol);
            alreadyRead += numberRead;
          #endif


          if ((foundSymbol == TERMINATE_CHAR) || (foundSymbol == '\0')) {
            std::cerr << "--> dopo rank New posN[k]=" << (unsigned int)vectTriple[k].posN << " New pileN[k]=" << (unsigned int)vectTriple[k].pileN << " posSymb= " << posSymb << " foundSymbol= " << (unsigned int)foundSymbol << std::endl;
          }

          assert (toRead == numberRead);
          cont += numberRead;
          //I have to update the value in vectTriple[k].posN, it must contain the position of the new symbol

          #if BUILD_BCR_ALTERNATE == 0  //We compute straightforward order the BWT of the sequences

            vectTriple[k].posN = counters[alpha[(unsigned int)foundSymbol]];
            //std::cerr << "--counters[alpha[(unsigned int)foundSymbol]=" << (unsigned int)counters[alpha[(unsigned int)foundSymbol]] <<std::endl;
            //std::cerr << "--New posN[k]=" << (unsigned int)posN[k] <<std::endl;
            /*
            if (verboseEncode == 1)
              std::cerr << "\nInit New P["<< k <<"]= " << vectTriple[k].posN <<std::endl;
            */
            for (dataTypedimAlpha g = 0 ; g < currentPile; g++) {  //I have to count in each pile g= 0... (currentPile-1)-pile
              vectTriple[k].posN = vectTriple[k].posN + tableOcc[g][alpha[(unsigned int)foundSymbol]];

              //std::cerr << "--New posN[k]=" << (unsigned int)posN[k] << " tableOcc[g][alpha[(unsigned int)symbol]] " << tableOcc[g][alpha[(unsigned int)symbol]] <<std::endl;
              /*
              if (verboseEncode == 1) {
                std::cerr << "g= " << (unsigned int)g << " symbol= " << (unsigned int)foundSymbol << " alpha[symbol]= "<< (unsigned int)alpha[(unsigned int)foundSymbol] <<std::endl;
                std::cerr << "Add New posN[k]=" << vectTriple[k].posN << " tableOcc[g][alpha[(unsigned int)symbol]] " << tableOcc[g][alpha[(unsigned int)foundSymbol]] <<std::endl;
              }
              */
            }
          #else //We compute alternate order the BWT of the sequences
            if (tableOcc[currentPile][alpha[(unsigned int)foundSymbol]] > counters[alpha[(unsigned int)foundSymbol]])
              vectTriple[k].posN = tableOcc[currentPile][alpha[(unsigned int)foundSymbol]] - counters[alpha[(unsigned int)foundSymbol]]+1;
            else
              vectTriple[k].posN = 1;
  /*
            if (verboseEncode == 1)
              std::cerr << "\nInit New P["<< k <<"]= " << vectTriple[k].posN <<std::endl;
  */
            for (dataTypedimAlpha g = currentPile+1 ; g <sizeAlpha ; g++) {
                //I have to count in each pile g= (currentPile+1)...sizeAlpha-pile
              vectTriple[k].posN = vectTriple[k].posN + tableOcc[g][alpha[(unsigned int)foundSymbol]];
              //std::cerr << "--New posN[k]=" << (unsigned int)vectTriple[k].posN << " tableOcc[g][alpha[(unsigned int)symbol]] " << tableOcc[g][alpha[(unsigned int)symbol]] <<std::endl;
              /*
              if (verboseEncode == 1) {
                std::cerr << "g= " << (unsigned int)g << " symbol= " << (unsigned int)foundSymbol << " alpha[symbol]= "<< (unsigned int)alpha[(unsigned int)foundSymbol] <<std::endl;
                std::cerr << "Add New posN[k]=" << vectTriple[k].posN << " tableOcc[g][alpha[(unsigned int)symbol]] " << tableOcc[g][alpha[(unsigned int)foundSymbol]] <<std::endl;
              }
  */
            }

          #endif

          //I have to insert the new symbol in the symbol-pile

          int tmpOLDpile=vectTriple[k].pileN;

          vectTriple[k].pileN=alpha[(unsigned int)foundSymbol];

            if (vectTriple[k].pileN == (unsigned int)TERMINATE_CHAR) {
              std::cerr << "-->New posN[k]=" << (unsigned int)vectTriple[k].posN << " New pileN[k]=" << (unsigned int)vectTriple[k].pileN << " posSymb= " << posSymb << " tmpOLDpile= " << tmpOLDpile << std::endl;
              std::cerr << "*** " << " foundSymbol= " << (unsigned int)foundSymbol << " alpha[symbol]= "<< (unsigned int)alpha[(unsigned int)foundSymbol] <<std::endl;
            }

          //if (verboseEncode == 1)
          //  std::cerr << "j  : Q[q]=" << (unsigned int)vectTriple[k].pileN << " P[q]=" << (dataTypeNChar)vectTriple[k].posN <<  " N[q]=" << (dataTypeNSeq)vectTriple[k].seqN << std::endl;
         #if (BCR_SET==1)
        }  //end if (newSymb[vectTriple[k].seqN] != TERMINATE_CHAR_LEN)
         #endif
      k++;
      }

      #if BUILD_EXT_MEMORY==1
      //fclose(InFileBWT);
        closeFilePartial(InFileBWT);
      #endif

      j=k;
    }
    else
      j++;
  }  //end while (old sequence)

  //Insert the symbols beloging to new strings
  #if (BCR_SET==1)
    if (nExamedTexts < nText) {
      for (j=0; j < nText; j++) {
        //if ((vectInsTexts[j]==0) && (newSymb[j] == TERMINATE_CHAR_LEN))  //nothing to do
        //if ((vectInsTexts[j]==0) && (newSymb[j] != TERMINATE_CHAR_LEN))  //The new symbol have to be inserted in some y-pile with y>0

        if ((vectInsTexts[j]==0) && (newSymb[j] != TERMINATE_CHAR_LEN)) {   //The new symbol have to be inserted in 0-pile
          //std::cerr << "Insert new symbols!" << std::endl;
          //std::cerr << "+++j= " << j << " vectInsTexts[j]= " << vectInsTexts[j] << " newSymb[j]= " << newSymb[j]  << std::endl;
          dataTypeNSeq rankInserted=0;
          for (dataTypeNSeq x = 0 ; x < j; x++)
            rankInserted +=  vectInsTexts[x];
          sortElement tripla;
          tripla.posN = rankInserted+1;  // position of the suffix (1-based)
          tripla.seqN = j;    // number of the sequence
          tripla.pileN = 0;    //The first symbols are in $-pile
          #if BUILD_LCP == 1
            tripla.lcpCurN = 0;
            tripla.lcpSucN = 0;
          #endif
          vectTriple.push_back(tripla);
          vectInsTexts[j]=1;
          nExamedTexts++;

        }  //end if
      }  //end-for
    }  //end if new strings
  #endif

  delete [] counters;
  delete [] filenameIn;
  delete [] filename;

  quickSort(vectTriple);

  #if verboseEncode==1
    std::cerr << "U  ";
    for (dataTypeNSeq g = 0 ; g < nText; g++)
      if (vectInsTexts[g]==1)
        std::cerr << newSymb[g] << " ";
    std::cerr << std::endl;
    std::cerr << "After Sorting" << std::endl;
    std::cerr << "Q  ";
    for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
      std::cerr << (unsigned int)vectTriple[g].pileN << " ";
    }
    std::cerr << std::endl;
    std::cerr << "P  ";
    for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
      std::cerr << vectTriple[g].posN  << " ";
    }
    std::cerr << std::endl;
    std::cerr << "N  ";
    for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
      std::cerr << vectTriple[g].seqN  << " ";
    }
    std::cerr << std::endl;
    #if BUILD_LCP == 1
      std::cerr << "C  ";             //LCP current
      for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
        std::cerr << (unsigned int)vectTriple[g].lcpCurN  << " ";
      }
      std::cerr << std::endl;
      std::cerr << "S  ";                     //LCP successive
      for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
        std::cerr << (unsigned int)vectTriple[g].lcpSucN  << " ";
      }
      std::cerr << std::endl;
    #endif
  #endif


  #if BUILD_LCP == 0
    #if BUILD_EXT_MEMORY == 1
      storeBWTFilePartial(newSymb, posSymb);
    #else
      storeBWTIntMem(newSymb, posSymb);
    #endif
  #else
    storeBWTandLCP(newSymb, posSymb);  //Now it also compute the generalized suffix array (position and number of sequence)?
  #endif

  #if verboseEncode==1
    std::cerr << "tableOcc: after storeBWT or storeBWTandLCP" << std::endl;
    for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
      for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++)
        std::cerr << tableOcc[j][h] << " ";
      std::cerr << std::endl;
    }

    std::cerr << "Q  ";
    for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
      std::cerr << (unsigned int)vectTriple[g].pileN << " ";
    }
    std::cerr << std::endl;
    std::cerr << "P  ";
    for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
      std::cerr << vectTriple[g].posN  << " ";
    }
    std::cerr << std::endl;
    std::cerr << "N  ";
    for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
      std::cerr << vectTriple[g].seqN  << " ";
    }
    std::cerr << std::endl;
    #if BUILD_LCP == 1
      std::cerr << "C  ";
      for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
        std::cerr << vectTriple[g].lcpCurN  << " ";
      }
      std::cerr << std::endl;
      std::cerr << "S  ";
      for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
        std::cerr << vectTriple[g].lcpSucN  << " ";
      }
    #endif

    /* #if BUILD_SA == 1
      std::cerr << "\nPartial Generalized Suffix array: " << std::endl;
      dataTypeNChar numcharSA=0;
      static FILE *InFileSA;                  // output and input file SA;
      char *filenameInSA = new char[12];
      char *filenameSA = new char[8];
      ElementType *bufferGSA = new ElementType[SIZEBUFFER];
      mmm=0;
      while (mmm < sizeAlpha) {
        numcharSA=sprintf (filenameSA, "sa_%d", mmm);
        numcharSA=sprintf (filenameInSA,"%s%s",filenameSA,ext);
        InFileSA = fopen(filenameInSA, "rb");
        if (InFileSA==NULL) {
          std::cerr << "In SA file " << (unsigned int)mmm <<": Error opening: " << filenameInSA << std::endl;
          exit (EXIT_FAILURE);
        }
        for (dataTypeNChar g = 0 ; g < SIZEBUFFER; g++) {
          bufferGSA[g].sa = 0;
          bufferGSA[g].numSeq  = 0;
        }
        numcharSA = fread(bufferGSA,sizeof(ElementType),SIZEBUFFER,InFileSA);
        std::cerr << "SA[" << (unsigned int)mmm << "]:\t";
        if (numcharSA==0)
          std::cerr  << "empty";
        else {
          for (dataTypeNChar g = 0 ; g < numcharSA; g++)
            std::cerr << "(" << bufferGSA[g].sa  << " " << bufferGSA[g].numSeq << ")\t";
        }
        while (numcharSA!=0) {
          for (dataTypeNChar g = 0 ; g < SIZEBUFFER; g++) {
            bufferGSA[g].sa = 0;
            bufferGSA[g].numSeq  = 0;
          }
          numcharSA = fread(bufferGSA,sizeof(ElementType),SIZEBUFFER,InFileSA);
          if (numcharSA!=0)
            for (dataTypeNChar g = 0 ; g < numcharSA; g++)
              std::cerr << "(" << bufferGSA[g].sa  << " " << bufferGSA[g].numSeq << ")\t";
        }
        std::cerr << std::endl;
        fclose(InFileSA);
        mmm++;
      }
      delete [] filenameInSA;
      delete [] filenameSA;
      delete [] bufferGSA;
    #endif  */

  #endif


  /*
    if it is only useful for decoding
  */

  ////std::cerr << "End storing BWT" << std::endl;
  //if (posSymb == lengthRead) { //We are storing the last simbols of sequences: $
  //  std::cerr << "Stores the 'end positions' of the $!"<< std::endl;
  //  char *fileEndPos="outFileEndPos.bwt";
  //  static FILE *OutFileEndPos;                  // output file of the end positions;
  //  OutFileEndPos = fopen(fileEndPos, "wb");
  //  if (OutFileEndPos==NULL) {
  //      std::cerr << "Error opening \"" << fileEndPos << "\" file"<< std::endl;
  //      exit (EXIT_FAILURE);
  //  }

  //  /*
  //  //Each symbol newSymb[seqN[i]] has been inserted in position posN[i] into the pile pileN[i]
  //  //We have to store the absolute positions in the entire BWT
  //  //So we need to use the tableOcc.
  //  //The symbol $ of the sequence i is in the position endPos[SeqN[i]]
  //  for (dataTypeNSeq i = 0; i < nText; i++) {
  //    for (dataTypedimAlpha r = 0; r < vectTriple[i].pileN; r++) {
  //      for (dataTypedimAlpha t = 0; t < sizeAlpha; t++) {
  //        vectTriple[i].posN += tableOcc[r][t];
  //      }
  //    }
  //  }
  //
  //  std::cerr << "Positions of the EOF into BWT" << std::endl;
  //  for (dataTypeNSeq i = 0; i < nText; i++) {
  //    std::cerr << posN[i] << " ";
  //  }
  //  std::cerr << std::endl;
  //  */
  //
  //  numchar = fwrite (&nText, sizeof(dataTypeNChar), 1 , OutFileEndPos);
  //  assert( numchar == 1); // we should always read the same number of characters

  //  for (dataTypeNSeq i = 0; i < nText; i++) {
  //    if (verboseEncode == 1)
  //      std::cerr << "Triple: " << vectTriple[i].seqN << " " << vectTriple[i].posN << " " << (unsigned int)vectTriple[i].pileN << std::endl;
  //    numchar = fwrite (&vectTriple[i].seqN, sizeof(dataTypeNSeq), 1 , OutFileEndPos);
  //    assert( numchar == 1); // we should always read the same number of characters
  //    numchar = fwrite (&vectTriple[i].posN, sizeof(dataTypeNChar), 1 , OutFileEndPos); //here vectTriple[i].posN is the relative position of $ in the partial BWT
  //    assert( numchar == 1); // we should always read the same number of characters
  //    numchar = fwrite (&vectTriple[i].pileN, sizeof(dataTypedimAlpha), 1 , OutFileEndPos);
  //    assert( numchar == 1); // we should always read the same number of characters
  //  }

  //  fclose(OutFileEndPos);
  //  std::cerr << "The'end positions' stored!"<< std::endl;
  //}
}

void BCRexternalBWT::storeBWTFilePartial(uchar const * newSymb, dataTypelenSeq posSymb) {

  //I have found the position where I have to insert the chars in the position t of the each text
  //Now I have to update the BWT in each file.
  const char *ext = ".aux";

  #if BUILD_EXT_MEMORY==1
    static FILE *OutFileBWT, *InFileBWT;                  // output and input file BWT;
    char *filenameOut = new char[16];
    char *filenameIn = new char[12];
    char *filename = new char[8];
    dataTypeNChar numcharWrite=0;
  #endif



  #if BUILD_SA == 1
    static FILE *OutFileSA, *InFileSA;                  // output and input file SA;
    char *filenameOutSA = new char[16];
    char *filenameInSA = new char[12];
    char *filenameSA = new char[8];
    ElementType *bufferSA = new ElementType[SIZEBUFFER];
    dataTypeNChar numcharSA=0;
    dataTypeNChar numcharWriteSA=0;
  #endif

  dataTypeNChar numchar=0;
  uchar *buffer = new uchar[SIZEBUFFER];
  dataTypeNChar toRead = 0;

  dataTypeNSeq j;
  dataTypedimAlpha currentPile=vectTriple[0].pileN;

  //std::cerr << "storeBWTFilePartial: posSymb= " << posSymb << std::endl;

  //uchar symbol='\0';
  j=0;
  while (j < nExamedTexts) {

    currentPile = vectTriple[j].pileN;

    #if BUILD_EXT_MEMORY==1
      numchar=sprintf (filename, "bwt_%d", currentPile);

      /*
      numchar=sprintf (filenameIn,"%s%s",filename,ext);
      InFileBWT = fopen(filenameIn, "rb");
      if (InFileBWT==NULL) {
        std::cerr << "storeBWTFilePartial: In BWT file, j= " << (unsigned int)j <<": Error opening " << std::endl;
        exit (EXIT_FAILURE);
      } */
      InFileBWT = openFilePartialIn (currentPile);
      numchar=sprintf (filenameOut,"new_%s%s",filename,ext);
      OutFileBWT = fopen(filenameOut, "wb");
      if (OutFileBWT==NULL) {
          std::cerr << "storeBWTFilePartial: Out BWT file, j= " << (unsigned int)j <<": Error opening " << std::endl;
          exit (EXIT_FAILURE);
      }
      //OutFileBWT = openFilePartialOut  (currentPile);
    #endif

    #if BUILD_SA == 1
      numcharSA=sprintf (filenameSA, "sa_%d", currentPile);
      numcharSA=sprintf (filenameInSA,"%s%s",filenameSA,ext);
      InFileSA = fopen(filenameInSA, "rb");
      if (InFileSA==NULL) {
        std::cerr << "storeBWTFilePartial: In SA file, j= " << (unsigned int)j <<": Error opening: " << filenameInSA << std::endl;
        exit (EXIT_FAILURE);
      }
      numcharSA=sprintf (filenameOutSA,"new_%s%s",filenameSA,ext);
      OutFileSA = fopen(filenameOutSA, "wb");
      if (OutFileSA==NULL) {
        std::cerr << "storeBWTFilePartial: Out SA file, j= " << (unsigned int)j <<": Error opening: " << filenameOutSA << std::endl;
        exit (EXIT_FAILURE);
      }
    #endif

    //For each new symbol in the same pile
    dataTypeNSeq k=j;
    dataTypeNChar cont = 0;
    while ((k< nExamedTexts) && (vectTriple[k].pileN == currentPile)) {

       #if (BCR_SET==1)
      if (newSymb[vectTriple[k].seqN] != TERMINATE_CHAR_LEN) {
         #endif
        //if (verboseEncode==1)
         //   std::cerr << "k= " << k << " Q[k]= " << (unsigned int)vectTriple[k].pileN << " P[k]= " << vectTriple[k].posN << " cont = "<< cont << std::endl;
        //So I have to read the k-BWT and I have to count the number of the symbols up to the position posN.
        //symbol = '\0';
        //As PosN starts to the position 1 and I have to insert the new symbol in position posN[k]
        // I have to read posN[k]-1 symbols
        //cont is the number of symbols already read!
        toRead = (vectTriple[k].posN-1) - cont;
        while (toRead > 0) {            //((numchar!=0) && (toRead > 0)) {
          if (toRead < SIZEBUFFER) { //The last reading for this sequence
            //numchar = fread(buffer,sizeof(uchar),toRead,InFileBWT);
            #if BUILD_EXT_MEMORY==1
              numchar =  readOnFilePartial(buffer, toRead, InFileBWT) ;
              assert(numchar == toRead); // we should always read/write the same number of characters
              numcharWrite =  writeOnFilePartial(buffer, numchar, OutFileBWT) ;
              //numcharWrite = fwrite (buffer, sizeof(uchar), numchar , OutFileBWT);
              assert(numchar == numcharWrite); // we should always read/write the same number of characters
            #endif

            #if BUILD_SA == 1
                numcharSA = fread(bufferSA,sizeof(ElementType),toRead,InFileSA);
                assert(numcharSA == toRead);
                numcharWriteSA = fwrite (bufferSA, sizeof(ElementType), numcharSA , OutFileSA);
                assert(numcharSA == numcharWriteSA);
            #endif
          }
          else {
            #if BUILD_EXT_MEMORY==1
              numchar =  readOnFilePartial(buffer, SIZEBUFFER, InFileBWT) ;
              //numchar = fread(buffer,sizeof(uchar),SIZEBUFFER,InFileBWT);
              assert(numchar == SIZEBUFFER); // we should always read/write the same number of characters
              numcharWrite =  writeOnFilePartial(buffer, numchar, OutFileBWT) ;
              //numcharWrite = fwrite (buffer, sizeof(uchar), numchar , OutFileBWT);
              assert(numchar == numcharWrite); // we should always read/write the same number of characters
            #endif

            #if BUILD_SA == 1
                numcharSA = fread(bufferSA,sizeof(ElementType),SIZEBUFFER,InFileSA);
                assert(numcharSA == SIZEBUFFER); // we should always read/write the same number of characters
                numcharWriteSA = fwrite (bufferSA, sizeof(ElementType), numcharSA , OutFileSA);
                assert(numcharSA == numcharWriteSA); // we should always read/write the same number of characters
            #endif
          }

          cont   += numchar;  //number of read symbols
          toRead -= numchar;
          if ((numchar == 0) && (toRead > 0)) {  //it means that we have read 0 character, but there are still toRead characters to read
            std::cerr << "storeBWTFilePartial: sequence number" << (unsigned int)k << " read 0 character, but there are still " << toRead << " characters to read  " << std::endl;
            exit (EXIT_FAILURE);
          }

        }
        //Now I have to insert the new symbol associated with the suffix of the sequence k
        //And I have to update the number of occurrences of each symbol
        if (toRead==0) {
          //numchar = fwrite (&newSymb[vectTriple[k].seqN], sizeof(uchar), 1, OutFileBWT);
          #if BUILD_EXT_MEMORY==1
            numchar =  writeSymbolOnFilePartial(newSymb[vectTriple[k].seqN], 1, OutFileBWT) ;
            assert(numchar == 1); // we should always read/write the same number of characters
          #endif

          tableOcc[currentPile][alpha[(unsigned int)newSymb[vectTriple[k].seqN]]]++;       //update the number of occurrences in BWT of the pileN[k]


          #if BUILD_SA == 1
            ElementType newEle;
            newEle.sa= posSymb;
            newEle.numSeq=vectTriple[k].seqN;
            numcharSA = fwrite (&newEle, sizeof(ElementType), 1, OutFileSA);
            assert(numcharSA == 1);
          #endif
          cont++;    //number of read symbols
          toRead--;
        }
       #if (BCR_SET==1)
      }   //end if (newSymb[vectTriple[k].seqN] != TERMINATE_CHAR_LEN) {
       #endif
      k++;   //  I changed the number of the sequence. New iteration.
    }
    //it means that posN[k]<>currentPile, so I have to change BWT-file
    //But before, I have to copy the remainder symbols from the old BWT to new BWT
    while (numchar!=0) {

      #if BUILD_EXT_MEMORY==1
        numchar =  readOnFilePartial(buffer, SIZEBUFFER, InFileBWT) ;
        //numchar = fread(buffer,sizeof(uchar),SIZEBUFFER,InFileBWT);
        if (numchar > 0) {
          numcharWrite =  writeOnFilePartial(buffer, numchar, OutFileBWT) ;
          //numcharWrite = fwrite (buffer, sizeof(uchar), numchar , OutFileBWT);
          assert(numchar == numcharWrite); // we should always read/write the same number of characters
        }
      #endif


      #if BUILD_SA == 1
        numcharSA = fread(bufferSA,sizeof(ElementType),SIZEBUFFER,InFileSA);
        numcharWriteSA = fwrite (bufferSA, sizeof(ElementType), numcharSA , OutFileSA);
        assert(numcharSA == numcharWriteSA);
      #endif
    }

    #if BUILD_EXT_MEMORY==1
      /* fclose(InFileBWT);
      fclose(OutFileBWT); */
      closeFilePartial(InFileBWT);
      closeFilePartial(OutFileBWT);
    #endif


    #if BUILD_SA == 1
      fclose(InFileSA);
      fclose(OutFileSA);
    #endif
    j=k;
  }

  //static FILE *tmpFile;
  //tmpFile = fopen("sizeFileBwt.txt", "a");

  #if BUILD_EXT_MEMORY==1
    //Renaming new to old (bwt files)
    for (dataTypedimAlpha g = 0 ; g < sizeAlpha; g++) {
      /* numchar=sprintf (filename, "bwt_%d", g);
      numchar=sprintf (filenameIn,"%s%s",filename,ext);
      numchar=sprintf (filenameOut,"new_%s%s",filename,ext);
      //std::cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << std::endl;
      OutFileBWT = fopen(filenameOut, "rb");
      if (OutFileBWT!=NULL) { //If it exists
        fclose(OutFileBWT);
        if (remove(filenameIn)!=0)
          std::cerr << filenameIn <<": Error deleting file" << std::endl;
        else
          if(rename(filenameOut,filenameIn))
            std::cerr << filenameOut <<": Error renaming " << std::endl;
      } */
      renameFilePartial(g);



  /*    #if verboseEncode==1
        struct stat results;
        if (stat(filenameIn, &results) == 0)
          // The size of the file in bytes is in results.st_size
          //fprintf(tmpFile,"%s\t%u\n", filenameIn, results.st_size);
          std::cerr << filenameIn <<"\t" << results.st_size << std::endl;
        else
          //fprintf(tmpFile,"An error occurred %s\n", filenameIn);
          std::cerr << "An error occurred" << std::endl;
      #endif
  */
    }

    #if BUILD_SA == 1
        //Renaming new to old
        for (dataTypedimAlpha g = 0 ; g < sizeAlpha; g++) {
          numchar=sprintf (filenameSA, "sa_%d", g);
          numchar=sprintf (filenameInSA,"%s%s",filenameSA,ext);
          numchar=sprintf (filenameOutSA,"new_%s%s",filenameSA,ext);
          //std::cerr << "Filenames:" << filenameInSA << "\t" <<filenameOutSA << std::endl;
          OutFileSA = fopen(filenameOutSA, "rb");

          if (OutFileSA!=NULL) { //If it exists
            fclose(OutFileSA);
            if (remove(filenameInSA)!=0)
              std::cerr << filenameInSA <<": Error deleting SA file" << std::endl;
            else
              if(rename(filenameOutSA,filenameInSA))
                std::cerr << filenameOutSA <<": Error renaming SA file " << std::endl;
          }


          if (verboseEncode == 1) {
            struct stat results;
            if (stat(filenameInSA, &results) == 0)
              // The size of the file in bytes is in results.st_size
              //fprintf(tmpFile,"%s\t%u\n", filenameInSA, results.st_size);
              std::cerr << filenameInSA <<"\t" << results.st_size << std::endl;
            else
              //fprintf(tmpFile,"An error occurred %s\n", filenameInSA);
              std::cerr << "An error occurred" << std::endl;
          }
        }
    #endif

    #if BUILD_EXT_MEMORY==1
      delete [] filenameIn;
      delete [] filename;
      delete [] filenameOut;
    #endif


    #if BUILD_SA == 1
      delete [] filenameInSA;
      delete [] filenameSA;
      delete [] filenameOutSA;
    #endif

  #endif

  delete [] buffer;


  #if BUILD_SA == 1
    delete [] bufferSA;
  #endif
}

#if BUILD_EXT_MEMORY == 0
void BCRexternalBWT::storeBWTIntMem(uchar const * newSymb, dataTypelenSeq posSymb) {
  //I have found the position where I have to insert the chars in the position t of the each text
  //Now I have to update the partial BWT
  dataTypeNSeq j=0;
  dataTypedimAlpha currentPile=vectTriple[0].pileN;
  std::vector<char> vectBWTcurrentPile;
  dataTypeNChar toRead = 0;

  while (j < nExamedTexts) {
    currentPile = vectTriple[j].pileN;
    vectBWTcurrentPile.clear();

    //For each new symbol in the same pile
    dataTypeNSeq k=j;
    dataTypeNChar cont = 0;
    dataTypeNChar eleCurrentPile=0;
    while ((k< nExamedTexts) && (vectTriple[k].pileN == currentPile)) {
       #if (BCR_SET==1)
      if (newSymb[vectTriple[k].seqN] != TERMINATE_CHAR_LEN) {
       #endif
        //As PosN starts to the position 1 and I have to insert the new symbol in position posN[k]
        // I have to read posN[k]-1 symbols
        //cont is the number of symbols already read!
        toRead = (vectTriple[k].posN-1) - cont;
        while (toRead > 0) {
          vectBWTcurrentPile.push_back ( vectVectBWT[currentPile][eleCurrentPile] );
          toRead--;
          cont++;    //number of read symbols
          eleCurrentPile++;
        }
        //Now I have to insert the new symbol associated with the suffix of the sequence k
        //And I have to update the number of occurrences of each symbol
        if (toRead==0) {
          vectBWTcurrentPile.push_back ( newSymb[vectTriple[k].seqN] );
          cont++;    //number of read symbols
          toRead--;
          tableOcc[currentPile][alpha[(unsigned int)newSymb[vectTriple[k].seqN]]]++;       //update the number of occurrences in BWT of the pileN[k]
        }
        else {
          std::cerr << "storeBWTIntMem: Error toRead>0, i.e. toRead= " << toRead <<  std::endl;
          exit (EXIT_FAILURE);
        }
       #if (BCR_SET==1)
      }
       #endif
      k++;   //  I changed the number of the sequence. New iteration.
    }



    //it means that posN[k]<>currentPile, so I have to change BWT-segment
    //But before, I have to copy the remainder symbols from the old BWT to new BWT
    while (eleCurrentPile < vectVectBWT[currentPile].size() ) {
      vectBWTcurrentPile.push_back ( vectVectBWT[currentPile][eleCurrentPile] );
      toRead--;
      cont++;    //number of read symbols
      eleCurrentPile++;
    }

    //Now I have to replace vectVectBWT[currentPile] with vectBWTcurrentPile
    vectVectBWT[currentPile].clear();
    vectVectBWT[currentPile] = vectBWTcurrentPile;
    j=k;
  }
}
#endif


#if BUILD_LCP == 1
void BCRexternalBWT::storeBWTandLCP(uchar const * newSymb, dataTypelenSeq posSymb) {
  //if (verboseEncode==1)
  //  std::cerr << "Store nText symbols in the old bwts" << std::endl;
  dataTypelenSeq maxValueLen = lengthRead+2;
  vector <dataTypelenSeq> minLCPcur;
  vector <bool> minLCPcurFound;
  vector <dataTypelenSeq> minLCPsuc;
  vector <dataTypeNSeq> minLCPsucText;
    vector <bool> minLCPsucToFind;
  minLCPcur.resize(sizeAlpha);       //for each symbol of the alphabet
  minLCPcurFound.resize(sizeAlpha);
  minLCPsuc.resize(sizeAlpha);
  minLCPsucText.resize(sizeAlpha);
    minLCPsucToFind.resize(sizeAlpha);

  //I have found the position where I have to insert the chars in the position t of the each text
  //Now I have to update the BWT in each file.
  static FILE *OutFileBWT, *InFileBWT;                  // output and input file BWT;
  char *filenameOut = new char[16];
  char *filenameIn = new char[12];
  char *filename = new char[8];
  static FILE *OutFileLCP, *InFileLCP;                  // output and input file LCP;
  char *filenameOutLCP = new char[16];
  char *filenameInLCP = new char[12];
  char *filenameLCP = new char[8];
  const char *ext = ".aux";



  #if BUILD_SA == 1
    static FILE *OutFileSA, *InFileSA;                  // output and input file SA;
    char *filenameOutSA = new char[16];
    char *filenameInSA = new char[12];
    char *filenameSA = new char[8];
    ElementType *bufferSA = new ElementType[SIZEBUFFER];
    dataTypeNChar numcharSA=0;
    dataTypeNChar numcharWriteSA=0;
  #endif

  dataTypeNChar numchar=0;
  dataTypeNChar numcharWrite=0;
  uchar *buffer = new uchar[SIZEBUFFER];
  dataTypelenSeq *bufferLCP = new dataTypelenSeq[SIZEBUFFER];
  dataTypeNChar toRead = 0;

  dataTypeNSeq j;
  dataTypedimAlpha currentPile;
  //uchar symbol='\0';

  j=0;
  while (j < nExamedTexts) {
    currentPile = vectTriple[j].pileN;
    //if (verboseEncode==1)
      //std::cerr << "\n---storeBWTandLCP: New Segment; index text j= " << j << " current BWT segment is " << (unsigned int)currentPile << std::endl;
    for (dataTypedimAlpha g = 0 ; g < sizeAlpha; g++) {
      minLCPcur[g]=maxValueLen;
      minLCPcurFound[g]=0;
      minLCPsuc[g]=maxValueLen;
      minLCPsucText[g]=0;     //denotes the number of the text associated with the symbol g
            minLCPsucToFind[g] = 0;
    }
    //std::cerr << "---Pile " << (unsigned int)currentPile << std::endl;
        assert(currentPile <= sizeAlpha);
    numchar=sprintf (filename, "bwt_%d", currentPile);
    numchar=sprintf (filenameIn,"%s%s",filename,ext);
    InFileBWT = fopen(filenameIn, "rb");
    if (InFileBWT==NULL) {
      std::cerr << "In BWT file " << (unsigned int)j <<": Error opening: " << filenameIn << std::endl;
      exit (EXIT_FAILURE);
    }
    numchar=sprintf (filenameOut,"new_%s%s",filename,ext);
    OutFileBWT = fopen(filenameOut, "wb");
    if (OutFileBWT==NULL) {
        std::cerr << "Out BWT file " << (unsigned int)j <<": Error opening: " << filenameOut << std::endl;
        exit (EXIT_FAILURE);
    }
    //std::cerr << "In File " << filenameIn << std::endl;
    //std::cerr << "Out File " << filenameOut << std::endl;
    numchar=sprintf (filenameLCP, "lcp_%d", currentPile);
    numchar=sprintf (filenameInLCP,"%s%s",filenameLCP,ext);
    InFileLCP = fopen(filenameInLCP, "rb");
    if (InFileLCP==NULL) {
      std::cerr << "In LCP file " << (unsigned int)j <<": Error opening: " << filenameInLCP << std::endl;
      exit (EXIT_FAILURE);
    }
    numchar=sprintf (filenameOutLCP,"new_%s%s",filenameLCP,ext);
    OutFileLCP = fopen(filenameOutLCP, "wb");
    if (OutFileLCP==NULL) {
        std::cerr << "Out LCP file " << (unsigned int)j <<": Error opening: " << filenameOutLCP << std::endl;
        exit (EXIT_FAILURE);
    }



    #if BUILD_SA == 1
      numcharSA=sprintf (filenameSA, "sa_%d", currentPile);
      numcharSA=sprintf (filenameInSA,"%s%s",filenameSA,ext);
      InFileSA = fopen(filenameInSA, "rb");
      if (InFileSA==NULL) {
        std::cerr << "storeBWTandLCP: In SA file " << (unsigned int)j <<": Error opening: " << filenameInSA << std::endl;
        exit (EXIT_FAILURE);
      }
      numcharSA=sprintf (filenameOutSA,"new_%s%s",filenameSA,ext);
      OutFileSA = fopen(filenameOutSA, "wb");
      if (OutFileSA==NULL) {
        std::cerr << "storeBWTandLCP: Out SA file " << (unsigned int)j <<": Error opening: " << filenameOutSA << std::endl;
        exit (EXIT_FAILURE);
      }
      //std::cerr << "SA In File " << filenameInSA << std::endl;
      //std::cerr << "SA Out File " << filenameOutSA << std::endl;
    #endif


    //For each new symbol in the same pile
    dataTypeNSeq k=j;
    dataTypeNChar cont = 0;
    while ((k< nExamedTexts) && (vectTriple[k].pileN == currentPile)) {
      //std::cerr << "Start k= "<< k <<": Symb="<< newSymb[vectTriple[k].seqN] << ", Seq=" << vectTriple[k].seqN << ", Pos=" << vectTriple[k].posN << ", cont="<< cont << "\n";
      #if (BCR_SET==1)
            if (newSymb[vectTriple[k].seqN] != TERMINATE_CHAR_LEN) {
      #endif
                //if (verboseEncode==1)
             //std::cerr << "k= " << k << " Q[k]= " << (unsigned int)vectTriple[k].pileN << " P[k]= " << vectTriple[k].posN << " cont = "<< cont << std::endl;
        //So I have to read the k-BWT and I have to count the number of the symbols up to the position posN.
        //symbol = '\0';
        //As PosN starts to the position 1 and I have to insert the new symbol in position posN[k]
        // I have to read posN[k]-1 symbols
        //cont is the number of symbols already read!
        toRead = (vectTriple[k].posN-1) - cont;
        //if (verboseEncode == 1)
        //std::cerr << "Before of (toRead > 0) - Start: symb="<< newSymb[vectTriple[k].seqN] <<", vectTriple["<< k << "].posN=" << vectTriple[k].posN << ", cont="<< cont <<", toRead= " << toRead << "\n";

        while (toRead > 0) {            //((numchar!=0) && (toRead > 0)) {
          if (toRead < SIZEBUFFER) { //The last reading for this sequence
            numchar = fread(buffer,sizeof(uchar),toRead,InFileBWT);
            //if (verboseEncode == 1)
            //std::cerr << "BWT: number read " << numchar << " to Read " << toRead << "\n";
            assert(numchar == toRead); // we should always read/write the same number of characters
            numcharWrite = fwrite (buffer, sizeof(uchar), numchar , OutFileBWT);
            assert(numchar == numcharWrite); // we should always read/write the same number of characters
            //std::cerr << "toread number write " << numcharWrite << "\n";
            numchar = fread(bufferLCP,sizeof(dataTypelenSeq),toRead,InFileLCP);
            //if (verboseEncode == 1)
              //std::cerr << "In LCP: number read " << numchar << " to Read " << toRead << "\n";
            assert(numchar == toRead); // we should always read/write the same number of characters
            numcharWrite = fwrite (bufferLCP,sizeof(dataTypelenSeq), numchar , OutFileLCP);
            //if (numchar != numcharWrite)   //***********************************
              //std::cerr << "In LCP: numchar " << numchar << " numchar " << numcharWrite << " to Read " << toRead << "\n";
            assert(numchar == numcharWrite); // we should always read/write the same number of characters

            #if BUILD_SA == 1
              numcharSA = fread(bufferSA,sizeof(ElementType),toRead,InFileSA);
              //std::cerr << "In SA: numcharSA " << numcharSA << " to Read " << toRead << "\n";
              assert(numcharSA == toRead);
              //std::cerr << "In SA: numcharSA " << numcharSA << " numcharWriteSA " << numcharWriteSA << " to Read " << toRead << "\n";
              numcharWriteSA = fwrite (bufferSA, sizeof(ElementType), numcharSA , OutFileSA);
              assert(numcharSA == numcharWriteSA);
            #endif
          }
          else {
            numchar = fread(buffer,sizeof(uchar),SIZEBUFFER,InFileBWT);
            //if (verboseEncode == 1)
            //std::cerr << "number read " << numchar << "\n";
            assert(numchar == SIZEBUFFER); // we should always read/write the same number of characters
            numcharWrite = fwrite (buffer, sizeof(uchar), numchar , OutFileBWT);
            assert(numchar == numcharWrite); // we should always read/write the same number of characters
            //std::cerr << "sizebuffer number write " << numcharWrite << "\n";
            numchar = fread(bufferLCP,sizeof(dataTypelenSeq),SIZEBUFFER,InFileLCP);
            assert(numchar == SIZEBUFFER); // we should always read/write the same number of characters
            numcharWrite = fwrite (bufferLCP,sizeof(dataTypelenSeq), numchar , OutFileLCP);
            assert(numchar == numcharWrite); // we should always read/write the same number of characters.

            #if BUILD_SA == 1
              numcharSA = fread(bufferSA,sizeof(ElementType),SIZEBUFFER,InFileSA);
              assert(numcharSA == SIZEBUFFER); // we should always read/write the same number of characters
              numcharWriteSA = fwrite (bufferSA, sizeof(ElementType), numcharSA , OutFileSA);
              assert(numcharSA == numcharWriteSA); // we should always read/write the same number of characters
            #endif
          }
          //I must to compute the minimum LCP. It is need to compute the lcpValue for the next iteration
          //std::cerr << "For each letter in the buffer before of the position where I have to insert the new symbol\n";
          for (dataTypeNChar bb = 0 ; bb < numcharWrite; bb++) {
            //Update the min1 for each letter of the alphabet, for which I have already met the symbol
            for (dataTypedimAlpha gg = 0 ; gg < sizeAlpha; gg++) {
              if (minLCPcurFound[gg]==1)   //I already met the symbol gg. So, I must compute the minimum
                if (bufferLCP[bb] < minLCPcur[gg])  //comparison with the last inserted lcp
                  minLCPcur[gg] = bufferLCP[bb];
            }

            minLCPcur[alpha[(unsigned int)buffer[bb]]] = maxValueLen; //For each occurrence of buffer[bb], I have to set the min1 (the interval starts from the next symbol)
            minLCPcurFound[alpha[(unsigned int)buffer[bb]]] = 1;  //So I open the LCP interval for buffer[bb] (for the first occurrence of buffer[bb])
            //std::cerr << "minLCPcur[" << buffer[bb]<<  "]= " << minLCPcur[alpha[(unsigned int)buffer[bb]]]  << " minLCPcurFound[" << buffer[bb]<<  "]= " << minLCPcurFound[alpha[(unsigned int)buffer[bb]]] << "\n";

            //First, I check if the symbol buffer[bb] closes a previous interval or it is in the middle or no.
            //In any case, for each symbol, we have to update the minimum
            for (dataTypedimAlpha gg = 0 ; gg < sizeAlpha; gg++) {
                  if (minLCPsucToFind[(unsigned int)gg] == 1) {    //We have to compute the minimum for the lcp of the symbol gg
                    //std::cerr << "---1) (dovrebbe essere <MAX) minLCPsuc[" << (unsigned int)gg <<  "]= " << minLCPsuc[(unsigned int)gg] << "\n";
                  if (bufferLCP[bb] < minLCPsuc[(unsigned int)gg])  //comparison with the last inserted lcp
                  minLCPsuc[(unsigned int)gg] = bufferLCP[bb];
                     }
            }
            //If the symbol buffer[bb] does not close an LCP interval, ok!
            if (minLCPsucToFind[alpha[(unsigned int)buffer[bb]]] == 1) {   //The symbol buffer[bb] closed an interval LCP
              //We have to compute the minimum for the lcp of the symbol buffer[bb]
                //std::cerr << "1a) (dovrebbe essere < MAX) minLCPsuc[" << (unsigned int)buffer[bb] <<  "]= " << minLCPsuc[alpha[(unsigned int)buffer[bb]]] << "\n";
              //if (minLCPsuc[alpha[(unsigned int)buffer[bb]]] != maxValueLen) {       //it means that there is an interval lcp to close for this symbol
              //The symbol buffer[bb] closes a LCP interval
              //I have already computed the minimum (close the previous lcp interval).
              //I can set lcpSucN of minLCPsucText[alpha[(unsigned int)[buffer[bb]]]
              vectTriple[minLCPsucText[alpha[(unsigned int)buffer[bb]]]].lcpSucN  = minLCPsuc[alpha[(unsigned int)buffer[bb]]];
              //std::cerr << "---Set: lcpSucN[" << (unsigned int)buffer[bb] <<  "]= " << minLCPsuc[alpha[(unsigned int)buffer[bb]]];
              //It closes the LCP interval, so
              minLCPsuc[alpha[(unsigned int)buffer[bb]]] = maxValueLen;
              minLCPsucToFind[alpha[(unsigned int)buffer[bb]]] = 0;
              minLCPsucText[alpha[(unsigned int)buffer[bb]]] = 0;
              //std::cerr << "---minLCPsuc[" << buffer[bb]<<  "]= " << minLCPsuc[alpha[(unsigned int)buffer[bb]]]  << " minLCPsucText[" << buffer[bb]<<  "]= " << minLCPsucText[alpha[(unsigned int)buffer[bb]]] << "\n";
              //}
            }
          }
          cont   += numchar;  //number of read symbols
          toRead -= numchar;
          if ((numchar == 0) && (toRead > 0)) {  //it means that we have read 0 character, but there are still toRead characters to read
            std::cerr << "storeBWTandLCP: sequence number" << (unsigned int)k << " read 0 character, but there are still " << toRead << " characters to read  " << std::endl;
            exit (EXIT_FAILURE);
          }
        }
        //Now I have to insert the new symbol associated with the suffix of the sequence k
        //And I have to update the number of occurrences of each symbol
        //And I have to insert the valueLCP store in lcpCurN + 1 in the previous iteration
        if (toRead==0) {
          //std::cerr << "\n---Now I can insert the new symbol and lcp, indeed toRead= " << toRead << std::endl;
          //std::cerr << "---STO inserendo newSymb[vectTriple[k].seqN]= " << newSymb[vectTriple[k].seqN] << " ";
          //std::cerr << "---per la sequenza  " << vectTriple[k].seqN << " in the segment " << currentPile << "\n";
          numchar = fwrite (&newSymb[vectTriple[k].seqN], sizeof(uchar), 1, OutFileBWT);
          assert(numchar == 1); // we should always read/write the same number of characters
          tableOcc[currentPile][alpha[(unsigned int)newSymb[vectTriple[k].seqN]]]++;       //update the number of occurrences in BWT of the pileN[k]
          //std::cerr << "new number write " << numchar << "\n";

          #if BUILD_SA == 1
            ElementType newEle;
            newEle.sa= posSymb;
            newEle.numSeq=vectTriple[k].seqN;
            numcharSA = fwrite (&newEle, sizeof(ElementType), 1, OutFileSA);
            assert(numcharSA == 1);
          #endif

          dataTypelenSeq lcpValueNow;
          if (vectTriple[k].pileN == 0) { //In 0-pile the value LCP is 0
            lcpValueNow = vectTriple[k].lcpCurN;
          }
          else {
            if (vectTriple[k].posN == 1) {   //it is the first symbol of the segment
              lcpValueNow = vectTriple[k].lcpCurN;
              //std::cerr << "\nFirst position of the segment so lcp should be 0, indeed is = " << vectTriple[k].lcpCurN << std::endl;
            }
            else
              lcpValueNow = vectTriple[k].lcpCurN + 1;
          }

          numchar = fwrite (&lcpValueNow, sizeof(dataTypelenSeq), 1, OutFileLCP);   //Insert the lcp for the new symbol
          assert(numchar == 1); // we should always read/write the same number of characters
          //std::cerr << "I insert the symbol= " << newSymb[vectTriple[k].seqN] <<  " and lcp " << lcpValueNow << std::endl;
          //Update the lcpCurN for the next iteration
          if (minLCPcurFound[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] == 0) {
            if (minLCPcur[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] == maxValueLen) {
              //it means that we have not met the symbol before of the position posN because minLCPcurFound is 0.
              //if minLCPcurFound is 0, then minLCPcur is maxValueLen
              vectTriple[k].lcpCurN = 0;         //The next suffix has suffix 0+1=1
              //std::cerr << "Found is 0 --> vectTriple["<<k<<"].lcpCurN = " << vectTriple[k].lcpCurN  << std::endl;
            }
          }
          else {  //it means that (minLCPcurFound[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] == 1)
            if (minLCPcur[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] == maxValueLen) {
              //it means that we have met the symbol before of the position posN because minLCPcurFound is 1.
              //But minLCPcur is maxValueLen, this means that the previous occurrences of new symbol is the previous position
              vectTriple[k].lcpCurN = lcpValueNow;         //The next suffix has suffix lcpValueNow+1
              //std::cerr << "Found is 1 and minLCPcur is max --> vectTriple["<<k<<"].lcpCurN = " << vectTriple[k].lcpCurN << "\n";
            }
            else {
              if (lcpValueNow < minLCPcur[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]]) {
                //comparison with the last inserted lcp. It means that the previous occurrence of new symbol is not the previous symbol
                vectTriple[k].lcpCurN = lcpValueNow;
                //std::cerr << "Found is 1 and lcpValueNow is lower --> vectTriple["<<k<<"].lcpCurN = " << vectTriple[k].lcpCurN << "\n";
              }
              else {
                vectTriple[k].lcpCurN  = minLCPcur[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]];
                //std::cerr << "Found is 1 and minimum is lower --> vectTriple["<<k<<"].lcpCurN = " << vectTriple[k].lcpCurN << "\n";
              }
            }
          }

          //std::cerr << "*We have computed minLCPcur[" << newSymb[vectTriple[k].seqN] <<  "] is " << minLCPcur[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] << std::endl;
          //std::cerr << "*vectTriple[k].lcpCurN is " << vectTriple[k].lcpCurN << std::endl;

          //it may happen that the new symbol closes a previous interval or it is in the middle or no.
          //In any case, for each symbol, we have to update the minimum
          for (dataTypedimAlpha gg = 0 ; gg < sizeAlpha; gg++) {
            if (minLCPsucToFind[(unsigned int)gg] == 1) {    //We have to compute the minimum for the lcp of the symbol gg
              //std::cerr << "2) (dovrebbe essere <MAX) minLCPsuc[" << (unsigned int)gg <<  "]= " << minLCPsuc[(unsigned int)gg]  << "\n";
                        //std::cerr << "     minLCPsucText[" << (unsigned int)gg <<  "]= " << minLCPsucText[(unsigned int)gg]  << "\n";
                      //std::cerr << "     minLCPsucToFind[" << (unsigned int)gg <<  "]= " << minLCPsucToFind[(unsigned int)gg]  << "\n";
                    //if (minLCPsuc[gg] != maxValueLen) {      //There are an LCP interval opened for the symbol c_g
              if (lcpValueNow < minLCPsuc[(unsigned int)gg])  //comparison with the last inserted lcp
                minLCPsuc[(unsigned int)gg] = lcpValueNow;
              //std::cerr << "---Double check minLCPsuc[" << (unsigned int)gg <<  "]= " << minLCPsuc[(unsigned int)gg]  << "\n";
                        //std::cerr << "---     minLCPsucText[" << (unsigned int)gg <<  "]= " << minLCPsucText[(unsigned int)gg]  << "\n";
                      //std::cerr << "---     minLCPsucToFind[" << (unsigned int)gg <<  "]= " << minLCPsucToFind[(unsigned int)gg]  << "\n";
              //}
              }

            if (minLCPcurFound[(unsigned int)gg] == 1) {    //We have to compute the minimum for the lcp of the symbol gg
              if (lcpValueNow < minLCPcur[(unsigned int)gg])  //comparison with the last inserted lcp
                minLCPcur[(unsigned int)gg] = lcpValueNow;
            }
          }

          //I have to re-set
          minLCPcur[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = maxValueLen;    //I initialized the minLCPcur
          minLCPcurFound[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = 1;    // I set the fact that I met the new symbol

///////////////////////////////////////////////////////////////////
/*
          //If the new symbol does not close an LCP interval
          if (minLCPsuc[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] == maxValueLen) {       //it means that there is not an interval lcp to compute for the new symbol
            //The new symbol opens an LCP interval
//            minLCPsuc[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = vectTriple[k].lcpCurN+1;   //It sets the min_2 for successive symbol with the current of the new symbol (next text)

            minLCPsucText[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = 0;
//            std::cerr << "minSuc was the maxValue. The new value for minLCPsuc[" << newSymb[vectTriple[k].seqN] << "] is " << minLCPsuc[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] << "\n";
            //std::cerr << "minSuc is the maxValue yet. minLCPsucText is " << minLCPsucText[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] << "\n";
          }
          else {
*/
          if (minLCPsucToFind[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] == 1) { //If the new symbol closes a LCP interval
                    //std::cerr << "---2a) (dovrebbe essere <MAX) minLCPsuc[" << newSymb[vectTriple[k].seqN] << "]= "<< minLCPsuc[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] << "\n";
                  //std::cerr << "---     minLCPsucText[" << newSymb[vectTriple[k].seqN] <<  "]= " << minLCPsucText[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]]  << "\n";
                //std::cerr << "---     minLCPsucToFind[" << newSymb[vectTriple[k].seqN] <<  "]= " << minLCPsucToFind[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]]  << "\n";
          //if (minLCPsuc[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] != maxValueLen) {
            //I have already computed the minimum in the previous FOR (close the previous lcp interval).
          //I can set lcpSucN of minLCPsucText[alpha[(unsigned int)[vectTriple[k].seqN]]]
            vectTriple[minLCPsucText[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]]].lcpSucN  = minLCPsuc[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]];
          }
          //It set minLCP for the new LCP interval for the new symbol
//          minLCPsuc[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = vectTriple[k].lcpCurN+1;   //It sets the min_2 for successive symbol with the current of the new symbol (next text)
          minLCPsuc[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = maxValueLen;
          minLCPsucText[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = k;
          minLCPsucToFind[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = 1;  //I have to open the lcp succ for this sequence
          //std::cerr << "minSuc was not the maxValue. The new value for minLCPsuc[" << newSymb[vectTriple[k].seqN] << "] is " << minLCPsuc[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] << "\n";

          //Since we have inserted, we must update them
          cont++;    //number of read symbols
          toRead--;

                 
          //Now I have to update the lcp of the next symbol, if it exists.
          //std::cerr << "----*We have computed minLCPsuc. Now we verify that the new next symbol is in the successive position.\n";
          //std::cerr << "----*newSymb[vectTriple[k].seqN= " << (char)newSymb[vectTriple[k].seqN] << ".\n";

                    // The next symbol must be newSymb[vectTriple[k+1].seqN] != TERMINATE_CHAR_LEN, otherwise we must ignore it
//          if(newSymb[vectTriple[k+1].seqN] != TERMINATE_CHAR_LEN) { 

            //std::cerr << "----==CurrentPile " << (unsigned int)currentPile << "\n";
            //std::cerr << "---- k=" << k << " triple["<<k<<"].lcpSucN="<< vectTriple[k].lcpSucN << "\n";
            if ((k + 1 < nExamedTexts) && (vectTriple[k+1].pileN == currentPile) && (vectTriple[k+1].posN == vectTriple[k].posN + 1)) {
              //std::cerr << "----2 k=" << k << " triple["<<k<<"].lcpSucN="<< vectTriple[k].lcpSucN << "\n";
              //If the next symbol is the new symbol associated with another text (in the same segment), that I have to insert yet
              //I can ignored this updating, because the lcp of the new symbol is already computed
              //We set the minLCPsuc with the value LCP associated with the next symbol that we have to insert in it.
              //it should be vectTriple[k].lcpSucN + 1 == vectTriple[k+1].lcpCurN
              if (vectTriple[k].lcpSucN + 1 != vectTriple[k+1].lcpCurN +1) {
                std::cerr << "-->Warning! Should be the same? posSymb = " << (unsigned int)posSymb;
                std::cerr << " Segment: " << (unsigned int)vectTriple[k+1].pileN << " == " << (unsigned int)currentPile;
                std::cerr <<" Position: " << vectTriple[k+1].posN << " == " << vectTriple[k].posN + 1 << "\n";
                std::cerr << " triple["<<k<<"].lcpSucN(="<< vectTriple[k].lcpSucN << ") + 1= " << vectTriple[k].lcpSucN + 1 << " == triple["<<k+1<<"].lcpCurN(="<< vectTriple[k+1].lcpCurN << ")+1= " << vectTriple[k+1].lcpCurN+1<< " ";
                std::cerr << ", Seq k N. " << vectTriple[k].seqN << " and Seq k+1 N. " << vectTriple[k+1].seqN << "\n";
              }

              //Hence, at the next step, I have to insert the symbol newSymb[vectTriple[k+1].seqN]
              //I check if newSymb[vectTriple[k+1].seqN] is equal to the inserted symbol now, that is newSymb[vectTriple[k].seqN]
              if (newSymb[vectTriple[k].seqN] == newSymb[vectTriple[k+1].seqN]) {
                //std::cerr << "----A) newSymb[vectTriple[k].seqN]=" << newSymb[vectTriple[k].seqN] << " newSymb[vectTriple[k+1].seqN="<< newSymb[vectTriple[k+1].seqN] << "\n";
                //std::cerr << "----A) k=" << k << " triple["<<k<<"].lcpSucN="<< vectTriple[k].lcpSucN << "\n";
                //In this case, I can set the lcpSuc of newSymb[vectTriple[k].seqN] to vectTriple[k+1].lcpCurN + 1
                
                //WARNING: vectTriple[k].lcpSucN = vectTriple[k+1].lcpCurN +1; 
                if (currentPile == 0)
                  vectTriple[k].lcpSucN = vectTriple[k+1].lcpCurN;       //I set the lcpSucN of the current symbol (in position k)
                else
                  vectTriple[k].lcpSucN = vectTriple[k+1].lcpCurN+1;       //I set the lcpSucN of the current symbol (in position k)
                
                //std::cerr << "----AA) k=" << k << " triple["<<k<<"].lcpSucN="<< vectTriple[k].lcpSucN << "\n";
                //I close the LCP interval for newSymb[vectTriple[k].seqN]], so
                minLCPsuc[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = maxValueLen;
                minLCPsucToFind[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = 0;  //closes the LCP interval
                minLCPsucText[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = 0;
              }
              else  {
                //std::cerr << "----B) k=" << k << " triple["<<k<<"].lcpSucN="<< vectTriple[k].lcpSucN << "\n";
                //In this case, I cannot set the lcpSuc of newSymb[vectTriple[k].seqN], because the symbol corresponding to k+1 is different
                //I set minLCPsuc of newSymb[vectTriple[k].seqN] to vectTriple[k+1].lcpCurN +1, and I search the symbol newSymb[vectTriple[k].seqN]
                minLCPsuc[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = vectTriple[k+1].lcpCurN +1; //set the lcp interval
                minLCPsucToFind[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = 1;
                minLCPsucText[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = k;
              }
            }
            else {
              //std::cerr << "----3 k=" << k << " triple["<<k<<"].lcpSucN="<< vectTriple[k].lcpSucN << "\n";
            //The next symbol is in the current segment, if there exists, it is an old symbol
              uchar sucSymbol='\0';
              //it there exists another symbol in the segment file then I have to update it.
              numchar = fread(&sucSymbol,sizeof(uchar),1,InFileBWT);
              if (numchar == 1) {  //There exists at least a symbol in the current segment
                //std::cerr << "I read the symbol sucSymbol is " << (unsigned int)alpha[(unsigned int)sucSymbol] << "\n";

                numcharWrite = fwrite (&sucSymbol, sizeof(uchar), numchar, OutFileBWT);
                //if (verboseEncode == 1)
                //  std::cerr << "In BWT: number read " << numchar << " ==1 \n";
                assert(numchar == numcharWrite); // we should always read/write the same number of characters

                #if BUILD_SA == 1
                  numcharSA = fread(bufferSA,sizeof(ElementType),1,InFileSA);
                  assert(numcharSA == 1); // we should always read/write the same number of characters
                  numcharWriteSA = fwrite (bufferSA, sizeof(ElementType), numcharSA , OutFileSA);
                  assert(numcharSA == numcharWriteSA); // we should always read/write the same number of characters
                #endif

                //I have to update the lcp of the next symbol
                dataTypelenSeq sucLCP = 0;
                numchar = fread(&sucLCP,sizeof(dataTypelenSeq),1,InFileLCP);   //I have to change it
                //if (verboseEncode == 1)
                //  std::cerr << "In LCP: number read " << numchar << " ==1 \n";
                assert(numchar == 1); // we should always read/write the same number of characters
                //I have to update the lcp of this symbol and I have to copy it into the new bwt segment

                dataTypelenSeq lcpValueNow=0;
                if (vectTriple[k].pileN == 0) { //In 0-pile the value LCP is 0
                  lcpValueNow = 0;
                }
                else {
                  lcpValueNow = vectTriple[k].lcpSucN + 1;
                }
                //std::cerr << "---#############lcp was " << sucLCP << " and now is vectTriple["<< k << "].lcpSucN + 1=lcpValueNow= " << lcpValueNow << "\n";
                numcharWrite = fwrite (&lcpValueNow ,sizeof(dataTypelenSeq), numchar , OutFileLCP);    //Updated the lcpSuc
                assert(numchar == numcharWrite); // we should always read/write the same number of characters


                //Now, I have to check if the symbol sucSymbol close the LCP interval the new symbol
                if (newSymb[vectTriple[k].seqN] == sucSymbol) {
                  //std::cerr << "---The (read) succSymb is equal to the new symbol\n";
                  //If it is equal to newSymb[vectTriple[k].seqN] then I can set lcpSucN of newSymb[vectTriple[k].seqN]
                  vectTriple[k].lcpSucN = lcpValueNow;
                  minLCPsucToFind[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = 0;  //Close the LCP interval
                  minLCPsuc[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = maxValueLen;
                  //std::cerr << "---#######vectTriple["<< k << "].lcpSucN = maxValueLen = " << maxValueLen << "\n";
                  minLCPsucText[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = 0;
                }
                else {
                  //std::cerr << "The succSymb is not equal to the new symbol\n";
                  minLCPsucToFind[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = 1;  //I have to search the symbol newSymb[vectTriple[k].seqN]]
                  minLCPsuc[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = lcpValueNow;  //I set the minLCPsuc for the new symbol
                  minLCPsucText[alpha[(unsigned int)newSymb[vectTriple[k].seqN]]] = k;
                  //std::cerr << "The new symbol is " << newSymb[vectTriple[k].seqN] << " and vectTriple.sucN[" << k <<  "] is " << vectTriple[k].lcpSucN << std::endl;
                  //It could close the LCP interval for the symbol sucSymb if it is opened
                  //If the symbol sucSymbol does not close an LCP interval, ok!
                  if (minLCPsucToFind[alpha[(unsigned int)sucSymbol]] == 1) {    //We have to compute the minimum for the lcp of the symbol (unsigned int)sucSymbol
                    //it means that there is an interval lcp to close for the symbol (unsigned int)sucSymbol
                    //std::cerr << "4) (dovrebbe essere <MAX) minLCPsuc[" << sucSymbol <<  "]= " << minLCPsuc[alpha[sucSymbol]]  << "\n";
                    //if (minLCPsuc[alpha[(unsigned int)sucSymbol]] != maxValueLen) {
                    //std::cerr << "The symbol " << (unsigned int)sucSymbol << " closes the LCP interval of the symbol. So \n";
                    //The symbol sucSymbol closes a LCP interval
                    if (lcpValueNow < minLCPsuc[alpha[(unsigned int)sucSymbol]])  //comparison with the last inserted lcp
                      minLCPsuc[alpha[(unsigned int)sucSymbol]] = lcpValueNow;
                    vectTriple[minLCPsucText[alpha[(unsigned int)sucSymbol]]].lcpSucN  = minLCPsuc[alpha[(unsigned int)sucSymbol]];  //I can set lcpSucN of minLCPsucText[alpha[(unsigned int)[sucSymbol]]
                    //std::cerr << "---We update sucN to " << vectTriple[minLCPsucText[alpha[(unsigned int)sucSymbol]]].lcpSucN << std::endl;
                    //It closes the LCP interval, so
                    minLCPsucToFind[alpha[(unsigned int)sucSymbol]] = 0;   //Close the LCP interval for the symbol sucSymbol
                    minLCPsuc[alpha[(unsigned int)sucSymbol]] = maxValueLen;
                    minLCPsucText[alpha[(unsigned int)sucSymbol]] = 0;
                  }
                }


                //Now, I have to update the lcpSucc of the opened interval lcp
                for (dataTypedimAlpha gg = 0 ; gg < sizeAlpha; gg++) {
                  //std::cerr << "---+++ Before: Symbol " << (unsigned int)gg << " minSuc is " << minLCPsuc[gg] << std::endl;
                  //Update the minLCPsuc
                  if (minLCPsucToFind[(unsigned int)gg] == 1) {    //We have to compute the minimum for the lcp of the symbol gg
                    //std::cerr << "---3) (dovrebbe essere <MAX) minLCPsuc[" << (unsigned int)gg <<  "]= " << minLCPsuc[(unsigned int)gg]  << "\n";
                    //if (minLCPsuc[gg] != maxValueLen) {      //There are an LCP interval opened for the symbol c_g
                    if (lcpValueNow < minLCPsuc[(unsigned int)gg]) { //comparison with the last inserted lcp
                      minLCPsuc[(unsigned int)gg] = lcpValueNow;
                      //std::cerr << "---New min2: The symbol is " << (unsigned int)gg << " and minLCPsuc[" << gg <<  "] is " << minLCPsuc[(unsigned int)gg] << std::endl;
                    }
                  }
                  //Update the minLCPCur
                  if (minLCPcurFound[gg]==1)   //I already met the symbol gg. So, I must compute the minimum
                    //std::cerr << "Before: Symbol " << (unsigned int)gg << " minCur is " << minLCPcur[gg] << std::endl;
                    if (lcpValueNow < minLCPcur[gg]) { //comparison with the last inserted lcp for the symbol gg
                      minLCPcur[gg] = lcpValueNow;
                      //std::cerr << "New min1: The symbol is " << (unsigned int)gg << " and minLCPcur[" << gg <<  "] is " << minLCPcur[(unsigned int)gg] << std::endl;
                    }
                }

                //For the current LCP
                minLCPcur[alpha[(unsigned int)sucSymbol]] = maxValueLen; //For each occurrence of sucSymbol, I have to set the min1(the interval starts from the next symbol)
                minLCPcurFound[alpha[(unsigned int)sucSymbol]] = 1;  //So I open the LCP interval for sucSymbol (for the first occurrence of sucSymbol)

                //We have read another symbol of bwt and its associated lcp
                cont++;    //number of read symbols
                toRead--;
              }
              else {  //Then there are not other symbols.
                //it means that the file does not contain other symbols and we have inserted the symbol in the last position
                //all lcp intervals have to be close and initializate
                for (dataTypedimAlpha gg = 0 ; gg < sizeAlpha; gg++) {
                  if (minLCPsucToFind[(unsigned int)gg] == 1) {    //We have to close the lcp interval of the symbol gg

                    //I have to set the lcpSuc of the text minLCPsucText[(unsigned int)gg] to 0
                    vectTriple[minLCPsucText[(unsigned int)gg]].lcpSucN  = 0;
                    minLCPsucToFind[(unsigned int)gg] = 0;
                    minLCPsuc[(unsigned int)gg] = maxValueLen;
                    minLCPsucText[(unsigned int)gg] = 0;
                    //std::cerr << "---++ Symbol " << (unsigned int)gg << " minCur is " << minLCPsuc[gg] << std::endl;
                  }
                }
              } //end else
            }   //end else

//          } 

                }   //end if (toRead==0)
        //}

          if (verboseEncode==1) {
            std::cerr << "\nFine inserimento del simbolo= " << newSymb[vectTriple[k].seqN] << " della sequenza " << vectTriple[k].seqN << std::endl;
                    std::cerr << std::endl;
          std::cerr << "Q  ";
          for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
            std::cerr << (unsigned int)vectTriple[g].pileN << " ";
            }
          std::cerr << std::endl;
          std::cerr << "P  ";
                  for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
                    std::cerr << vectTriple[g].posN  << " ";
              }
            std::cerr << std::endl;
          std::cerr << "N  ";
          for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
                      std::cerr << vectTriple[g].seqN  << " ";
                }
              std::cerr << std::endl;
            std::cerr << "C  ";
          for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
            std::cerr << vectTriple[g].lcpCurN  << " ";
          }
                  std::cerr << std::endl;
                std::cerr << "S  ";
              for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
                std::cerr << vectTriple[g].lcpSucN  << " ";
          }
          std::cerr << std::endl;
        }

       #if (BCR_SET==1)
      }  //close if (newSymb[vectTriple[k].seqN] != TERMINATE_CHAR_LEN)
       #endif
            //}


      k++;   //  I changed the number of the sequence. New iteration.
    }    //close while ((k< nText) && (vectTriple[k].pileN == currentPile))

    //std::cerr << "End symbols to insert in the current segment\n" << std::endl;
    //it means that posN[k]<>currentPile, so I have to change BWT-file
    //But before, I have to copy the remainder symbols from the old BWT to new BWT
    //We could need to compute the minLCPsuc for some text

        
    while (numchar!=0) {
            //std::cerr << "Controlla se ci sono altri simboli nel segmento (while)!\n";
      //For BWT
      numchar = fread(buffer,sizeof(uchar),SIZEBUFFER,InFileBWT);
      numcharWrite = fwrite (buffer, sizeof(uchar), numchar , OutFileBWT);
      assert(numchar == numcharWrite); // we should always read/write the same number of characters

      #if BUILD_SA == 1
        numcharSA = fread(bufferSA,sizeof(ElementType),SIZEBUFFER,InFileSA);
        numcharWriteSA = fwrite (bufferSA, sizeof(ElementType), numcharSA , OutFileSA);
        assert(numcharSA == numcharWriteSA);
      #endif


        //For LCP
      numchar = fread(bufferLCP,sizeof(dataTypelenSeq),SIZEBUFFER,InFileLCP);
      numcharWrite = fwrite (bufferLCP, sizeof(dataTypelenSeq), numchar , OutFileLCP);
      assert(numchar == numcharWrite); // we should always read/write the same number of characters
        //Compute lcpSucN for the other texts
        //For each symbol in the buffer, we check it it close any interval, while each entry in minLcpSuc is maxValue
        //TO DO: TO OPTIMIZE. IT CAN END BEFORE. IT DOES NOT NEED TO READ THE ENTIRE BUFFER

      for (dataTypeNChar bb = 0 ; bb < numchar; bb++) {
        //std::cerr << "Symbol in the buffer is " << (unsigned int)buffer[bb] << "\n";
        //First, I check if the symbol bb closes a previous interval or it is in the middle or no.
        //In any case, for each symbol, we have to update the minimum
        for (dataTypedimAlpha gg = 0 ; gg < sizeAlpha; gg++) {
                    if (minLCPsucToFind[(unsigned int)gg] == 1) {    //We have to compute the minimum for the lcp of the symbol gg
                        //std::cerr << "6) (dovrebbe essere <MAX) minLCPsuc[" << (unsigned int)gg <<  "]= " << minLCPsuc[(unsigned int)gg]  << "\n";
                        //std::cerr << "+Symbol " << (unsigned int)gg << " minSuc is " << minLCPsuc[gg] << " and bufferLCP[bb] is " << bufferLCP[bb] << std::endl ;
              //if (minLCPsuc[gg] != maxValueLen) {      //There are an LCP interval apened for the symbol c_g
            if (bufferLCP[bb] < minLCPsuc[(unsigned int)gg]) { //comparison with the last inserted lcp
              //std::cerr << "bufferLCP["<< (unsigned int)bb <<"] is "<< bufferLCP[bb] << " minLCPsuc["<< (unsigned int)gg << "] is " << minLCPsuc[gg] << std::endl;
              minLCPsuc[(unsigned int)gg] = bufferLCP[bb];
              //std::cerr << "bufferLCP["<< (unsigned int)bb <<"] is "<< bufferLCP[bb] << " minLCPsuc["<< (unsigned int)gg << "] is " << minLCPsuc[gg] << std::endl;
            }
          }
        }
        //std::cerr << "Considered Symbol " << (unsigned int)bb << ". minSuc is " << minLCPsuc[alpha[(unsigned int)bb]] << std::endl;

        //If the symbol buffer[bb] does not close an LCP interval, ok!
                if (minLCPsucToFind[alpha[(unsigned int)buffer[bb]]] == 1) {    //We have to compute the minimum for the lcp of the symbol gg
                    //std::cerr << "---7) (dovrebbe essere <MAX) minLCPsuc[" << (unsigned int)buffer[bb] <<  "]= " << minLCPsuc[alpha[(unsigned int)buffer[bb]]]  << "\n";
                    //if (minLCPsuc[alpha[(unsigned int)buffer[bb]]] != maxValueLen) {       //it means that there is an interval lcp to close for this symbol
            //The symbol bb closes a LCP interval
            //I have already computed the minimum (close the previous lcp interval).
            //I can set lcpSucN of minLCPsucText[alpha[(unsigned int)[bb]]
            vectTriple[minLCPsucText[alpha[(unsigned int)buffer[bb]]]].lcpSucN  = minLCPsuc[alpha[(unsigned int)buffer[bb]]];
            //std::cerr << "---Updated sucN: Symbol " << (unsigned int)buffer[bb] << " minSuc is " << minLCPsuc[alpha[(unsigned int)buffer[bb]]] << std::endl;
            //It close the LCP interval, so
                        minLCPsucToFind[alpha[(unsigned int)buffer[bb]]] = 0;
                        minLCPsuc[alpha[(unsigned int)buffer[bb]]] = maxValueLen;
            minLCPsucText[alpha[(unsigned int)buffer[bb]]] = 0;
        }
       }  //For
//           } //If       //***********************************
    }   //close while (numchar!=0)

        //The file is finished! but some interval lcp could be opened
        //In this case, we have to set the lcpSuc to 0
    for (dataTypedimAlpha gg = 0 ; gg < sizeAlpha; gg++) {
      //std::cerr << "++++++++++++++++Symbol " << (unsigned int)gg << " minSuc is " << minLCPsuc[(unsigned int)gg] << std::endl;
            if (minLCPsucToFind[(unsigned int)gg] == 1) {    //We have to close the lcp interval of the symbol gg
                //std::cerr << "8) (dovrebbe essere <MAX) minLCPsuc[" << (unsigned int)gg <<  "]= " << minLCPsuc[(unsigned int)gg] << "\n";
                //if (minLCPsuc[gg] != maxValueLen) {      //There are an LCP interval opened for the symbol c_g
            //std::cerr << " Lcp interval opened for Seq " << (unsigned int)minLCPsucText[(unsigned int)gg] << " minSuc is " << minLCPsuc[(unsigned int)gg] << " minSucToFind is " << minLCPsucToFind[(unsigned int)gg] << std::endl;
           vectTriple[minLCPsucText[(unsigned int)gg]].lcpSucN = 0;
               minLCPsucToFind[(unsigned int)gg] = 0;
               minLCPsuc[(unsigned int)gg] = maxValueLen;
           minLCPsucText[(unsigned int)gg] = 0;
         //std::cerr << " Now Lcp Succ for Sequence " << (unsigned int)minLCPsucText[gg] << "  is " << minLCPsuc[gg] << std::endl;
           }
        }

    fclose(InFileBWT);
    fclose(OutFileBWT);
    fclose(InFileLCP);
    fclose(OutFileLCP);

    #if BUILD_SA == 1
      fclose(InFileSA);
      fclose(OutFileSA);
    #endif

    //Rename new file in old file
    int numchar1=sprintf (filename, "bwt_%d", currentPile);
    numchar1=sprintf (filenameIn,"%s%s",filename,ext);
    numchar1=sprintf (filenameOut,"new_%s%s",filename,ext);
    //std::cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << std::endl;
    OutFileBWT = fopen(filenameOut, "rb");
    if (OutFileBWT!=NULL) { //If it exists
      fclose(OutFileBWT);
      if (remove(filenameIn)!=0)
        std::cerr << filenameIn <<": Error deleting file" << std::endl;
      else
        if(rename(filenameOut,filenameIn))
          std::cerr << filenameOut <<": Error renaming " << std::endl;
    }

    numchar1=sprintf (filenameLCP, "lcp_%d", currentPile);
    numchar1=sprintf (filenameInLCP,"%s%s",filenameLCP,ext);
    numchar1=sprintf (filenameOutLCP,"new_%s%s",filenameLCP,ext);
    //std::cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << std::endl;
    OutFileLCP = fopen(filenameOutLCP, "rb");
    if (OutFileLCP!=NULL) { //If it exists
      fclose(OutFileLCP);
      if (remove(filenameInLCP)!=0)
        std::cerr << filenameInLCP <<": Error deleting file" << std::endl;
      else
        if(rename(filenameOutLCP,filenameInLCP))
          std::cerr << filenameOutLCP <<": Error renaming " << std::endl;
    }

    #if BUILD_SA == 1
      //Renaming new to old
      for (dataTypedimAlpha g = 0 ; g < sizeAlpha; g++) {
        numchar=sprintf (filename, "sa_%d", g);
        numchar=sprintf (filenameIn,"%s%s",filename,ext);
        numchar=sprintf (filenameOut,"new_%s%s",filename,ext);
        //std::cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << std::endl;
        OutFileSA = fopen(filenameOut, "rb");

        if (OutFileSA!=NULL) { //If it exists
          fclose(OutFileSA);
          if (remove(filenameIn)!=0)
            std::cerr << filenameIn <<": Error deleting file" << std::endl;
          else
            if(rename(filenameOut,filenameIn))
              std::cerr << filenameOut <<": Error renaming " << std::endl;
        }


      }
    #endif

    j=k;
  }      // close while (j < nText)

  #if verboseEncode==1
    std::cerr << "After the computation of LCP for the next iteration" << std::endl;
    std::cerr << "Q  ";
    for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
      std::cerr << (unsigned int)vectTriple[g].pileN << " ";
    }
    std::cerr << std::endl;
    std::cerr << "P  ";
    for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
      std::cerr << vectTriple[g].posN  << " ";
    }
    std::cerr << std::endl;
    std::cerr << "N  ";
    for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
      std::cerr << vectTriple[g].seqN  << " ";
    }
    std::cerr << std::endl;
    std::cerr << "C  ";
    for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
      std::cerr << (unsigned int)vectTriple[g].lcpCurN  << " ";
    }
    std::cerr << std::endl;
    std::cerr << "S  ";
    for (dataTypeNSeq g = 0 ; g < nExamedTexts; g++) {
      std::cerr << (unsigned int)vectTriple[g].lcpSucN  << " ";
    }
    std::cerr << std::endl;
  #endif

  //Renaming new to old
  for (dataTypedimAlpha g = 0 ; g < sizeAlpha; g++) {
    numchar=sprintf (filename, "bwt_%d", g);
    numchar=sprintf (filenameIn,"%s%s",filename,ext);
    numchar=sprintf (filenameOut,"new_%s%s",filename,ext);
    //std::cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << std::endl;
    OutFileBWT = fopen(filenameOut, "rb");
    if (OutFileBWT!=NULL) { //If it exists
      fclose(OutFileBWT);
      if (remove(filenameIn)!=0)
        std::cerr << filenameIn <<": Error deleting file" << std::endl;
      else
        if(rename(filenameOut,filenameIn))
          std::cerr << filenameOut <<": Error renaming " << std::endl;
    }

    numchar=sprintf (filenameLCP, "lcp_%d", g);
    numchar=sprintf (filenameInLCP,"%s%s",filenameLCP,ext);
    numchar=sprintf (filenameOutLCP,"new_%s%s",filenameLCP,ext);
    //std::cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << std::endl;
    OutFileLCP = fopen(filenameOutLCP, "rb");
    if (OutFileLCP!=NULL) { //If it exists
      fclose(OutFileLCP);
      if (remove(filenameInLCP)!=0)
        std::cerr << filenameInLCP <<": Error deleting file" << std::endl;
      else
        if(rename(filenameOutLCP,filenameInLCP))
          std::cerr << filenameOutLCP <<": Error renaming " << std::endl;
    }


  }

  delete [] buffer;
  delete [] filenameIn;
  delete [] filename;
  delete [] filenameOut;
  delete [] bufferLCP;
  delete [] filenameInLCP;
  delete [] filenameLCP;
  delete [] filenameOutLCP;

  #if BUILD_SA == 1
    delete [] bufferSA;
    delete [] filenameInSA;
    delete [] filenameSA;
    delete [] filenameOutSA;
  #endif
}
#endif

#if OUTPUT_FORMAT_EGSA == 1

int BCRexternalBWT::storeEGSAoutput( const char* fn ) {
  
  
  static FILE  *InFileBWT;                  // input file BWT;
  static FILE  *InFileLCP;                  // input file LCP;
  static FILE  *InFilePairSA;                  //  input file SA;
  
  char *filenameIn = new char[1024];
  char *filename = new char[1024];
  const char *ext = ".aux";
  
  //Open EGSA file
  
    //dataTypeNSeq numTotalTexts;
    //#if BUILD_BCR_FROM_EGSA == 1
    //  numTotalTexts = nText+nAddedTextEGSA;
    //#else
    //  numTotalTexts = nText;
    //#endif
  
  FILE *f_ESA;      // pointer to the ESA input file
  char c_aux[1024];
//  sprintf(c_aux, "%s.%d.gesa", fn, numTotalTexts);  //
  sprintf(c_aux, "%s.%d.gesa", fn, 0);
  f_ESA = fopen(c_aux, "wb");//rb
  if (f_ESA == NULL) {
    std::cerr << "EGSA: Error opening: " << c_aux << std::endl;
    exit (EXIT_FAILURE);
  } 
  t_GSA GSA;
  
  std::cerr << "EGSA: sizeof(type of t_GSA): " << sizeof(t_GSA) << "\n";
  std::cerr << "file EGSA: " << lengthTot_plus_eof * sizeof(t_GSA) << "\n";

  dataTypeNChar numcharBWT, numcharPairSA, numcharLCP, numEle=0;
  //dataTypeNChar numcharWrite=0;
  
  uchar *bufferBWT = new uchar[SIZEBUFFER];
  dataTypelenSeq *bufferLCP = new dataTypelenSeq[SIZEBUFFER];
  ElementType *bufferSA = new ElementType[SIZEBUFFER];


  std::cerr << "Build the entire BWT/LCP/SA file (eGSA format).\n";

  #if (printEGSA == 1)
    dataTypeNChar *freqOut = new dataTypeNChar [SIZE_ALPHA];
    for (dataTypedimAlpha i = 0; i < SIZE_ALPHA-1; ++i)
         freqOut[i] = 0;
    freqOut[SIZE_ALPHA-1] = 0;

    cout << "bwt\t" << "gsa(pos,id)\t" << "lcp\t" << endl; 
  #endif
    
  
  dataTypeNChar numberBWTwritten = 0;
  
  for (dataTypedimAlpha g = 0 ; g < sizeAlpha; g++) {
    #if (printEGSA == 1)
      cout << endl; 
    #endif
  
    //BWT
    numcharBWT=sprintf (filename, "bwt_%d", g);
    numcharBWT=sprintf (filenameIn,"%s%s",filename,ext);
    InFileBWT = fopen(filenameIn, "rb");
    if (InFileBWT==NULL) {
      std::cerr << "BWT file " << (unsigned int)g <<": Error opening " << std::endl;
      exit (EXIT_FAILURE);
    }
    //LCP
    numcharLCP=sprintf (filename, "lcp_%d", g);
    numcharLCP=sprintf (filenameIn,"%s%s",filename,ext);
    InFileLCP = fopen(filenameIn, "rb");
    if (InFileLCP==NULL) {
      std::cerr << "storeEGSAoutput: LCP file " << (unsigned int)g <<": Error opening " << std::endl;
      exit (EXIT_FAILURE);
    }
    //GSA
    numcharPairSA=sprintf (filename, "sa_%d", g);
    numcharPairSA=sprintf (filenameIn,"%s%s",filename,ext);
    InFilePairSA = fopen(filenameIn, "rb");
    if (InFilePairSA==NULL) {
      std::cerr << "storeEGSAoutput: GSA file " << (unsigned int)g <<": Error opening " << std::endl;
      exit (EXIT_FAILURE);
    }
    
    while ((!feof(InFileBWT))  && (!feof(InFilePairSA)) && (!feof(InFileLCP)) )  {
      numcharBWT = fread(bufferBWT,sizeof(uchar),SIZEBUFFER,InFileBWT);
      numcharPairSA = fread(bufferSA,sizeof(ElementType),SIZEBUFFER,InFilePairSA);
      numcharLCP = fread(bufferLCP,sizeof(dataTypelenSeq),SIZEBUFFER,InFileLCP);
      if ((numcharLCP != numcharPairSA) || (numcharLCP != numcharBWT))
        std::cerr << "Error: number  in BWT in Pair SA  and in LCP\n";
      else {
        //std::cerr << "bwt\tlcp\tpos\tnumSeq\tSA\tQS\n";
        for (dataTypeNChar i=0; i < numcharBWT; i++) {
          numEle++;

          GSA.suff = bufferSA[i].sa;
          GSA.text = bufferSA[i].numSeq;
          
          GSA.lcp = bufferLCP[i];

          if (bufferBWT[i] == TERMINATE_CHAR)
            GSA.bwt='\0';
          else
            GSA.bwt=bufferBWT[i];

          fwrite(&GSA.text, sizeof(dataTypeNSeq), 1, f_ESA);  
          fwrite(&GSA.suff, sizeof(dataTypelenSeq), 1, f_ESA);  
          fwrite(&GSA.lcp, sizeof(dataTypelenSeq), 1, f_ESA); 
          fwrite(&GSA.bwt, sizeof(uchar), 1, f_ESA);  
          numberBWTwritten++;
          
          #if (printEGSA == 1)
            cout << (char)bufferBWT[i] << "\t" << (int)bufferSA[i].sa << "\t" << (int)bufferSA[i].numSeq << "\t" << (int)bufferLCP[i] << "\t" << std::endl;
            freqOut[(unsigned int)(bufferBWT[i])]++;
          #endif
        }
      }
    }  //end-while
    
        
    fclose(InFilePairSA);
    fclose(InFileLCP);
    fclose(InFileBWT);
      
      //Delete or rename the partial BWT bwt_%d????
      #if deletePartialBWT == 1
        if (remove(filenameIn)!=0)
          std::cerr << "storeEGSAoutput: Error deleting file" << std::endl;

        /*
        #else //rename the aux bwt file
          int lung = strlen(fn) + strlen(filenameIn)+1;
          char *newfilename = new char[lung];
          numchar=sprintf (newfilename,"%s%s",fn,filenameIn);
          //std::cerr  << newfilename << " " << filenameIn << std::endl;
          if(rename(filenameIn, newfilename))
            std::cerr  <<"storeEntireBWTFilePartial: Error renaming file" << std::endl;
        delete[] newfilename;
        */
      #endif
      
  } // end-for
  std::cerr  <<"Number records in EGSA file: " << numberBWTwritten << std::endl;
  assert ( numberBWTwritten== lengthTot_plus_eof);
  
  #if (printEGSA == 1)
    int check=0;
    std::cerr << "Distribution in BWT\n";
    for (dataTypedimAlpha g = 0 ; g < SIZE_ALPHA-1; g++) {
      if (freqOut[g] > 0){
        InFileBWT = openFilePartialIn(alpha[(unsigned int)g]);
        fseek(InFileBWT,0,SEEK_END);
        dataTypeNChar lengthBWTPartial=ftell(InFileBWT);
        if (freqOut[g] == lengthBWTPartial)
          std::cerr << (unsigned int)g << " " << g << " freq: " << freqOut[g] << "\n";
        else {
          std::cerr << (unsigned int)g << " " << g << " freq: " << freqOut[g] << " freq in BWT: " << lengthBWTPartial << "***********************\n";
          check = 1;
        }
        closeFilePartial(InFileBWT);
      }
    }
    if (freqOut[(SIZE_ALPHA-1)] > 0){
      InFileBWT = openFilePartialIn(alpha[(unsigned int)(SIZE_ALPHA-1)]);
      fseek(InFileBWT,0,SEEK_END);
      dataTypeNChar lengthBWTPartial=ftell(InFileBWT);
      if (freqOut[(SIZE_ALPHA-1)] == lengthBWTPartial)
        std::cerr << (unsigned int)(SIZE_ALPHA-1) << " " << SIZE_ALPHA-1 << " freq: " << freqOut[(SIZE_ALPHA-1)] << "\n";
      else {
        std::cerr << (unsigned int)(SIZE_ALPHA-1) << " " << SIZE_ALPHA-1 << " freq: " << freqOut[(SIZE_ALPHA-1)] << " freq in BWT: " << lengthBWTPartial << "***********************\n";
        check = 1;
      }
      closeFilePartial(InFileBWT);
    }

    if (check == 1)
      std::cerr << "WARNING! Some length partial BWT != the number of occorrences! " << std::endl;

    delete [] freqOut;
  #endif

  delete [] filenameIn;
  delete [] filename;
  delete [] bufferSA; 
  delete [] bufferBWT;
  delete [] bufferLCP; 
  fclose(f_ESA);
  
  return 0;
}


int BCRexternalBWT::storeEGSAoutputFromEntireFiles (string input) {

  string fnBWT, fnPairSA, fnLCP, fnEGSA;

    fnBWT = input + ".out\0";
  fnLCP = input + ".out.lcp\0";
  fnPairSA = input + ".out.pairSA\0";
  fnEGSA = input + "." + "0" + ".gesa\0";
  
  std::cout << "file BWT : "  << fnBWT <<  "." << std::endl;
  std::cout << "file LCP: "  << fnLCP <<  "." << std::endl;
  std::cout << "file PairSA: "  << fnPairSA <<  "." << std::endl;
  std::cout << "file EGSA: "  << fnEGSA <<  "." << std::endl;       



  //Open EGSA file
  
  FILE *f_ESA;      // pointer to the ESA input file

  
  f_ESA = fopen(fnEGSA.c_str(), "wb");
  if (f_ESA == NULL) {
    std::cerr << "EGSA: Error opening: " << fnEGSA << std::endl;
    exit (EXIT_FAILURE);
  } 
  t_GSA GSA;

  //Open BCR files
  FILE *InFileBWT = fopen(fnBWT.c_str(), "rb");
  if (InFileBWT==NULL) {
    std::cerr << "Entire BWT file: Error opening "  << fnBWT <<  "." << std::endl;
    exit (EXIT_FAILURE);
  }


  FILE *InFileLCP = fopen(fnLCP.c_str(), "rb");
  if (InFileLCP==NULL) {
    std::cerr << "Entire LCP file: Error opening "  << fnLCP << "." <<std::endl;
    exit (EXIT_FAILURE);
  }     

  FILE *InFilePairSA = fopen(fnPairSA.c_str(), "rb");
  if (InFilePairSA==NULL) {
    std::cerr << "Entire Pairs SA file: Error opening " << fnPairSA <<  "." << std::endl;
    exit (EXIT_FAILURE);
  }

  dataTypeNChar numcharBWT, numcharPairSA, numcharLCP, numEle=0;

  uchar *bufferBWT = new uchar[SIZEBUFFER];
  ElementType *buffer = new ElementType[SIZEBUFFER];
  dataTypelenSeq *bufferLCP = new dataTypelenSeq[SIZEBUFFER];
        
  while ((!feof(InFileBWT))  && (!feof(InFilePairSA)) && (!feof(InFileLCP)) )  {
    numcharBWT = fread(bufferBWT,sizeof(uchar),SIZEBUFFER,InFileBWT);
    numcharPairSA = fread(buffer,sizeof(ElementType),SIZEBUFFER,InFilePairSA);
    numcharLCP = fread(bufferLCP,sizeof(dataTypelenSeq),SIZEBUFFER,InFileLCP);
    if ((numcharLCP != numcharPairSA) || (numcharLCP != numcharBWT))
      std::cerr << "Error: number in BWT in Pair SA  and in LCP\n";
    else {
      //std::cerr << "bwt\tlcp\tpos\tnumSeq\tSA\tQS\n";
      for (dataTypeNChar i=0; i < numcharBWT; i++) {
        numEle++;

        GSA.suff = buffer[i].sa;
        GSA.text = buffer[i].numSeq;
        
        GSA.lcp = bufferLCP[i];

        if (bufferBWT[i] == '$')
          GSA.bwt='\0';
        else
          GSA.bwt=bufferBWT[i];


        fwrite(&GSA.text, sizeof(dataTypeNSeq), 1, f_ESA);  
        fwrite(&GSA.suff, sizeof(dataTypelenSeq), 1, f_ESA);  
        fwrite(&GSA.lcp, sizeof(dataTypelenSeq), 1, f_ESA); 
        fwrite(&GSA.bwt, sizeof(uchar), 1, f_ESA);  


      }
    }
  }
  
  std::cerr <<  "The total number of elements is " << numEle << "\n\n";
  delete [] buffer; 
  delete [] bufferBWT;
  delete [] bufferLCP; 
  fclose(InFilePairSA);
  fclose(InFileLCP);
  fclose(InFileBWT);
  fclose(f_ESA);
  
  return 0;
}


#endif
  
void BCRexternalBWT::storeEntireBWTFilePartial( const char* fn ) {

  static FILE *OutFileBWT, *InFileBWT;                  // output and input file BWT;
  char *filenameIn = new char[12];
  char *filename = new char[8];
  const char *ext = ".aux";


  dataTypeNChar numchar=0;
  dataTypeNChar numcharWrite=0;

  uchar *buffer = new uchar[SIZEBUFFER];

  dataTypeNChar *freqOut = new dataTypeNChar [SIZE_ALPHA];
  for (dataTypedimAlpha i = 0; i < SIZE_ALPHA-1; ++i)
        freqOut[i] = 0;
  freqOut[SIZE_ALPHA-1] = 0;

  OutFileBWT = fopen(fn, "wb");
  if (OutFileBWT==NULL) {
    std::cerr << "Entire BWT file: Error opening " << std::endl;
    exit (EXIT_FAILURE);
  }



  std::cerr << "Build the entire BWT file and compute the distribution of chars.\n";

  for (dataTypedimAlpha g = 0 ; g < sizeAlpha; g++) {
    numchar=sprintf (filename, "bwt_%d", g);
    numchar=sprintf (filenameIn,"%s%s",filename,ext);
    InFileBWT = fopen(filenameIn, "rb");
    if (InFileBWT==NULL) {
      std::cerr << "BWT file " << (unsigned int)g <<": Error opening " << std::endl;
      exit (EXIT_FAILURE);
    }


    //std::cerr << "BWT file " << (unsigned int)g << "= ";
    numchar = fread(buffer,sizeof(uchar),SIZEBUFFER,InFileBWT);
    numcharWrite = fwrite (buffer, sizeof(uchar), numchar , OutFileBWT);
    //std::cerr << "numchar= " << numchar << std::endl;
    //std::cerr << "numcharWrite= " << numcharWrite << std::endl;
    assert(numchar == numcharWrite); // we should always read/write the same number of characters

    for (dataTypeNChar j = 0 ; j < numchar; j++)
          freqOut[(unsigned int)(buffer[j])]++;



    while (numchar!=0) {
      numchar = fread(buffer,sizeof(uchar),SIZEBUFFER,InFileBWT);
      numcharWrite = fwrite (buffer, sizeof(uchar), numchar , OutFileBWT);
      assert(numchar == numcharWrite); // we should always read/write the same number of characters

      for (dataTypeNChar j = 0 ; j < numchar; j++)
        freqOut[(unsigned int)(buffer[j])]++;


    }

    fclose(InFileBWT);


    //Delete or rename the partial BWT bwt_%d????
    #if deletePartialBWT == 1
      if (remove(filenameIn)!=0)
        std::cerr << "storeEntireBWTFilePartial: Error deleting file" << std::endl;

      /*
    #else //rename the aux bwt file
      int lung = strlen(fn) + strlen(filenameIn)+1;
      char *newfilename = new char[lung];
      numchar=sprintf (newfilename,"%s%s",fn,filenameIn);
      //std::cerr  << newfilename << " " << filenameIn << std::endl;
      if(rename(filenameIn, newfilename))
          std::cerr  <<"storeEntireBWTFilePartial: Error renaming file" << std::endl;
      delete[] newfilename;
    */
    #endif
  }
  fclose(OutFileBWT);

  #if verboseEncode==1
    std::cerr << "\nThe Entire BWT:"<< std::endl;
    OutFileBWT = fopen(fn, "rb");
    if (OutFileBWT==NULL) {
      std::cerr << "Entire BWT file: Error opening " << std::endl;
      exit (EXIT_FAILURE);
    }

    numchar = 1;
    for (dataTypeNSeq g = 0 ; g < SIZEBUFFER; g++)
      buffer[g] = '\0';
    numchar = fread(buffer,sizeof(uchar),SIZEBUFFER,OutFileBWT);
      if (numchar==0)
         std::cerr  << "empty\n";
    else
          std::cerr  << buffer;
        while (numchar!=0) {
          for (dataTypeNChar g = 0 ; g < SIZEBUFFER; g++)
            buffer[g] = '\0';
         numchar = fread(buffer,sizeof(uchar),SIZEBUFFER,OutFileBWT);
         if (numchar!=0)
            std::cerr  << buffer;
        }
    std::cerr << std::endl;
    fclose(OutFileBWT);
  #endif

  delete [] buffer;

  int check=0;
  std::cerr << "Distribution in BWT\n";
  for (dataTypedimAlpha g = 0 ; g < SIZE_ALPHA-1; g++) {
    if (freqOut[g] > 0){
      InFileBWT = openFilePartialIn(alpha[(unsigned int)g]);
      fseek(InFileBWT,0,SEEK_END);
      dataTypeNChar lengthBWTPartial=ftell(InFileBWT);
      if (freqOut[g] == lengthBWTPartial)
        std::cerr << (unsigned int)g << " freq: " << freqOut[g] << "\n";
      else {
        std::cerr << (unsigned int)g << " freq: " << freqOut[g] << " freq in BWT: " << lengthBWTPartial << "***********************\n";
        check = 1;
      }
      closeFilePartial(InFileBWT);
    }
  }
  if (freqOut[(SIZE_ALPHA-1)] > 0){
    InFileBWT = openFilePartialIn(alpha[(unsigned int)(SIZE_ALPHA-1)]);
    fseek(InFileBWT,0,SEEK_END);
    dataTypeNChar lengthBWTPartial=ftell(InFileBWT);
    if (freqOut[(SIZE_ALPHA-1)] == lengthBWTPartial)
      std::cerr << (unsigned int)(SIZE_ALPHA-1) << " freq: " << freqOut[(SIZE_ALPHA-1)] << "\n";
    else {
      std::cerr << (unsigned int)(SIZE_ALPHA-1) << " freq: " << freqOut[(SIZE_ALPHA-1)] << " freq in BWT: " << lengthBWTPartial << "***********************\n";
      check = 1;
    }
    closeFilePartial(InFileBWT);
  }

  if (check == 1)
    std::cerr << "WARNING! Some length partial BWT != the number of occorrences! " << std::endl;

  delete [] freqOut;
  delete [] filenameIn;
  delete [] filename;


}

#if BUILD_LCP == 1
void BCRexternalBWT::storeEntireLCP( const char* fn ) {

  static FILE *OutFileLCP, *InFileLCP;                  // output and input file LCP;
  char *filenameIn = new char[12];
  char *filename = new char[8];
  const char *ext = ".aux";
  dataTypeNChar numchar=0;
  dataTypeNChar numcharWrite=0;

  dataTypelenSeq *buffer = new dataTypelenSeq[SIZEBUFFER];

  int lung = strlen(fn);
  char *fnLCP = new char[lung+7];
  numchar=sprintf (fnLCP,"%s%s",fn,".lcp");
  OutFileLCP = fopen(fnLCP, "wb");
  if (OutFileLCP==NULL) {
    std::cerr << "Entire LCP file: Error opening " << std::endl;
    exit (EXIT_FAILURE);
  }

  std::cerr << "Entire LCP file" << std::endl;
  dataTypeNChar numTotLcp = 0;

  for (dataTypedimAlpha g = 0 ; g < sizeAlpha; g++) {
    numchar=sprintf (filename, "lcp_%d", g);
    numchar=sprintf (filenameIn,"%s%s",filename,ext);
    InFileLCP = fopen(filenameIn, "rb");
    if (InFileLCP==NULL) {
      std::cerr << "LCP file " << (unsigned int)g <<": Error opening " << std::endl;
      exit (EXIT_FAILURE);
    }
    //std::cerr << "LCP file " << (unsigned int)g << "= ";
    numchar = fread(buffer,sizeof(dataTypelenSeq),SIZEBUFFER,InFileLCP);
    numcharWrite = fwrite (buffer, sizeof(dataTypelenSeq), numchar , OutFileLCP);
    assert(numchar == numcharWrite); // we should always read/write the same number of characters
    numTotLcp += numcharWrite;
    while (numchar!=0) {
      numchar = fread(buffer,sizeof(dataTypelenSeq),SIZEBUFFER,InFileLCP);
      numcharWrite = fwrite (buffer, sizeof(dataTypelenSeq), numchar , OutFileLCP);
      assert(numchar == numcharWrite); // we should always read/write the same number of characters
      numTotLcp += numcharWrite;
    }

    fclose(InFileLCP);


    //Delete or renome the partial LCP lcp_%d????
    #if (deletePartialLCP == 1)

      if (remove(filenameIn)!=0)
        std::cerr << "storeEntireLCP: Error deleting file" << std::endl;
    /*
    #else //renome the aux bwt file

      int lung = strlen(fn) + strlen(filenameIn)+1;
      char *newfilename = new char[lung];
      numchar=sprintf (newfilename,"%s%s",fn,filenameIn);
      //std::cerr  << newfilename << " " << filenameIn << std::endl;
      if(rename(filenameIn, newfilename))
          std::cerr  <<"storeEntireLCP: Error renaming file" << std::endl;
    */
    #endif
  }
  fclose(OutFileLCP);

  std::cerr << "LCP file contains " << numTotLcp << " values\n";

  #if verboseEncode==1
    std::cerr << "\nThe Entire LCP:"<< std::endl;
    OutFileLCP = fopen(fnLCP, "rb");
    if (OutFileLCP==NULL) {
      std::cerr << "Entire LCP file: Error opening " << std::endl;
      exit (EXIT_FAILURE);
    }
    for (dataTypeNChar g = 0 ; g < SIZEBUFFER; g++)
        buffer[g] = '\0';
    numchar = fread(buffer,sizeof(dataTypelenSeq),SIZEBUFFER,OutFileLCP);
    if (numchar==0)
      std::cerr  << "empty\n";
    else
      for (dataTypeNChar g = 0 ; g < numchar; g++)
        std::cerr  << (unsigned int)buffer[g] << " ";
    while (numchar!=0) {
      for (dataTypeNChar g = 0 ; g < SIZEBUFFER; g++)
        buffer[g] = '\0';
      numchar = fread(buffer,sizeof(dataTypelenSeq),SIZEBUFFER,OutFileLCP);
      if (numchar!=0)
        for (dataTypeNChar g = 0 ; g < numchar; g++)
          std::cerr  << (unsigned int)buffer[g] << " ";
    }
    std::cerr << std::endl;
    fclose(OutFileLCP);
  #endif

  delete [] buffer;
  delete [] filenameIn;
  delete [] filename;

  delete[] fnLCP;
}
#endif

#if BUILD_SA == 1
void BCRexternalBWT::storeEntirePairSA( const char* fn ) {

  static FILE *OutFileSA, *InFileSA;                  // output and input file SA;
  char *filenameIn = new char[13];
  char *filename = new char[9];
  const char *ext = ".aux";
  dataTypeNChar numcharWrite, numcharRead;
  ElementType *buffer = new ElementType[SIZEBUFFER];

  std::cerr << "\nEntire Pairs SA file (position, number of sequence)" << std::endl;

  int lung = strlen(fn);
  char *fnSA = new char[lung+8];
  numcharRead=sprintf (fnSA,"%s%s",fn,".pairSA");

  OutFileSA = fopen(fnSA, "wb");
  if (OutFileSA==NULL) {
    std::cerr << "Entire Pairs SA file: Error opening " << fnSA << std::endl;
    exit (EXIT_FAILURE);
  }


  dataTypeNSeq numTotalTexts;
  #if BUILD_BCR_FROM_EGSA == 1
    numTotalTexts = nText+nAddedTextEGSA;
  #else
    numTotalTexts = nText;
  #endif

  std::vector <dataTypelenSeq> vectLen;
  vectLen.resize(numTotalTexts);

  char fileLen[512];
  sprintf(fileLen, "%s.lenSeqs.aux", fn);
  //sprintf(fileLen, "outFileLen.aux");
  static FILE *InFileLen;                  // file of the lengths;
  InFileLen = fopen(fileLen, "rb");
  if (InFileLen==NULL) {
      std::cerr << "storeEntireSAfromPairSA: could not open file \"" << fileLen << "\"!"<< std::endl;
      exit (EXIT_FAILURE);
  }

  numcharRead = fread (&vectLen[0], sizeof(dataTypelenSeq), vectLen.size() , InFileLen);
  assert(numcharRead == numTotalTexts); // we should always read the same number of characters

  fclose(InFileLen);

  dataTypeNChar numTotPairSA = 0;


  for (dataTypedimAlpha g = 0 ; g < sizeAlpha; g++) {
    numcharRead=sprintf (filename, "sa_%d", g);
    numcharRead=sprintf (filenameIn,"%s%s",filename,ext);
    InFileSA = fopen(filenameIn, "rb");
    if (InFileSA==NULL) {
      std::cerr << "SA file " << (unsigned int)g <<": Error opening " << std::endl;
      exit (EXIT_FAILURE);
    }

    numcharRead = fread(buffer,sizeof(ElementType),SIZEBUFFER,InFileSA);
    numcharWrite = fwrite (buffer, sizeof(ElementType), numcharRead , OutFileSA);
    assert (numcharRead == numcharWrite);
    numTotPairSA += numcharWrite;

    while (numcharRead!=0) {
      numcharRead = fread(buffer,sizeof(ElementType),SIZEBUFFER,InFileSA);

      numcharWrite = fwrite (buffer, sizeof(ElementType), numcharRead , OutFileSA);
      assert (numcharRead == numcharWrite);
      numTotPairSA += numcharWrite;
    }

    fclose(InFileSA);

    //if (remove(filenameIn)!=0)
    //  std::cerr << filenameIn <<": Error deleting file" << std::endl;
    //Delete or rename the partial SA sa_%d????
    #if deletePartialSA == 1
      if (remove(filenameIn)!=0)
        std::cerr << "storeEntirePairSA: Error deleting file" << std::endl;
    /*#else //rename the aux sa file
      int lung = strlen(fn) + strlen(filenameIn)+1;
      char *newfilename = new char[lung];
      sprintf (newfilename,"%s%s",fn,filenameIn);
      //std::cerr  << newfilename << " " << filenameIn << std::endl;
      if(rename(filenameIn, newfilename))
          std::cerr  <<"storeEntirePairSA: Error renaming file" << std::endl;
      delete[] newfilename;
      //------------------
    */
    #endif
    //-------------------------------
  }
  fclose(OutFileSA);
  #if deletePartialSA == 0
    #if verboseEncode == 1
      std::cerr << "Final Generalized Suffix array in segments: " << std::endl;
      dataTypeNChar numcharSA=0;
      char *filenameInSA = new char[12];
      char *filenameSA = new char[8];
      ElementType *bufferGSA = new ElementType[SIZEBUFFER];
      dataTypeNChar mmm=0;
      while (mmm < sizeAlpha) {
        numcharSA=sprintf (filenameSA, "sa_%d", mmm);
        numcharSA=sprintf (filenameInSA,"%s%s",filenameSA,ext);
        InFileSA = fopen(filenameInSA, "rb");
        if (InFileSA==NULL) {
          std::cerr << "In SA file " << (unsigned int)mmm <<": Error opening: " << filenameInSA << std::endl;
          exit (EXIT_FAILURE);
        }
        for (dataTypeNChar g = 0 ; g < SIZEBUFFER; g++) {
          bufferGSA[g].sa = 0;
          bufferGSA[g].numSeq  = 0;
        }
        numcharSA = fread(bufferGSA,sizeof(ElementType),SIZEBUFFER,InFileSA);
        std::cerr << "SA[" << (unsigned int)mmm << "]:\t";
        if (numcharSA==0)
          std::cerr  << "empty";
        else {
          for (dataTypeNChar g = 0 ; g < numcharSA; g++)
            std::cerr << "(" << bufferGSA[g].sa  << " " << bufferGSA[g].numSeq << ")\t";
        }
        while (numcharSA!=0) {
          for (dataTypeNChar g = 0 ; g < SIZEBUFFER; g++) {
            bufferGSA[g].sa = 0;
            bufferGSA[g].numSeq  = 0;
          }
          numcharSA = fread(bufferGSA,sizeof(ElementType),SIZEBUFFER,InFileSA);
          if (numcharSA!=0)
            for (dataTypeNChar g = 0 ; g < numcharSA; g++)
              std::cerr << "(" << bufferGSA[g].sa  << " " << bufferGSA[g].numSeq << ")\t";
        }
        std::cerr << std::endl;
        fclose(InFileSA);
        mmm++;
      }
      delete [] filenameInSA;
      delete [] filenameSA;
      delete [] bufferGSA;
    #endif
  #endif

  #if verboseEncode==1
    std::cerr << "Final Generalized Suffix array: " << std::endl;
    OutFileSA = fopen(fnSA, "rb");
    if (OutFileSA==NULL) {
      std::cerr << "Entire SA file: Error opening (in VerboseEncode) " << fnSA << std::endl;
      exit (EXIT_FAILURE);
    }

    numcharRead = fread(buffer,sizeof(ElementType),SIZEBUFFER,OutFileSA);
    if (numcharRead==0)
      std::cerr  << "empty\n";
    else
      for (dataTypeNChar g = 0 ; g < numcharRead; g++) {
        std::cerr  << "(" << (unsigned int)buffer[g].sa << "," << buffer[g].numSeq << ") ";
      }
    while (numcharRead!=0) {
      numcharRead = fread(buffer,sizeof(ElementType),SIZEBUFFER,OutFileSA);
      for (dataTypeNChar g = 0 ; g < numcharRead; g++) {
          std::cerr  << "(" << buffer[g].sa << "," << buffer[g].numSeq << ") ";
      }
    }
    std::cerr << std::endl;

    fclose(OutFileSA);
  #endif

  std::cerr << "Pair SA file contains " << numTotPairSA << " values\n";

  delete [] buffer;

  delete [] filenameIn;
  delete [] filename;

  //delete[] fileLen;
  delete[] fnSA;
  //------------------
}
#endif

#if BUILD_EXT_MEMORY==0
void BCRexternalBWT::storeEntireBWTIntMem( const char* fn ) {
  static FILE *OutFileBWT;                  // output and input file BWT;
  dataTypeNChar numcharWrite=0;

  dataTypeNChar *freqOut = new dataTypeNChar [SIZE_ALPHA];
  for (dataTypedimAlpha i = 0; i < SIZE_ALPHA-1; ++i)
        freqOut[i] = 0;
  freqOut[SIZE_ALPHA-1] = 0;

  OutFileBWT = fopen(fn, "wb");
  if (OutFileBWT==NULL) {
    std::cerr << "Entire BWT file: Error opening " << std::endl;
    exit (EXIT_FAILURE);
  }

  std::cerr << "Build the entire BWT file and compute the distribution of chars.\n";
  std::cerr << "Distribution in BWT\n";
  for (dataTypedimAlpha g = 0 ; g < sizeAlpha; g++) {
      std::cerr << "BWT segment: " << (unsigned int)g << " size " << vectVectBWT[g].size() << std::endl;
      for (dataTypeNChar mm = 0 ; mm < vectVectBWT[g].size(); mm++) {
        numcharWrite = fwrite (&vectVectBWT[g][mm], sizeof(uchar), 1 , OutFileBWT);
      }
  }

  fclose(OutFileBWT);
  int check = 0;
  for (dataTypedimAlpha g = 0 ; g < SIZE_ALPHA-1; g++) {
    if (freqOut[g] > 0) {
      std::cerr << "BWT segment: " << (unsigned int)alpha[(unsigned int)g] << " size " << vectVectBWT[(unsigned int)alpha[(unsigned int)g]].size() << std::endl;
      if (freqOut[g] == vectVectBWT[alpha[(unsigned int)g].size())
        std::cerr << (unsigned int)g << " freq: " << freqOut[g] << "\n";
      else {
        std::cerr << (unsigned int)g << " freq: " << freqOut[g] << " freq in BWT: " << vectVectBWT[(unsigned int)alpha[(unsigned int)g]].size() << "***********************\n";
        check = 1;
      }
    }
  }

  if (check == 1)
    std::cerr << "WARNING! Some length partial BWT != the number of occorrences! " << std::endl;

  #if verboseEncode==1
    uchar *buffer = new uchar[SIZEBUFFER];
    std::cerr << "\nThe Entire BWT:"<< std::endl;
    OutFileBWT = fopen(fn, "rb");
    if (OutFileBWT==NULL) {
      std::cerr << "Entire BWT file: Error opening " << std::endl;
      exit (EXIT_FAILURE);
    }

    numcharWrite = 1;
    for (dataTypeNSeq g = 0 ; g < SIZEBUFFER; g++)
      buffer[g] = '\0';
    numcharWrite = fread(buffer,sizeof(uchar),SIZEBUFFER,OutFileBWT);
      if (numcharWrite==0)
         std::cerr  << "empty\n";
    else
          std::cerr  << buffer;
        while (numcharWrite!=0) {
          for (dataTypeNChar g = 0 ; g < SIZEBUFFER; g++)
            buffer[g] = '\0';
         numcharWrite = fread(buffer,sizeof(uchar),SIZEBUFFER,OutFileBWT);
         if (numcharWrite!=0)
            std::cerr  << buffer;
        }
    std::cerr << std::endl;
    fclose(OutFileBWT);
    delete [] buffer;
  #endif



  delete [] freqOut;
}
#endif

#if BUILD_SA == 1
void BCRexternalBWT::storeEntireSAfromPairSA( const char* fn ) {
  static FILE *OutFileSA, *InFilePairSA;                  // output and input file SA;

  std::cerr << "\nSA file from pair SA file" << std::endl;

  dataTypeNChar numchar, numcharWrite;

  std::vector <dataTypelenSeq> vectSumCumLen;

  char fileLen[512];
  sprintf(fileLen, "%s.lenSeqs.aux", fn);
  //sprintf(fileLen, "outFileLen.aux");
  static FILE *InFileLen;                  // file of the lengths;
  InFileLen = fopen(fileLen, "rb");
  if (InFileLen==NULL) {
      std::cerr << "storeEntireSAfromPairSA: could not open file \"" << fileLen << "\"!"<< std::endl;
      exit (EXIT_FAILURE);
  }
  dataTypelenSeq lenSeq=0;
  numchar = fread (&lenSeq, sizeof(dataTypelenSeq), 1 , InFileLen);
  assert( numchar == 1); // we should always read the same number of characters


  dataTypeNSeq numTotalTexts;
  #if BUILD_BCR_FROM_EGSA == 1
    numTotalTexts = nText+nAddedTextEGSA;
  #else
    numTotalTexts = nText;
  #endif
  vectSumCumLen.resize(numTotalTexts+1);
  vectSumCumLen[0] = 0;
  vectSumCumLen[1] = lenSeq + 1;   //Plus $

  std::cerr << "Total number of sequences = " << numTotalTexts << std::endl;
  //#if verboseEncode == 1
  //  std::cerr << "storeEntireSAfromPairSA: seq 1 lenSeq= " << lenSeq << std::endl;
  //#endif
  for (dataTypeNSeq num = 2; num < numTotalTexts+1; num++) {
    numchar = fread (&lenSeq, sizeof(dataTypelenSeq), 1 , InFileLen);
    assert(numchar == 1); // we should always read the same number of characters
    //#if verboseEncode == 1
    //  std::cerr << "storeEntireSAfromPairSA: seq " << num <<  " lenSeq= " << lenSeq << std::endl;
    //#endif
    vectSumCumLen[num] = vectSumCumLen[num-1] + lenSeq + 1;  //Plus $
    //std::cerr << "storeEntireSAfromPairSA: seq " << num <<  " vectSumCumLen[num]= " << vectSumCumLen[num] << std::endl;
  }
  fclose(InFileLen);

  int lung = strlen(fn);
  char *fnEntireSA = new char[lung+4];
  char *fnPairSA = new char[lung+8];
  numchar=sprintf (fnEntireSA,"%s%s",fn,".sa");
  numchar=sprintf (fnPairSA,"%s%s",fn,".pairSA");

  InFilePairSA = fopen(fnPairSA, "rb");
  if (InFilePairSA==NULL) {
    std::cerr << "storeEntireSAfromPairSA: Entire Pairs SA file: Error opening " << fnPairSA << std::endl;
    exit (EXIT_FAILURE);
  }

  OutFileSA = fopen(fnEntireSA, "wb");
  if (OutFileSA==NULL) {
    std::cerr << "storeEntireSAfromPairSA: Entire SA file: Error opening " << fnEntireSA << std::endl;
    exit (EXIT_FAILURE);
  }

  ElementType *buffer = new ElementType[SIZEBUFFER];
  dataTypeNChar *bufferNChar = new dataTypeNChar[SIZEBUFFER];

  while (!feof(InFilePairSA)) {
    numchar = fread(buffer,sizeof(ElementType),SIZEBUFFER,InFilePairSA);
    //std::cerr << "number read " << numchar << "\n";
    if (numchar > 0) {
      for (dataTypeNChar i=0; i < numchar; i++) {
        bufferNChar[i] = (dataTypeNChar)(vectSumCumLen[buffer[i].numSeq] + buffer[i].sa);                                         //buffer[i].numSeq * (lengthRead+1) + buffer[i].sa;
        //std::cerr << "vectSumCumLen["<< buffer[i].numSeq<< "]= " << (unsigned int)vectSumCumLen[buffer[i].numSeq] << " + " << (unsigned int)buffer[i].sa << " --> " << (unsigned int)bufferNChar[i] << "\n";
        //std::cerr << buffer[i].numSeq << " " << (unsigned int)lengthRead << " " << (unsigned int)buffer[i].sa << " --> " << (unsigned int)bufferNChar[i] << "\n";
      }
      numcharWrite = fwrite (bufferNChar, sizeof(dataTypeNChar), numchar, OutFileSA);
      assert(numcharWrite == numchar);    //added
      //std::cerr << "number write " << numcharWrite << "\n";
    }
  }
  fclose(InFilePairSA);
  fclose(OutFileSA);

  #if verboseEncode==1
    std::cerr << "\nThe Entire SA. The file is "<< fnEntireSA << std::endl;
    OutFileSA = fopen(fnEntireSA, "rb");
    if (OutFileSA==NULL) {
      std::cerr << "storeEntireSAfromPairSA: Entire SA file (in verboseEncode): Error opening " << fnEntireSA << std::endl;
      exit (EXIT_FAILURE);
    }

    numchar = fread(bufferNChar,sizeof(dataTypeNChar),SIZEBUFFER,OutFileSA);
    if (numchar==0)
      std::cerr  << "empty\n";
    else
      for (dataTypeNChar g = 0 ; g < numchar; g++) {
        std::cerr  << bufferNChar[g] << " ";
      }
    while (numchar!=0) {
      numchar = fread(buffer,sizeof(ElementType),SIZEBUFFER,OutFileSA);
      for (dataTypeNChar g = 0 ; g < numchar; g++) {
        std::cerr  << bufferNChar[g] << " ";
      }
    }
    std::cerr << std::endl;

    fclose(OutFileSA);
  #endif

  std::cerr << std::endl;
  delete [] buffer;
  delete [] bufferNChar;
  delete [] fnEntireSA;
  delete [] fnPairSA;

}
#endif


void BCRexternalBWT::printSegments()
{
  dataTypeNChar numchar;
  const char *ext = ".aux";
  dataTypedimAlpha mmm=0;
  uchar *buffer = new uchar[SIZEBUFFER];

  std::cerr  << "\nPartial BWT segments:\n";

  #if BUILD_EXT_MEMORY==1
    static FILE *InFileBWT;
    char *filenameIn = new char[12];
    char *filename = new char[8];
    while (mmm < sizeAlpha) {
      numchar=sprintf (filename, "bwt_%d", mmm);
      numchar=sprintf (filenameIn,"%s%s",filename,ext);
      //printf("===currentPile= %d\n",mmm);
      InFileBWT = fopen(filenameIn, "rb");
      for (dataTypeNChar g = 0 ; g < SIZEBUFFER; g++)
        buffer[g] = '\0';
      numchar = fread(buffer,sizeof(uchar),SIZEBUFFER,InFileBWT);
      std::cerr << "B[" << (unsigned int)mmm << "]:\t";
      if (numchar==0)
        std::cerr  << "empty\n";
      else
        std::cerr  << buffer << "\n";
      fclose(InFileBWT);
      mmm++;
    }
    delete [] filenameIn;
    delete [] filename;
  #else
    for (mmm=0; mmm<sizeAlpha; mmm++) {
      std::cerr << "B[" << (unsigned int)mmm << "]:\t";
      if (vectVectBWT[mmm].size() == 0)
        std::cerr  << "empty\n";
      else {
        for (dataTypeNChar r=0; r<vectVectBWT[mmm].size(); r++) {
          std::cerr  << vectVectBWT[mmm][r];
        }
        std::cerr  << "\n";
      }
    }
  #endif


  delete [] buffer;

    #if BUILD_LCP == 1
    static FILE *InFileLCP;                  // output and input file LCP;
    char *filenameInLCP = new char[12];
    char *filenameLCP = new char[8];
    //const char *ext = ".aux";
    dataTypelenSeq *bufferLCP = new dataTypelenSeq[SIZEBUFFER];
    std::cerr << "\nPartial LCP array: " << std::endl;
    mmm=0;
    while (mmm < sizeAlpha) {
      numchar=sprintf (filenameLCP, "lcp_%d", mmm);
      numchar=sprintf (filenameInLCP,"%s%s",filenameLCP,ext);
      //printf("===currentPile= %d\n",mmm);
      InFileLCP = fopen(filenameInLCP, "rb");
      for (dataTypeNChar g = 0 ; g < SIZEBUFFER; g++)
        bufferLCP[g] = 0;
      numchar = fread(bufferLCP,sizeof(dataTypelenSeq),SIZEBUFFER,InFileLCP);
      std::cerr << "L[" << (unsigned int)mmm << "]:\t";
      if (numchar==0)
        std::cerr  << "empty";
      else
        for (dataTypeNSeq g = 0 ; g < numchar; g++)
          std::cerr  << (unsigned int)bufferLCP[g] << " ";
      std::cerr  << "\n";
      fclose(InFileLCP);
      mmm++;
    }
    delete [] bufferLCP;
    delete [] filenameInLCP;
    delete [] filenameLCP;
  #endif

  #if BUILD_SA == 1
    std::cerr << "\nPartial Generalized Suffix array: " << std::endl;
    dataTypeNChar numcharSA=0;
    static FILE *InFileSA;                  // output and input file SA;
    char *filenameInSA = new char[12];
    char *filenameSA = new char[8];
    ElementType *bufferGSA = new ElementType[SIZEBUFFER];
    mmm=0;
    while (mmm < sizeAlpha) {
      numcharSA=sprintf (filenameSA, "sa_%d", mmm);
      numcharSA=sprintf (filenameInSA,"%s%s",filenameSA,ext);
      InFileSA = fopen(filenameInSA, "rb");
      if (InFileSA==NULL) {
        std::cerr << "In SA file " << (unsigned int)mmm <<": Error opening: " << filenameInSA << std::endl;
        exit (EXIT_FAILURE);
      }
      for (dataTypeNChar g = 0 ; g < SIZEBUFFER; g++) {
        bufferGSA[g].sa = 0;
        bufferGSA[g].numSeq  = 0;
      }
      numcharSA = fread(bufferGSA,sizeof(ElementType),SIZEBUFFER,InFileSA);
      std::cerr << "SA[" << (unsigned int)mmm << "]:\t";
      if (numcharSA==0)
        std::cerr  << "empty";
      else {
        for (dataTypeNChar g = 0 ; g < numcharSA; g++)
          std::cerr << "(" << bufferGSA[g].sa  << " " << bufferGSA[g].numSeq << ")\t";
      }
      while (numcharSA!=0) {
        for (dataTypeNChar g = 0 ; g < SIZEBUFFER; g++) {
          bufferGSA[g].sa = 0;
          bufferGSA[g].numSeq  = 0;
        }
        numcharSA = fread(bufferGSA,sizeof(ElementType),SIZEBUFFER,InFileSA);
        if (numcharSA!=0)
          for (dataTypeNChar g = 0 ; g < numcharSA; g++)
            std::cerr << "(" << bufferGSA[g].sa  << " " << bufferGSA[g].numSeq << ")\t";
      }
      std::cerr << std::endl;
      fclose(InFileSA);
      mmm++;
    }
    delete [] filenameInSA;
    delete [] filenameSA;
    delete [] bufferGSA;
  #endif
}

void BCRexternalBWT::printOutput(char *fileOutput)
{
  if ((BUILD_LCP==1) && (BUILD_SA==1) && (OUTPUT_FORMAT_EGSA)==0) {
    std::cerr << "Store the text file containing the BWT, LCP and SA\n";

    int lung = strlen(fileOutput);
    char *fnSA = new char[lung+4];
    char *fnPairSA = new char[lung+8];
    char *fileOutRes = new char[lung+5];
    char *fnLCP = new char[lung+8];
    sprintf (fnLCP,"%s%s",fileOutput,".lcp");
    sprintf (fnSA,"%s%s",fileOutput,".sa");
    sprintf (fnPairSA,"%s%s",fileOutput,".pairSA");
    sprintf (fileOutRes,"%s%s",fileOutput,".txt");

    std::cerr << "BCRexternalBWT: fileoutput: "  << fileOutput <<  "." << std::endl;
    std::cerr << "BCRexternalBWT: fnLCP: "  << fnLCP <<  "." << std::endl;
    std::cerr << "BCRexternalBWT: fnPairSA: "  << fnPairSA <<  "." << std::endl;
    std::cerr << "BCRexternalBWT: fileOutRes: "  << fileOutRes <<  "." << std::endl;

    FILE *InFileBWT = fopen(fileOutput, "rb");
    if (InFileBWT==NULL) {
      std::cerr << "BCRexternalBWT: Entire BWT file: Error opening "  << fileOutput <<  "." << std::endl;
      exit (EXIT_FAILURE);
    }


    FILE *InFileLCP = fopen(fnLCP, "rb");
    if (InFileLCP==NULL) {
      std::cerr << "BCRexternalBWT: Entire LCP file: Error opening "  << fnLCP << "." <<std::endl;
      exit (EXIT_FAILURE);
    }

    FILE *InFilePairSA = fopen(fnPairSA, "rb");
    if (InFilePairSA==NULL) {
      std::cerr << "BCRexternalBWT: Entire Pairs SA file: Error opening " << fnPairSA <<  "." << std::endl;
      exit (EXIT_FAILURE);
    }

    FILE* InFileSA = fopen(fnSA, "rb");
    if (InFileSA==NULL) {
      std::cerr << "BCRexternalBWT: Entire SA file: Error opening " << fnSA << std::endl;
      exit (EXIT_FAILURE);
    }

    FILE *OutFile = fopen(fileOutRes, "w");
    if (OutFile==NULL) {
      std::cerr << "BCRexternalBWT: Error opening output \"" << fileOutRes << "\" file"<< std::endl;
      exit (EXIT_FAILURE);
    }

    dataTypeNChar numcharBWT, numcharPairSA, numcharSA, numcharLCP;

    uchar *bufferBWT = new uchar[SIZEBUFFER];
    ElementType *buffer = new ElementType[SIZEBUFFER];
    dataTypeNChar *bufferNChar = new dataTypeNChar[SIZEBUFFER];
    dataTypelenSeq *bufferLCP = new dataTypelenSeq[SIZEBUFFER];

    while ((!feof(InFileBWT)) && (!feof(InFileSA)) && (!feof(InFilePairSA)) && (!feof(InFileLCP)) )  {
      numcharBWT = fread(bufferBWT,sizeof(uchar),SIZEBUFFER,InFileBWT);
      numcharPairSA = fread(buffer,sizeof(ElementType),SIZEBUFFER,InFilePairSA);
      numcharSA = fread(bufferNChar,sizeof(dataTypeNChar),SIZEBUFFER,InFileSA);
      numcharLCP = fread(bufferLCP,sizeof(dataTypelenSeq),SIZEBUFFER,InFileLCP);

      fprintf(OutFile, "bwt\tlcp\tpos\tnumSeq\tSA\r\n");


      if ((numcharPairSA != numcharSA) || (numcharLCP != numcharSA) || (numcharLCP != numcharBWT))
        std::cerr << "Error: number  in BWT in Pair SA in SA and in LCP\n";
      else {
        //std::cerr << "bwt\tlcp\tpos\tnumSeq\tSA\n";
        for (dataTypeNChar i=0; i < numcharBWT; i++) {
          std::cerr << (unsigned int)buffer[i].sa << "\t"<< buffer[i].numSeq << "\t" << bufferNChar[i] << "\n";

          fprintf(OutFile, "%c\t%d\t%d\t%d\t%lu\n", bufferBWT[i], bufferLCP[i],buffer[i].sa, buffer[i].numSeq, bufferNChar[i]);
        }
      }
    }
    delete[] fnSA;
    delete[] fnPairSA;
    delete[] fileOutRes;
    delete[] buffer;
    delete[] bufferNChar;
    delete[] fnLCP;
    delete[] bufferBWT;
    delete[] bufferLCP;
    fclose(InFileBWT);
    fclose(InFilePairSA);
    fclose(InFileSA);
    fclose(OutFile);
    fclose(InFileLCP);
    //controllare
    delete[] alphaInverse;
    //------------------


  }
}

#if (BUILD_BCR_FROM_EGSA == 1)
dataTypeNChar BCRexternalBWT::readEGSA(char const * fileinput)
{
    #if (verboseEncode == 1)
      std::cerr << "TableOcc after the first insertions \n";
      for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
        for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++)
          std::cerr << tableOcc[j][h] << "\t";
        std::cerr << "\n";
      }
    #endif


    std::cerr << "readEGSA for building BCR. Now we append the output of EGSA. The number of (old) sequences is: " << nAddedTextEGSA << "\n";

    if ((BUILD_LCP == 0) || (BUILD_SA== 0)) {
      std::cerr << "We do not append the output of EGSA, because (BUILD_LCP == 0) or (BUILD_SA== 0) \n";
      exit (EXIT_FAILURE);
    }
    if (sizeAlpha > 6) {
        std::cerr << "readEGSA - BWT - Error sizeAlpha > 6: " << sizeAlpha << std::endl;
      }
    assert(sizeAlpha <= 6);


    char fnOutOcc [512];
    sprintf(fnOutOcc, "%s.occurrences.aux", fileinput);

    FILE* OutFileOcc = fopen(fnOutOcc, "rb");
    if (OutFileOcc==NULL) {
      std::cerr << "readEGSA - occurrences - Error opening file > 6: " << fileinput << std::endl;
      exit (EXIT_FAILURE);
    }

    unsigned long countChar[6];
    fread (countChar, sizeof(unsigned long), 6, OutFileOcc);

    #if (verboseEncode == 1)
      printf( "%lu\t", countChar[0]);
      printf( "%lu\t", countChar[1]);
      printf( "%lu\t", countChar[2]);
      printf( "%lu\t", countChar[3]);
      printf( "%lu\t", countChar[4]);
      printf( "%lu\n", countChar[5]);
    #endif

    std::cerr << "Number of added symbols (EGSA): " ;
    std::cerr << countChar[0]+countChar[1]+countChar[2]+countChar[3]+countChar[4]+countChar[5] << "\n";

    fclose(OutFileOcc);

    dataTypeNChar totSymbols=0, totSymbolsInBWT=0;

    static FILE *OutFileBWT;                  // output  file BWT;
    char *filenameOutBWT = new char[12];
    char *filenameBWT = new char[8];

    static FILE *OutFileSA;                  // output  file GSA;
    char *filenameOutSA = new char[12];
    char *filenameSA = new char[8];

    static FILE *OutFileLCP;                  // output  file LCP;
    char *filenameOutLCP = new char[12];
    char *filenameLCP = new char[8];



    const char *ext = ".aux";
    ElementType newEle;

    FILE    *f_ESA;     // pointer to the ESA input file
    char c_aux[512];
    //sprintf(c_aux, "7seqs.fasta.7.gesa");
    sprintf(c_aux, "%s.%d.gesa", fileinput, nAddedTextEGSA);
    std::cerr << "WARNING!!! readEGSA: File: " << c_aux << std::endl;

    f_ESA = fopen(c_aux, "rb");//rb
    if (f_ESA == NULL) {
      std::cerr << "readEGSA: Error opening: " << c_aux << std::endl;
      exit (EXIT_FAILURE);
    }
    fseek(f_ESA , 0, SEEK_SET);
    t_GSA GSA;

/*    uchar *bufferBWTegsa = new uchar[SIZEBUFFER];
    ElementType *bufferPAIRSAegsa = new ElementType[SIZEBUFFER];
    dataTypelenSeq *bufferLCPegsa = new dataTypelenSeq[SIZEBUFFER]; */

    for (dataTypedimAlpha i = 0; i < sizeAlpha; i++) {
      sprintf (filenameBWT, "bwt_%d", i);
      sprintf (filenameOutBWT,"%s%s",filenameBWT,ext);
      OutFileBWT = fopen(filenameOutBWT, "a+b");
      if (OutFileBWT==NULL) {
        std::cerr << "readEGSA: BWT file " << (unsigned int)i <<" : Error opening " << std::endl;
        exit (EXIT_FAILURE);
      }

      sprintf (filenameSA, "sa_%d",i);
      sprintf (filenameOutSA,"%s%s",filenameSA,ext);
      OutFileSA = fopen(filenameOutSA, "a+b");
      if (OutFileSA==NULL) {
        std::cerr << "readEGSA: SA file " << (unsigned int)i <<" : Error opening " << std::endl;
        exit (EXIT_FAILURE);
      }

      sprintf (filenameLCP, "lcp_%d",i);
      sprintf (filenameOutLCP,"%s%s",filenameLCP,ext);
      OutFileLCP = fopen(filenameOutLCP, "a+b");
      if (OutFileLCP==NULL) {
        std::cerr << "readEGSA: LCP file " << (unsigned int)i <<" : Error opening " << std::endl;
        exit (EXIT_FAILURE);
      }


      dataTypeNChar m=0;
      uchar c;
      //std::cerr << "i="<< (unsigned int)i << " e countChar[i]= " << countChar[i] << std::endl;
      for (dataTypeNChar h = 0 ; h < countChar[i]; h++) {
        fread(&GSA.text, sizeof(dataTypeNSeq), 1, f_ESA);
        fread(&GSA.suff, sizeof(dataTypelenSeq), 1, f_ESA);
        fread(&GSA.lcp, sizeof(dataTypelenSeq), 1, f_ESA);
        fread(&GSA.bwt, sizeof(uchar), 1, f_ESA);

        newEle.sa=GSA.suff;
        newEle.numSeq= GSA.text + nText;
        //cerr <<  " text=" << GSA.text << " sa=" << GSA.suff <<  " " << GSA.lcp ;
        if (GSA.bwt == '\0')
          c='$';
        else
          c=GSA.bwt;
        //cerr << " bwt=" << c << endl;
        fwrite(&newEle, sizeof(ElementType), 1, OutFileSA);
        fwrite(&GSA.lcp, sizeof(dataTypelenSeq), 1, OutFileLCP);
        fwrite(&c, sizeof(uchar), 1, OutFileBWT);


        tableOcc[i][alpha[(unsigned int)c]]++;
        totSymbolsInBWT++;
      }

      fclose(OutFileBWT);
      fclose(OutFileSA);
      fclose(OutFileLCP);

    }

    if (totSymbolsInBWT < totSymbols) {
      std::cerr << "readEGSA: ERROR totSymbolsInBWT= " << totSymbolsInBWT << " < totSymbols= " << totSymbols << std::endl;
      exit (EXIT_FAILURE);
    }

    #if (verboseEncode == 1)
      std::cerr << "TableOcc after readEGSA \n";
      for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
        for (dataTypedimAlpha h = 0 ; h < sizeAlpha; h++)
          std::cerr << tableOcc[j][h] << "\t";
        std::cerr << "\n";
      }
    #endif


    delete [] filenameOutBWT;
    delete [] filenameBWT;
    delete [] filenameOutSA;
    delete [] filenameSA;
    delete [] filenameOutLCP;
    delete [] filenameLCP;



    fclose(f_ESA);

    return totSymbols;
}
#endif





int BCRexternalBWT::createFilePartialBWT() {
  //Creates one file (partial BWT) for each letter in the alphabet. From 1 to sizeAlpha-1
  static FILE *OutFileBWT;                  // output and input file BWT;
  char *filenameOut = new char[12];
  char *filename1 = new char[8];
  const char *ext = ".aux";
  for (dataTypedimAlpha i = 0; i < sizeAlpha; i++) {
    sprintf (filename1, "bwt_%d", i);
    sprintf (filenameOut,"%s%s",filename1,ext);

    OutFileBWT = fopen(filenameOut, "wb");
    if (OutFileBWT==NULL) {
      std::cerr << "BWT file " << (unsigned int)i <<" : Error opening " << std::endl;
      exit (EXIT_FAILURE);
    }
    fclose(OutFileBWT);
  }
  delete [] filename1;
  delete [] filenameOut;



  //std::cerr << "\nPartial File name for input: " << fileOut <<" \n\n";
  return 1;
}


FILE * BCRexternalBWT::openWriteFilePartialBWT_0( ) {
  static FILE *OutFileBWT;                  // output and input file BWT;
  char *filenameOut = new char[12];
  char *filename = new char[8];

  const char *ext = ".aux";
  sprintf (filename, "bwt_%d",0);
  sprintf (filenameOut,"%s%s",filename,ext);

  OutFileBWT = fopen(filenameOut, "wb");
  if (OutFileBWT==NULL) {
    std::cerr << "BWT file $: Error opening: " << filenameOut << std::endl;
    exit (EXIT_FAILURE);
  }
  delete [] filenameOut;
  delete [] filename;

  return OutFileBWT;
}


dataTypeNChar BCRexternalBWT::writeFilePartial(uchar * newSymb, FILE * OutFile) {
  dataTypeNChar num = fwrite (newSymb, sizeof(uchar), nExamedTexts , OutFile);
  return num;
}



FILE * BCRexternalBWT::openFilePartialIn(dataTypedimAlpha currentPile) {
  static FILE *inFile;
  char *filenameIn = new char[12];
  char *filename = new char[8];
  const char *ext = ".aux";
  sprintf (filename, "bwt_%d", currentPile);
  sprintf (filenameIn,"%s%s",filename,ext);
  inFile = fopen(filenameIn, "rb");
  if (inFile==NULL) {
    std::cerr << "openFilePartialIn: file currentPile=" << (unsigned int)currentPile << ": Error opening: " << filenameIn << std::endl;
    exit (EXIT_FAILURE);
  }
  delete [] filenameIn;
  delete [] filename;
  return inFile;
}


int BCRexternalBWT::closeFilePartial(FILE * pFile) {
  fclose(pFile);
  return 1;
}


FILE * BCRexternalBWT::openFilePartialOut(dataTypedimAlpha currentPile) {
  static FILE *outFile;
  char *filenameOut = new char[12];
  char *filename = new char[8];
  const char *ext = ".aux";
  sprintf (filename, "new_bwt_%d", currentPile);
  sprintf (filenameOut,"%s%s",filename,ext);
  outFile = fopen(filenameOut, "ab");
  if (outFile==NULL) {
    std::cerr << "openFilePartialOut: file currentPile= " << (unsigned int)currentPile << ": Error opening: " << filenameOut << std::endl;
    exit (EXIT_FAILURE);
  }
  delete [] filenameOut;
  delete [] filename;
  return outFile;
}


int BCRexternalBWT::renameFilePartial(dataTypedimAlpha currentPile) {
  static FILE *OutFileBWT;
  char *filenameIn = new char[12];
  char *filenameOut = new char[20];
  char *filename = new char[8];
  const char *ext = ".aux";
  sprintf (filename, "bwt_%d", currentPile);
  sprintf (filenameIn,"%s%s",filename,ext);
  sprintf (filenameOut,"new_%s%s",filename,ext);
  //std::cerr << "Filenames:" << filenameIn << "\t" <<filenameOut << std::endl;
  OutFileBWT = fopen(filenameOut, "rb");
  if (OutFileBWT!=NULL) { //If it exists
    fclose(OutFileBWT);
    if (remove(filenameIn)!=0)
      std::cerr << filenameIn <<": Error deleting file" << std::endl;
    else
      if(rename(filenameOut,filenameIn))
        std::cerr << filenameOut <<": Error renaming " << std::endl;
  }
  delete [] filenameIn;
  delete [] filename;
  delete [] filenameOut;
  return 1;
}


 dataTypeNChar BCRexternalBWT::readOnFilePartial(uchar *buffer, dataTypeNChar toRead, FILE * InFileBWT) {
  dataTypeNChar numchar;
  if (InFileBWT==NULL) {
      std::cerr << "readOnFile: file" << std::endl;
      exit (EXIT_FAILURE);
    }
  numchar = fread(buffer,sizeof(uchar),toRead,InFileBWT);

  return numchar;
}

 dataTypeNChar BCRexternalBWT::writeOnFilePartial(uchar *buffer, dataTypeNChar numchar, FILE * OutFileBWT) {
  dataTypeNChar numcharWrite;
  numcharWrite = fwrite (buffer, sizeof(uchar), numchar , OutFileBWT);
  return numcharWrite;
}

dataTypeNChar BCRexternalBWT::writeSymbolOnFilePartial(uchar symbol, dataTypeNChar numchar, FILE * OutFileBWT) {
  dataTypeNChar numcharWrite;
  numcharWrite = fwrite (&symbol, sizeof(uchar), numchar , OutFileBWT);
  return numcharWrite;
}


BCRexternalBWT::~BCRexternalBWT() {
   }
