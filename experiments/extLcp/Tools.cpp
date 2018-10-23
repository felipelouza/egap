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

#include "Tools.h"


double Tools::startTime;

void Tools::StartTimer()
{
    static const double in_seconds = 1.0/static_cast<double>(CLOCKS_PER_SEC);
    startTime = clock() * in_seconds;
}

double Tools::GetTime()
{
    static const double in_seconds = 1.0/static_cast<double>(CLOCKS_PER_SEC);
    double curTime = clock() * in_seconds;
    return curTime - startTime;
}

/*
uchar * Tools::GetRandomString(unsigned min, unsigned max, unsigned &alphabetSize)
{
    unsigned len = std::rand() % (max - min) + min;
    alphabetSize = std::rand() % 26 + 1;
    uchar* temp = new uchar[len + 2];
    for (unsigned i = 0; i < len; i++)
        temp[i] = 97 + std::rand() % alphabetSize;
    temp[len] = 0u ;temp[len+1] = '\0';
    return temp;
}
*/

unsigned Tools::FloorLog2(ulong i)
{
    ulong b = 0;
    if (i == 0)
        return 0;
    while (i)
    { 
        b++; 
        i >>= 1; 
    }
    return b - 1;
}

//Creating table to find logn in small time
unsigned * Tools::MakeTable()
{
    unsigned *table = new unsigned[512];
    for(unsigned i = 0; i < 9; i++)
    {
        if (i == 0)
           table[i] = 0;
        if (i >= 1 && i < (1 << 1 )) 
           table[i] = 1;
        if (i >= (1 << 1 ) && i < (1 << 2 ))
           table[i] = 2;
        if (i >= (1 << 2 ) && i < (1 << 3 ))
           table[i] = 3;
        if (i >= (1 << 3 ) && i < (1 << 4 ))
           table[i] = 4;
        if (i >= (1 << 4 ) && i < (1 << 5 ))
           table[i] = 5;
        if (i >= (1 << 5 ) && i < (1 << 6 ))
           table[i] = 6;
        if (i >= (1 << 6 ) && i < (1 << 7 ))
           table[i] = 6;
        if (i >= (1 << 7 ) && i < (1 << 8 ))
           table[i] = 7;
        if (i >= (1 << 8 ) && i < (1 << 9 ))
           table[i] = 8;
    }
        return table;
}

unsigned Tools::FastFloorLog2(unsigned i)
{
    unsigned *table = MakeTable(); unsigned u;
    if (i >> 24)    u = 22 + table[ i >> 24] ;
    if (i >> 16)    u = 14 + table[ i >> 16] ;
    if (i >> 8)     u = 6 + table[ i >> 8] ;
    u =  table[i] - 1;
    delete [] table;
    return u;
}

unsigned Tools::CeilLog2(ulong i)
{
    unsigned j = FloorLog2(i);
    if ((ulong)(1lu << j) != i)
        return j + 1;
        
    return j;
}

uchar * Tools::GetFileContents(char *filename, ulong maxSize)
{
    std::ifstream::pos_type posSize;
    std::ifstream file ((char *)filename, std::ios::in|std::ios::binary|std::ios::ate);
    if (file.is_open())
    {
        posSize = file.tellg();
        ulong size = posSize;
        if (maxSize != 0 && size > maxSize)
            size = maxSize;
        char *memblock = new char [size + 1];
        file.seekg (0, std::ios::beg);
        file.read (memblock, size);
        memblock[size] = '\0';
        file.close();
	return (uchar *)memblock;
    }
    else
        return 0;
}



