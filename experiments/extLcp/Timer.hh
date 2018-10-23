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
 
#ifndef INCLUDED_TIMER
#define INCLUDED_TIMER

#include <fstream>
#include <ctime>
#include <sys/resource.h>
#include <sys/time.h>

#include <cstdlib>

// Class Name : Timer
// Description: Maintains info on actual and processing time
class Timer
{

  public:
  Timer( void );
  std::ostream& print( std::ostream& os );
  // timeNow: returns current date and time as an ASCII string
  const char* timeNow( void ) const;

  private:
  
  rusage thisUsage_;
  rusage lastUsage_;
  timeval thisTime_;
  timeval lastTime_;
  //  timeb thisTime_;
  //  timeb lastTime_;

}; // Timer

std::ostream& operator<<( std::ostream& os, Timer& timer );

#endif
// end of Timer.hh
