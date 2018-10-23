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
 
#include "Timer.hh"



Timer::Timer(void)
{
    if (gettimeofday(&lastTime_, NULL) != 0)
        exit(-1);
    if (getrusage(RUSAGE_SELF, &lastUsage_) != 0)
        exit(-1);
} // ~Timer::Timer( void )


std::ostream& Timer::print(std::ostream& os)
{
    if (gettimeofday(&thisTime_, NULL) != 0)
        exit(-1);
    if (getrusage(RUSAGE_SELF, &thisUsage_) != 0)
        exit(-1);

    static double elapsedActual, elapsedUser, elapsedSystem;

    elapsedUser = thisUsage_.ru_utime.tv_sec - lastUsage_.ru_utime.tv_sec
            + (thisUsage_.ru_utime.tv_usec
                    - (double) lastUsage_.ru_utime.tv_usec) / 1000000;

    elapsedSystem = thisUsage_.ru_stime.tv_sec - lastUsage_.ru_stime.tv_sec
            + (thisUsage_.ru_stime.tv_usec
                    - (double) lastUsage_.ru_stime.tv_usec) / 1000000;

    elapsedActual = thisTime_.tv_sec - lastTime_.tv_sec + (thisTime_.tv_usec
            - (double) lastTime_.tv_usec) / 1000000;

    os << "User: " << elapsedUser << "s System: " << elapsedSystem
            << "s Actual: " << elapsedActual << "s Efficiency: "
            << ((elapsedActual == 0) ? 0 : ((elapsedUser + elapsedSystem)
                    * 100.0 / elapsedActual)) << '\045'; // only way to print % liked by both Intel and gcc!
    lastTime_ = thisTime_;
    lastUsage_ = thisUsage_;
    return os;
} // ~std::ostream& Timer::print( std::ostream& os )

std::ostream& operator<<(std::ostream& os, Timer& timer)
{
    return timer.print(os);
} // ~std::ostream& operator<<( std::ostream& os, Timer& timer )

const char* Timer::timeNow(void) const
{
    time_t tt(time( NULL));
    return ctime(&tt);
} // ~const char* Timer::timeNow( void ) const
