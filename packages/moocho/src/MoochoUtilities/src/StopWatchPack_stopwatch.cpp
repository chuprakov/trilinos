// //////////////////////////////////////////////////////////////////////////////////////
// StopWatchPack_stopwatch.cpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

//#include <iostream>
//#include <iosfwd>

#include "StopWatchPack_stopwatch.hpp"

#ifndef _INTEL_CXX

// Implementation using C standard library.
// In MS VC++ 6.0 the precision is only about 0.05 sec.

#include <time.h>

double StopWatchPack::seconds(void)
{
    static const double secs_per_tick = ((double)1.0) / CLOCKS_PER_SEC;
	const clock_t ticks = clock();
	const double sec = ( (double) ticks ) * secs_per_tick;
	//std::cout << "ticks = " << ticks << ", sec = " << sec << std::endl;
    return sec;
}

#else	// _INTEL_CXX implementation

// Windows implementation.

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <assert.h>

namespace {

bool seconds_initialized = false;
LARGE_INTEGER start_count, count_freq;	// counts per sec.

inline void seconds_initialize() {
	if( seconds_initialized ) return;
	// Figure out how often the performance counter increments
	::QueryPerformanceFrequency( &count_freq );
	// Set this thread's priority as high as reasonably possible to prevent
    // timeslice interruptions
    ::SetThreadPriority( ::GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL );
	// Get the first count.
	assert( QueryPerformanceCounter( &start_count ) );
	seconds_initialized = true;
}

}	// end namespace

double StopWatchPack::seconds(void)
{
	seconds_initialize();
	LARGE_INTEGER count;
	QueryPerformanceCounter( &count );
	// "QuadPart" is a 64 bit integer (__int64).  VC++ supports them!
	const double
		sec = (double)( count.QuadPart - start_count.QuadPart ) / count_freq.QuadPart;
	//std::cout << "ticks = " << ticks << ", sec = " << sec << std::endl;
    return sec;
}

#endif	// _INTEL_CXX
