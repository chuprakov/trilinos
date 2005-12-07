// ///////////////////////////////////////////////////
// TestStopWatchMain.cpp
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

#include "StopWatchPack_stopwatch.hpp"

int main() {
	using std::cout;
	using std::endl;
	using StopWatchPack::stopwatch;

	stopwatch timer;

	// Finding min resolution.
	double min_resolution = 0.0; // in seconds
	int total_num_calls   = 0;
	{
		cout	<< "\n*** Measuring miminum resolution.\n";
		timer.start();
		double last_time = timer.read();
		const int max_num_samples = 20;
		int num_samples = 0;
		int num_calls = 0;
		while( num_samples < max_num_samples ) {
			double time = timer.read();
			num_calls++;
			if( time - last_time > 0.0 ) {
				cout	<< "time_diff = " << time - last_time
						<< ", num_calls = " << num_calls << endl;
				min_resolution += time - last_time;
				++total_num_calls;
				last_time = time;
				num_calls = 0;
				num_samples++;
			}
		}
		min_resolution /= total_num_calls;
	}

	std::cerr << "Minimum stopwatch resolution = " << min_resolution << " sec\n";

	// Finding increasing resolution.
	{
		cout	<< "\n*** Measuring increasing resolution.\n";
		timer.start();
		double start_time = timer.read(), last_time = start_time;
		const int max_num_samples = 20;
		int num_samples = 0;
		int num_calls = 0;
		while( num_samples < max_num_samples ) {
			double time = timer.read();
			num_calls++;
			if( time - last_time > 0.0 ) {
				cout	<< "time = " << time - start_time
						<< ", num_calls = " << num_calls << endl;
				last_time = time;
				num_calls = 0;
				num_samples++;
			}
		}
	}

	return 0;
}
