// //////////////////////////////////////////////////////////////
// StopWatchPack_stopwatch.hpp
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

#ifndef STPWATCH_H
#define STPWATCH_H

namespace StopWatchPack {

/** @name namespace StopWatchPack
  *
  * @memo Package for CPU timing.
  */
//@{

double seconds(void);

///
/** Simple stopwatch object.
  */
class stopwatch {
public:

	/// Initializes of not running.
	stopwatch() : running_(false), last_time_(0.0), total_(0.0)
	{}

	/// Returns true if <tt>this</tt> is currently timming.
	bool is_running() const {
		return running_;
	}

	/// Starts timing if it has already not been started.
	void start() {
		if (!running_) {
			last_time_ = seconds();
			running_ = true;
		}
	}

	/// Stops timing and returns the time (sec.) since start() was called
	double stop()  {
		if (running_) {
			total_ += seconds() - last_time_; 
			running_ = false;
		}
		//std::cout << "total = " << total_ << std::endl;
		return total_; 
	}

	/// Stops and resets the clock if it is running.
	void reset() {
		running_ = false;
		last_time_ = 0.0;
		total_ = 0.0;
	}

	/// Reads the elapsed time (sec.) and leaves the clock running.
	double read() {
		if (running_) {
			double curr_time = seconds();
			total_ += curr_time - last_time_;
			last_time_ = curr_time;
		}
		return total_;
	}
	
private:
	bool running_;
	double last_time_;
	double total_;
};

//	end namespace StopWatchPack 
//@}

}  // end namespace StopWatchPack 

#endif // STPWATCH_H
