#ifndef TSFTIMER_H
#define TSFTIMER_H

#include "TSFConfig.h"
#include <string>
#include <stdexcept>
#include <time.h>
#include "TSFWriterBase.h"
#include "TSFSmartPtr.h"
#include "TSFArray.h"


namespace TSF
{
	
	using std::string;
	

	/** \ingroup IO
	 *
	 */
	class TSFTimer
		{
		public:
			/** */
			TSFTimer(const string& name) 
				: total_(0.0), t0_(0.0), isRunning_(false), name_(name) {;}
			/** empty ctro for use in containers */
			TSFTimer() 
				: total_(0.0), t0_(0.0), isRunning_(false), name_() {;}

			/** */
			inline double getTime() const {return total_/((double) CLOCKS_PER_SEC);}

			/** */
			inline void reset() {total_=0.0; t0_=0.0; isRunning_=false;}

			/** */
			inline void start() 
				{
					t0_ = clock(); 
					isRunning_=true;
				}

			/** */
			inline void stop() 
				{
					total_ += clock() - t0_;
					isRunning_=false;
				}

			/** */
			inline bool isRunning() const {return isRunning_;}

			inline const string& name() const {return name_;}

			/** Print summary statistics for a group of timers. Timings are gathered
			 * from all processors */
			static void summarize();

			/** */
			static TSFSmartPtr<TSFTimer> getNewTimer(const string& name);
		private:

			static void gatherTimings(const TSFArray<double>& timings,
																TSFArray<double>& minTime,
																TSFArray<double>& avgTime,
																TSFArray<double>& maxTime);
			double total_;
			double t0_;
			bool isRunning_;
			string name_;

			static TSFArray<TSFSmartPtr<TSFTimer> > timers_;
		};

}
#endif
