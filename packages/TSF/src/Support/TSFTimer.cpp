#include "TSFTimer.h"
#include "TSFMPI.h"

#include "TSFOut.h"

using namespace TSF;

TSFArray<TSFSmartPtr<TSFTimer> > TSFTimer::timers_;

TSFSmartPtr<TSFTimer> TSFTimer::getNewTimer(const string& name)
{
	TSFSmartPtr<TSFTimer> rtn = new TSFTimer(name);
	timers_.append(rtn);
	return rtn;
}


void TSFTimer::summarize()
{
	TSFArray<string> names(timers_.size());
	TSFArray<double> timings(timers_.size());

	for (int i=0; i<timers_.size(); i++)
		{
			names[i] = timers_[i]->name();
			timings[i] = timers_[i]->getTime();
		}

	int np=1;
	int rank=0;
#ifdef HAVE_MPI
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	if (np==1)
		{
			for (int i=0; i<names.size(); i++)
				{
					TSFOut::printf("%-40s: %g\n", names[i].c_str(), timings[i]);
				}
		}
	else
		{
			TSFArray<double> minTime(timers_.size());
			TSFArray<double> maxTime(timers_.size());
			TSFArray<double> avgTime(timers_.size());
			gatherTimings(timings, minTime, avgTime, maxTime);
			if (rank==0)
				{
					for (int i=0; i<names.size(); i++)
						{
							TSFOut::printf("%-30s: %-12g %-12g %-12g\n", names[i].c_str(), 
														 minTime[i], avgTime[i], maxTime[i]);
						}
				}
		}

}

void TSFTimer::gatherTimings(const TSFArray<double>& timings,
														 TSFArray<double>& minTime,
														 TSFArray<double>& avgTime,
														 TSFArray<double>& maxTime)
{
#ifdef HAVE_MPI
	int np;
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	void* tPtr = (void*) &(timings[0]);
	void* minPtr = (void*) &(minTime[0]);
	void* avgPtr = (void*) &(avgTime[0]);
	void* maxPtr = (void*) &(maxTime[0]);

	int count = (int) timings.size();

	MPI_Reduce(tPtr, minPtr, count, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(tPtr, avgPtr, count, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(tPtr, maxPtr, count, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	for (int i=0; i<avgTime.size(); i++)
		{
			avgTime[i] = avgTime[i]/((double) np);
		}
#endif
}

