#include "TSFCommandLine.h"


#include <cstring>
#include "TSFError.h"
#include "TSFOut.h"

#if HAVE_MPI
#include "mpi.h"
#endif

using namespace TSF;


// define static vars
TSFArray<string> TSFCommandLine::tokens_;
TSFHashtable<string, int> TSFCommandLine::tokenMap_;
bool TSFCommandLine::frozen_ = false;


// set argc and argv. This should only be done once, so the class
// is frozen when done. 

void TSFCommandLine::init(int argc, void** argv)
{
	if (frozen_) TSFError::raise("TSFCommandLine::set should never be called more than once");
	
#if HAVE_MPI
	MPI_Bcast(&argc, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int rank  = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int count=0;
	for (int i=0; i<(argc-1); i++)
		{
			if (argv[i+1]==0) continue; 

			if (rank==0)
				{
					tokens_.append(string((char*)argv[i+1]));
				}
			else
				{
					tokens_.append(" ");
				}

			int len;
			MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&(tokens_[count][0]), len+1, MPI_CHAR, 0, MPI_COMM_WORLD);

			tokenMap_.put(tokens_[count], count);
			count++;
		}
#endif
	frozen_ = true;
}


void TSFCommandLine::checkInitialization()
{
	if (!frozen_) TSFError::raise("TSFCommandLine::init(argc, argv) was not called?");
}

// check for presence of string on command line

bool TSFCommandLine::find(const string& str)
{
	checkInitialization();
	return tokenMap_.containsKey(str);
}

bool TSFCommandLine::find(const string& str, int& position) 
{
	checkInitialization();
	if (!tokenMap_.containsKey(str)) return false;
	position = tokenMap_.get(str);
	return true;
	
}


// return value via reference. Returns a bool indicating success or failure.

bool TSFCommandLine::findDouble(const string& str, double& val) 
{
	int position;

	checkInitialization();

	bool rtn = tokenMap_.containsKey(str);
	if (rtn)
		{
			position = tokenMap_.get(str);
			val = atof(tokens_[position+1].c_str());
		}
	return rtn;
}

bool TSFCommandLine::findInt(const string& str, int& val) 
{
	int position;

	checkInitialization();

	bool rtn = tokenMap_.containsKey(str);
	if (rtn)
		{
			position = tokenMap_.get(str);
			val = atoi(tokens_[position+1].c_str());
		}
	return rtn;
}

bool TSFCommandLine::findString(const string& str, string& val) 
{
	int position;

	checkInitialization();


	bool rtn = tokenMap_.containsKey(str);
	if (rtn)
		{
			position = tokenMap_.get(str);
			val = tokens_[position+1];
		}
	return rtn;
}

// grab a list of arguments

bool TSFCommandLine::findDouble(const string& str, TSFArray<double>& val, int count) 
{
	int position;

	checkInitialization();
	
	bool rtn = tokenMap_.get(str);
	if (rtn)
		{
			val.resize(count);
			for (int i=0; i<count; i++) val[i] = atof(tokens_[position+1+i].c_str());
		}
	return rtn;
}

bool TSFCommandLine::findInt(const string& str, TSFArray<int>& val, int count)
{
	int position;
	
	checkInitialization();
		
	bool rtn = tokenMap_.get(str);
	if (rtn)
		{
			val.resize(count);
			for (int i=0; i<count; i++) val[i] = atoi(tokens_[position+1+i].c_str());
		}
	return rtn;

}

bool TSFCommandLine::findString(const string& str, TSFArray<string>& val, int count)
{
	int position;

	checkInitialization();
		
	bool rtn = tokenMap_.get(str);
	if (rtn)
		{
			val.resize(count);
			for (int i=0; i<count; i++) val[i] = tokens_[position+1+i];
		}
	return rtn;
	
}

void TSFCommandLine::print() 
{
	TSFOut::println(toString(tokens_));
}
			
