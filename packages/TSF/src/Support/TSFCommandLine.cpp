#include "TSFCommandLine.h"


#include <cstring>
#include "TSFError.h"
#include "TSFOut.h"

#ifdef HAVE_MPI
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
	if (frozen_) return;
	
#ifdef HAVE_MPI
	MPI_Bcast(&argc, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int rank  = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (verbose())
		{
			TSFOut::rootPrintf("number of arguments: %d\n", argc);
		}
	
	int count=0;
	for (int i=0; i<(argc-1); i++)
		{
			int len;
			if (argv[i+1]==0) 
				{
					len = -1;
				}
			else
				{
					len = strlen((char*)argv[i+1]);
				}

			if (rank==0 && len > 0)
				{
					tokens_.append(string((char*)argv[i+1]));
				}
			else
				{
					tokens_.append(" ");
				}

			if (verbose())
				{
					TSFOut::rootPrintf("root sees argument: %s\n", argv[i+1]);
					TSFOut::rootPrintf("root sees len: %d\n", len);
				}

			MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);

			if (verbose())
				{
					TSFOut::printf("len: %d\n", len);
				}

			if (len > 0)
				{
					char* tmp = new char[len+1];
					if (rank==0)
						{
							strncpy(tmp, tokens_[count].c_str(), (size_t) len+1);
						}
					MPI_Bcast((void*) tmp, len+1, MPI_CHAR, 0, MPI_COMM_WORLD);
					if (rank != 0)
						{
							tokens_[count] = string(tmp);
						}
					if (verbose())
						{
							TSFOut::printf("bcast'd argument: %s\n", 
														 tokens_[count].c_str());
						}
					tokenMap_.put(tokens_[count], count);
					delete [] tmp;
				}
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
			
