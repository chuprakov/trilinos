#include "TSFOut.h"
#include "TSFDefaultWriter.h"

#include "TSFError.h"
#include "TSFMPI.h"

using namespace TSF;



/* initialize the static raise handler object to the default handler. This
 * can be changed later with a call to setRaiseHandler() */

TSFSmartPtr<TSFWriterBase> TSFOut::writer_ = new TSFDefaultWriter();


void TSFOut::print(const std::string& msg)
{
	writer_->print(msg);
}

void TSFOut::println(const std::string& msg)
{
	writer_->println(msg);
}

void TSFOut::rootPrintln(const std::string& msg)
{
	int rank=0;
#if HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	if (rank==0)
		{
			writer_->println(msg);
		}
}

void TSFOut::rootPrintf(const char* format, ...)
{
	int rank=0;
#if HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	if (rank==0)
		{
			va_list args;
			va_start(args, format);
			TSFOut::vprintf(format, args);
		}
}

void TSFOut::printf(const char* format, ...)
{
	va_list args;
	va_start(args, format);

	int bufSize = 100;

	while (bufSize > 0)
		{
			char* str = new char[bufSize];
			int rtn = vsnprintf(str, bufSize, format, args);
			if (rtn > 0)
				{
					writer_->print(str);
					va_end(args);
					delete [] str;
					return;
				} 
			else if (rtn==0)
				{
					va_end(args);
					delete [] str;
					return;
				}
			else
				{
					bufSize *= 2;
				}
		}
	
	TSFError::raise("buffer overflow in TSFOut::printf()");
}

void TSFOut::vprintf(const char* format, va_list args)
{
	int bufSize = 100;

	while (bufSize > 0)
		{
			char* str = new char[bufSize];
			int rtn = vsnprintf(str, bufSize, format, args);
			if (rtn > 0)
				{
					writer_->print(str);
					va_end(args);
					delete [] str;
					return;
				} 
			else if (rtn==0)
				{
					va_end(args);
					delete [] str;
					return;
				}
			else
				{
					bufSize *= 2;
				}
		}
	
	TSFError::raise("buffer overflow in TSFOut::printf()");
}


