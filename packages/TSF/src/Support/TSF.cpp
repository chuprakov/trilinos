#include "TSF.h"

namespace TSF
{
	void init(int argc, void** argv)
	{
		TSFMPI::init(argc, argv);
	}

	void finalize()
	{
		TSFMPI::finalize();
	}

	void handleError(exception& e, const string& file)
	{
		TSFOut::println("exception " + string(e.what()) + " detected in " + file);
	}
}
