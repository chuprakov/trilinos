#ifndef TSFVERSION_H
#define TSFVERSION_H

#include "TSFConfig.h"
#include "TSFOut.h"
//bvbw added:
#include <iostream>

namespace TSF
{
	static const char versionID_[] = "Hilbert v. 1.0.3 (03/06/02)";

	inline void showVersion() {std::cerr << versionID_ << std::endl;}
}

/* Revision history */

/* 11/27/01 version 1.0.1 - base version, works on linux and solaris */

/* 01/03/02 version 1.0.2 - 
 * - updated to work with Epetra. 
 * - Petra distributed vectors now work correctly.
 * - fixed dtor-time bug in the use of Ifpack preconditioners.
 */

/* 03/06/02 version 1.0.3 - 
 * - port to gnu 3.0.4
 */

#endif
