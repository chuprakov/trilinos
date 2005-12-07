// ///////////////////////////////////////////////////////////
// TestIterationPack.cpp
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

#include <ostream>

#include "IterationPack_TestIterationPack.hpp"
#include "TestingHelperPack_update_success.hpp"

bool IterationPack::TestingPack::TestIterationPack(std::ostream* out)
{
	using TestingHelperPack::update_success;
	using namespace IterationPack::TestingPack;

	bool success = true;
	update_success( TestIterQuantityAccessContiguous(out), &success );
	update_success( TestAlgorithmState(out), &success );
	update_success( TestAlgorithm(out), &success );

	if(out) {
		if(success) {
			*out << "\n*** Congradulations, IterationPack seems to"
					" check out.\n";
		}
		else {
			*out << "\n*** Oops, at least one of the above tests for IterationPack "
					" did not return the expected results\n";
		}
	} 

	return success;
}
