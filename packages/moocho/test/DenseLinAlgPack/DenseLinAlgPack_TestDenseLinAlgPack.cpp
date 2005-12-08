// ////////////////////////////////////////////////////////////////////////////
// TestDenseLinAlgPack.cpp
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

#include "DenseLinAlgPack_TestDenseLinAlgPack.hpp"
#include "TestingHelperPack_update_success.hpp"

bool DenseLinAlgPack::TestingPack::TestDenseLinAlgPack(std::ostream* out)
{
	using TestingHelperPack::update_success;

	bool success = true;
//	bool result;

	if(out)
		*out
      << "\n*******************************"
      << "\n*** Testing DenseLinAlgPack ***"
      << "\n*******************************\n";

	update_success( TestVectorClass(out), &success );
	update_success( TestVectorOp(out), &success );
	update_success( TestGenMatrixClass(out), &success );
	update_success( TestGenMatrixOp(out), &success );

	if(out) {
		if(success)
			(*out)
				<< "\n*** Congradulations, DenseLinAlgPack seems to check out. ***\n";
		else
			(*out)
				<< "\n*** Oops, all of the tests for DenseLinAlgPack "
					"where not successful. ***\n";
	}

	return success;
}
