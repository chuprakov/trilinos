// ////////////////////////////////////////////////////////////////////
// IterationPack_TestIterationPack.hpp
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

#ifndef TEST_GENERAL_ITERATION_PACK_H
#define TEST_GENERAL_ITERATION_PACK_H

#include <iosfwd>

namespace IterationPack {
	namespace TestingPack {
		/// Test all of IterationPack
		bool TestIterationPack( std::ostream* out );
		///
		bool TestIterQuantityAccessContiguous( std::ostream* out );
		///
		bool TestAlgorithmState( std::ostream* out );
		///
		bool TestAlgorithm( std::ostream* out );
	}
}

#endif // TEST_GENERAL_ITERATION_PACK_H
