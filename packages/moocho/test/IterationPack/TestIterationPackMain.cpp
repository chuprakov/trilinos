// ///////////////////////////////////////////////////////////////////
// TestIterationPackMain.cpp
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

#include <iostream>

#include "IterationPack_TestIterationPack.hpp"
#include "TestingHelperPack_update_success.hpp"

int main() {
	using TestingHelperPack::update_success;
	using namespace IterationPack::TestingPack;

	std::ostream* out = &std::cout;

	bool success = true;
	update_success( TestIterationPack(out), &success );

	if( success )
		std::cerr << "IterationPack seems to check out!\n";
	else
		std::cerr << "Oops! At least one of the tests in IterationPack failed!\n";
		
	return success == true ? 0 : -1;
}
