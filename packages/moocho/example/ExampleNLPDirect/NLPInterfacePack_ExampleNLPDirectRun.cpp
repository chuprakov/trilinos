// //////////////////////////////////////////////////////////
// ExampleNLPFirstOrderDirectRun.cpp
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
//

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	
// disable VC 5.0 warnings about truncated identifier names (templates).
#pragma warning(disable : 4503)	

#include <assert.h>

#include <fstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "ExampleNLPFirstOrderDirectRun.h"
#include "ExampleNLPFirstOrderDirect.h"
#include "NLPInterfacePack/test/test_nlp_first_order_direct.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "OptionsFromStream.h"

bool NLPInterfacePack::ExampleNLPFirstOrderDirectRun(
	const VectorSpace&   vec_space
	,value_type          xo
	,bool                has_bounds
	,bool                dep_bounded
	,std::ostream*       out
	)
{
	using std::endl;
	using std::setw;
	namespace rcp = ReferenceCountingPack;
	using rcp::ref_count_ptr;
	namespace ofsp = OptionsFromStreamPack;
	using ofsp::OptionsFromStream;

	bool prog_return = true;

	int err = 0;
	
	int w = 15;
	int prec = 8;

	if(out)
		*out
			<< std::setprecision(prec)
			<< std::scientific
			<< "********************************************************\n"
			<< "*** Running Tests on Example First Order Dirrect NLP ***\n"
			<< "********************************************************\n";

	// Read in the options
	std::ifstream      options_in_file("ExampleNLPFirstOrderDirectRun.opt");	
	OptionsFromStream  options(options_in_file);

	// Create the nlp
	ExampleNLPFirstOrderDirect
		nlp(VectorSpace::space_ptr_t(&vec_space,false),xo,has_bounds,dep_bounded);

	// Test the NLPFirstOrderDirect interface!
	const bool
		result = NLPInterfacePack::test_nlp_first_order_direct(&nlp,&options,out);
	if(!result)
		prog_return = false;

	if(prog_return == true) {
		std::cerr   << "Congradulations!  The vector space and NLP class seems to check out!\n";
		if(out)
			*out    << "\nCongradulations!  The vector space and NLP class seems to check out!\n";
	}
	else {
		std::cerr   << "Oh No!  Something did not checkout!\n";
		if(out)
			*out    << "\nOh No!  Something did not checkout!\n";
	}

	return prog_return;
}
