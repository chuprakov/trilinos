// //////////////////////////////////////////////////////////
// ExampleNLPFirstOrderInfoRun.cpp
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

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	
// disable VC 5.0 warnings about truncated identifier names (templates).
#pragma warning(disable : 4503)	

#include <assert.h>

#include <fstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "ExampleNLPFirstOrderInfoRun.h"
#include "ExampleNLPFirstOrderInfo.h"
#include "ExampleBasisSystem.h"
#include "NLPInterfacePack/test/test_nlp_first_order_info.h"
#include "NLPInterfacePack/test/test_basis_system.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "OptionsFromStream.h"
#include "StringToBool.h"

bool NLPInterfacePack::ExampleNLPFirstOrderInfoRun(
	const VectorSpace&   vec_space
	,value_type          xo
	,bool                has_bounds
	,bool                dep_bounded
	,std::ostream*       out_in
	,std::ostream*       eout_in
	)
{
	using std::endl;
	using std::setw;
	namespace rcp = ReferenceCountingPack;
	using rcp::ref_count_ptr;
	namespace ofsp = OptionsFromStreamPack;
	using ofsp::OptionsFromStream;
	using ofsp::StringToBool;

	bool prog_return = true;

	int err = 0;
	
	int w = 15;
	int prec = 8;

	if(out_in)
		*out_in
			<< std::setprecision(prec)
			<< std::scientific
			<< "*****************************************************\n"
			<< "*** Running Tests on Example First Order Info NLP ***\n"
			<< "*****************************************************\n";

	// Read in the options
	std::ifstream      options_in_file("ExampleNLPFirstOrderInfoRun.opt");	
	OptionsFromStream  options(options_in_file);

	bool suppress_output = false;
	const std::string                    optgrp_name = "ExampleNLPFirstOrderInfoRun";
	OptionsFromStream::options_group_t   optgrp = options.options_group( optgrp_name );
	if( OptionsFromStream::options_group_exists( optgrp ) ) {
		const std::string val = optgrp.option_value("suppress_output");
		if(OptionsFromStream::options_group_t::option_exists(val))
			suppress_output = StringToBool( "suppress_output", val.c_str() );
	}

	std::ostream
		*out  = NULL,
		*eout = NULL;
	if(suppress_output) {
		if(out_in)
			*out_in << "\nOption ExampleNLPFirstOrdeInfoRun::suppress_out = true was set, suppressing all future output!\n";
		out  = NULL;
		eout = NULL;
	}
	else {
		out  = out_in;
		eout = eout_in;
	}

	// Create the nlp
	ExampleNLPFirstOrderInfo
		nlp(VectorSpace::space_ptr_t(&vec_space,false),xo,has_bounds,dep_bounded);

	// Test the NLPFirstOrderInfo interface!
	bool
		result = NLPInterfacePack::test_nlp_first_order_info(&nlp,&options,out);
	if(!result)
		prog_return = false;

	// Create the basis system
	ExampleBasisSystem
		basis_sys(rcp::rcp(&vec_space,false));

	// Test the basis system
	result = NLPInterfacePack::test_basis_system(&nlp,&basis_sys,&options,out);
	if(!result)
		prog_return = false;

	if(prog_return == true) {
		if(eout_in && eout_in != out_in)
			*eout_in   << "Congradulations! The VectorSpace, NLP and BasisSystem objects check out!\n";
		if(out_in)
			*out_in    << "\nCongradulations! The VectorSpace, NLP and BasisSystem objects check out!\n";
	}
	else {
		if(eout_in && eout_in!= out_in)
			*eout_in   << "Oh No!  Something did not check out!\n";
		if(out_in)
			*out_in    << "\nOh No!  Something did not check out!\n";
	}

	return prog_return;
}
