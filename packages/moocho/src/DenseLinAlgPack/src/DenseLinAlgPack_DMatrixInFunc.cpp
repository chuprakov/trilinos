// /////////////////////////////////////////////////////////////////////////////////
// GenMatrixInFunc.cpp
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

#include <sstream>

#include "Misc/include/EatInputComment.h"
#include "../include/GenMatrixInFunc.h"
#include "../include/VectorInFunc.h"
#include "../include/GenMatrixClass.h"

namespace {	// Local inplementation
std::istream& input_gms(std::istream& is, LinAlgPack::GenMatrixSlice* gms, const char func[]);
}

std::istream& LinAlgPack::input(std::istream& is, GenMatrix* gm, LinAlgPackIO::fmtflags extra_flags) {
	if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) ) {
		size_type m, n;
		is >> m >> n;
		if(is.fail())
			throw LinAlgPackIO::InputException( "LinAlgPack::input() {GenMatrix}: "
				"Input operation of matrix dimension failed.  Check that the constant n "
				"is a valid integer." );
		if(is.bad())
			throw std::ios_base::failure( "LinAlgPack::input() {GenMatrix}: "
				"Input operation failed because the stream became currupted." );
		gm->resize(m,n);
	}
	GenMatrixSlice gms = (*gm)();
	return input_gms(is,&gms,"LinAlgPack::input() {GenMatrix}");
}

std::istream& LinAlgPack::input(std::istream& is, GenMatrixSlice* gms, LinAlgPackIO::fmtflags extra_flags) {
	if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) ) {
		size_type m, n;
		is >> m >> n;
		if(is.fail())
			throw LinAlgPackIO::InputException( "LinAlgPack::input() {GenMatrixSlice}: "
				"Input operation of matrix dimension failed.  Check that the constant n "
				" is a valid integer.");
		if(is.bad())
			throw std::ios_base::failure( "LinAlgPack::input() {GenMatrixSlice}: "
				"Input operation failed because the stream became currupted." );
		LinAlgPack::assert_gms_lhs(*gms,m,n);
	}
	return input_gms( is, gms, "LinAlgPack::input() {GenMatrixSlice}" );
}

// //////////////////////
// Local implementation

namespace {

// Read in a specified number of elements into a GenMatrixSlice object.
// The dim of gms is not checked.  If an element input operation fails or the end of the file
// is reached before all of the elements are read in then a LinAlgPackIO::InputException is thrown.
// If the stream becomes currupted durring the input then a std::ios_base::failure exception
// is thrown.  The state of the input steam remains the same on return accept for the char's
// that have been extracted.
std::istream& input_gms(std::istream& is, LinAlgPack::GenMatrixSlice* gms, const char func[]) {
	using std::ios_base;
	using LinAlgPack::size_type;
	using LinAlgPack::VectorSlice;
	if(!gms->rows()) return is;	// If we are inputting an unsized matrix then there are no elements
								// to extract so just return.
	ios_base::iostate old_state = is.exceptions();		// save the old state
	is.exceptions(ios_base::badbit | ios_base::failbit | ios_base::eofbit);
	try {
		// Read in the rows
		for(size_type i = 1; i <= gms->rows(); ++i) {
			InputStreamHelperPack::eat_comment_lines(is,'*');
			VectorSlice gms_row_i = gms->row(i);
			LinAlgPack::input( is, &gms_row_i
				, (LinAlgPack::LinAlgPackIO::fmtflags)(LinAlgPack::LinAlgPackIO::ignore_dim_bit) );
		}
	}
	catch(...) {
		is.exceptions(old_state);
		throw;
	}
	is.exceptions(old_state);
	return is;
}

}	// end namespace
