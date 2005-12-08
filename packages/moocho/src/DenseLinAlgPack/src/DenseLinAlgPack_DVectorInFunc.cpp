// /////////////////////////////////////////////////////////////////////////////////
// DVectorInFunc.cpp
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

#include "DenseLinAlgPack_DVectorInFunc.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"

namespace {	// Local implementation
std::istream& input_vs(std::istream& is, DenseLinAlgPack::DVectorSlice* vs, const char func[]);
}

std::istream& DenseLinAlgPack::input(std::istream& is, DVector* v, LinAlgPackIO::fmtflags extra_flags) {
	if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) ) {
		size_type n;
		is >> n;
		if(is.fail())
			throw LinAlgPackIO::InputException("DenseLinAlgPack::input() {DVector}:  Input operation of vector dimension failed.  Check that the constant n is a valid integer.");
		if(is.bad())
			throw std::ios_base::failure("DenseLinAlgPack::input() {DVector}: Input operation failed because the stream became currupted.");
		v->resize(n);
	}
	return input_vs(is,&(*v)(),"DenseLinAlgPack::input() {DVector}");
}

std::istream& DenseLinAlgPack::input(std::istream& is, DVectorSlice* vs, LinAlgPackIO::fmtflags extra_flags) {
	if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) ) {
		size_type n;
		is >> n;
		if(is.fail())
			throw LinAlgPackIO::InputException("DenseLinAlgPack::input() {DVectorSlice}:  Input operation of vector dimension failed.  Check that the constant n is a valid integer.");
		if(is.bad())
			throw std::ios_base::failure("DenseLinAlgPack::input() {DVectorSlice}: Input operation failed because the stream became currupted.");
		DenseLinAlgPack::Vp_V_assert_sizes( vs->dim(), n );
	}
	return input_vs(is,vs,"DenseLinAlgPack::input() {DVectorSlice}");
}


// ///////////////////////////////////
// Local implementation

namespace {

// Read in a specified number of elements into a DVectorSlice object.
// The dim of vs is not checked.  If an element input operation fails or the end of the file
// is reached before all of the elements are read in then a LinAlgPackIO::InputException is thrown.
// If the stream becomes currupted durring the input then a std::ios_base::failure exception
// is thrown.  The state of the input steam remains the same on return accept for the char's
// that have been extracted.
std::istream& input_vs(std::istream& is, DenseLinAlgPack::DVectorSlice* vs, const char func[]) {
	using std::ios_base;
	using DenseLinAlgPack::DVectorSlice;
	if(!vs->dim()) return is;	// If there are no elements to read in just return
	ios_base::iostate old_state = is.exceptions();		// save the old state
	is.exceptions(ios_base::badbit | ios_base::failbit);
	try {
		// Read in the elements
		for(DVectorSlice::iterator itr = vs->begin(); itr != vs->end(); ++itr)
			is >> *itr;
	}
	catch(std::ios_base::failure& excpt) {
		is.exceptions(old_state);
		if(is.bad()) throw;	// The stream was bad so rethrow the exception
		if(is.fail()) {
			std::ostringstream os;
			os << func << ":  An vector element input failed.  Check that the vector element is a valid C number.  "
			   << excpt.what();
			throw DenseLinAlgPack::LinAlgPackIO::InputException(os.str());			
		}
		if(is.eof()) {
			std::ostringstream os;
			os << func << ":  DVector input failed.  The end of the file was found before all of the elements where read in.  "
			   << excpt.what();;
			throw DenseLinAlgPack::LinAlgPackIO::InputException(os.str());			
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
