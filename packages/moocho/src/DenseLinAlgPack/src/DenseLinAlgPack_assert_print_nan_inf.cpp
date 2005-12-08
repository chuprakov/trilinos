// ///////////////////////////////////////////////////////////
// assert_print_nan_inf.cpp
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
#include <sstream>
#include <iomanip>

#include "DenseLinAlgPack_assert_print_nan_inf.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "check_nan_inf.h"

bool DenseLinAlgPack::assert_print_nan_inf( const value_type& val, char name[]
	, bool throw_excpt, std::ostream* out )
{
	
	if( RTOp_is_nan_inf(val) ) {
		std::ostringstream omsg;
		omsg
			<< "The scalar \"" << name
			<< "\" = " << val << " is not a valid bounded number";
		if(out)
			*out << omsg.str() << std::endl;
		if( throw_excpt ) {
			if(out)
				out->flush();	
			throw NaNInfException( "assert_print_nan_inf(...) : Error, "
				+ omsg.str() );
		}
		return false;
	}
	return true;
}

bool DenseLinAlgPack::assert_print_nan_inf( const DVectorSlice& v, char name[]
	, bool throw_excpt, std::ostream* out )
{
	
	bool has_nan_or_inf = false;
	bool printed_header = false;

	for( DVectorSlice::const_iterator v_itr = v.begin(); v_itr != v.end(); ++v_itr ) {
		if( RTOp_is_nan_inf(*v_itr) ) {
			if(out) {
				if(!printed_header) {
					*out
						<< "The vector \"" << name
						<< "\" has the following NaN or Inf entries\n";
					printed_header = true;
				}
				*out
					<< name << "(" << v_itr - v.begin() + 1 << ") = "
					<< *v_itr << std::endl;
			}
			has_nan_or_inf = true;
		}
	}
	if( has_nan_or_inf && throw_excpt ) {
		if(out)
			out->flush();	
		std::ostringstream omsg;
		omsg
			<< "assert_print_nan_inf(...) : Error, the vector named "
			<< name << " has at least one element which is NaN or Inf";
		throw NaNInfException( omsg.str() );
	}

	return !has_nan_or_inf;
}

bool DenseLinAlgPack::assert_print_nan_inf( const DMatrixSlice& m, char name[]
	, bool throw_excpt, std::ostream* out )
{
	
	bool has_nan_or_inf = false;
	bool printed_header = false;

	for( size_type j = 1; j <= m.cols(); ++j ) {
		const DVectorSlice& v = m.col(j);
		for( DVectorSlice::const_iterator v_itr = v.begin(); v_itr != v.end(); ++v_itr ) {
			if( RTOp_is_nan_inf(*v_itr) ) {
				if(out) {
					if(!printed_header) {
						*out
							<< "The matrix \"" << name
							<< "\" has the following NaN or Inf entries\n";
						printed_header = true;
					}
					*out
						<< name << "(" << v_itr - v.begin() + 1 << "," << j << ") = "
						<< *v_itr << std::endl;
				}
				has_nan_or_inf = true;
			}
		}
	}
	
	if( has_nan_or_inf && throw_excpt ) {
		if(out)
			out->flush();	
		std::ostringstream omsg;
		omsg
			<< "assert_print_nan_inf(...) : Error, the matrix named "
			<< name << " has at least one element which is NaN or Inf";
		throw NaNInfException( omsg.str() );
	}

	return has_nan_or_inf;
}
