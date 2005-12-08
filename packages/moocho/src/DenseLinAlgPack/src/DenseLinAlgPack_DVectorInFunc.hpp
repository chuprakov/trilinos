// /////////////////////////////////////////////////////////////////////////////////
// DenseLinAlgPack_DVectorInFunc.hpp
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

#ifndef VECTOR_IN_FUNC_H
#define VECTOR_IN_FUNC_H

#include "DenseLinAlgPack_IOBasic.hpp"

namespace DenseLinAlgPack {

/* * @name DVector/DVectorSlice input stream functions.
  *
  * These are functions that are used to read a DVector
  * or DVectorSlice object from a formated input stream.
  *
  * The input format is diferent depending on the on whether the
  * bit #LinAlgPackIO::ignore_dim_bit# is set.  If #exta_flags & LinAlgPackIO::ignore_dim_bit != 0#
  * then the input format is:
  *
  * Case 1\\
  *	#n#\\
  * #v(1) v(2) v(3) ... v(n)#\\
  *
  * If #exta_flags & LinAlgPackIO::ignore_dim_bit == 0# then the input format is:
  *
  * Case 2\\
  * #v(1) v(2) v(3) ... v(v.size())#\\
  *
  * The numbers of the input must be seperated by white space and be valid
  * C numeric constants.  For example, the input format for the vector {1.1, 2.2, 3.3}
  * for case 1 is:
  *
  * #3#\\
  * #1.1	2.2		3.3#\\
  *
  * And for case 2 is:
  *
  * #1.1	2.2		3.3#\\
  *
  * It is permisible for the dimension #n# in case 1 to be 0.  In this case there will be
  * no elements.  So to input an empty vector you would use:
  *
  * #0#\\
  *
  * If any of the input operations fails then a LinAlgPackIO::InputException exception
  * is thrown.  In other words if #is.fail()# or #is.eof()# is true
  * before all of the elements have been read in then the exception is thrown.
  * Also if the stream becomes corrupted (#is.bad() == true#) then a #std::ios_base::failure#
  * exception is thrown. 
  */
// @{

///
/* * DVector input stream function.
  *
  * Inputs a DVector object from an input stream.  If #exta_flags & LinAlgPackIO::ignore_dim_bit != 0#
  * then #v# is resized to #n# given in the file.  If #exta_flags & LinAlgPackIO::ignore_dim_bit == 0#
  * then the number of elements read in depends on the current size of #v#.
  *
  */
std::istream& input(std::istream& is, DVector* v, LinAlgPackIO::fmtflags extra_flags);

///
/* * DVectorSlice input stream function.
  *
  * Inputs a DVectorSlice object from an input stream.  If #exta_flags & LinAlgPackIO::ignore_dim_bit != 0#
  * then the size (!= 0) of #vs# is compared to the #n# given in the file and if they are not equal
  * then a #LinAlgPackIO::InputException# is thrown.  If #vs# is unsized then it is resized to #n#.
  * If #exta_flags & LinAlgPackIO::ignore_dim_bit == 0# then the number of elements read in depends
  * on the current size of #vs#.
  */
std::istream& input(std::istream& is, DVectorSlice* vs, LinAlgPackIO::fmtflags extra_flags);

// @}

}	// end namespace DenseLinAlgPack

#endif	// VECTOR_IN_FUNC_H
