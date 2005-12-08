// /////////////////////////////////////////////////////////////////////////////////
// DenseLinAlgPack_DMatrixOutFunc.hpp
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

#ifndef GENMATRIX_OUT_FUNC_H
#define GENMATRIX_OUT_FUNC_H

#include "DenseLinAlgPack_IOBasic.hpp"

namespace DenseLinAlgPack {

/* * @name DMatrixSlice output stream function.
  *
  * This is a function that are used to output a DMatrixSlice object
  * to a char based output stream.  This format is ment to be both
  * human and machine readable.  In fact the \Ref{input} function
  * can be used to read in the output produced by this function.
  * This function is ment to be the machinary for an optional output stream
  * operator \Ref{operator<<}.
  *
  * The output format is diferent depending on the on whether the
  * bits #LinAlgPackIO::ignore_dim_bit# and #LinAlgPackIO::no_insert_newlines_bit# are set.
  * The default output format is (exta_flags == 0):
  *
  * #gms.rows()  gms.cols()#\\
  * #gms(1,1)           gms(1,2)          gms(1,3)          ...      gms(1,gms.cols())#\\
  * #gms(2,1)           gms(2,2)          gms(2,3)          ...      gms(2,gms.cols())#\\
  * #   .                   .                .                               .#\\
  * #gms(gms.rows(),1) gms(gms.rows(),2)  gms(gms.rows(),3)          gms(gms.rows(),gms.cols())#\\
  *
  * If #extra_flags & LinAlgPackIO::ignore_dim_bit != 0# then the dimensions of the
  * matrix are not output and if #extra_flags & LinAlgPackIO::no_insert_newlines_bit != 0#
  * then newline char ('\n') are not inserted after each row.
  *
  * Each of the elements are are put into columns according to the width set in the 
  * output stream #os# when it is called and other formating commands.  Even if the set
  * width is 0 or less than
  * the number of char's for the element a space ' ' will be inserted between them.
  * The elements are formated according to the format in the stream #os#.
  *
  * If #gms.rows() == 0# then no elements will be output and if
  * #extra_flags & LinAlgPackIO::ignore_dim_bit != 0# and only the operation
  * #os << gms.rows() << ' ' << gms.cols() << endl;# will be performed.
  *
  * If any of the output operations fails then a #std::ios_base::failure# exception is thrown. 
  */
std::ostream& output(std::ostream& os, const DMatrixSlice& gms, LinAlgPackIO::fmtflags extra_flags);

}	// end namespace DenseLinAlgPack

#endif	// GENMATRIX_OUT_FUNC_H
