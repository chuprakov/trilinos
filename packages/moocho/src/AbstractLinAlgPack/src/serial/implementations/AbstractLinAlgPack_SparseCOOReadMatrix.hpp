// //////////////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_SparseCOOReadMatrix.hpp
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
// This is a function that reads a sparse matrix in form an imput stream.
//

#ifndef SPARSECOOREADMATRIX_H
#define SPARSECOOREADMATRIX_H

#include <istream>
#include <valarray>

#include "AbstractLinAlgPack_Types.hpp"

namespace AbstractLinAlgPack {

///
/** Read in a Coordinate Matrix from a C++ input stream and store it in valarrays.
  *
  * The format for the imput is:
  *
  * #m  n  nz#\\
  * #a1:i1:j1  a2:i2:j2 .... anz:inz:jnz#\\
  *
  * In the above format, each non-zero element is given as a three item pair:
  * value of the non-zero element, row indice (1-based) of the non-zero element,
  * and the column indice (1-based) of the non-zero element.  There must
  * be no spaces between the numbers and the \':\' charachter and there must
  * be at least one whitespace character between elements.
  *
  * The vectors #a#, #ivect#, and #jvect# are resized to #nz#.
  * If any problem is found in the input an InputException will be throw that
  * will include a descriptive message about the error.
  *
  * @param	m		number of rows of the sparse matrix
  * @param	n		number of columns of sparse matrix
  * @param	nz		number of non-zero elements of the sparse matrix
  * @param	a		vector holding the non-zero elements
  * @param	ivect	vector holding the row indices (1-based)
  * @param	jvect	vector holding the column indices (1-based)
  */
void read_coo_into_valarrays(std::istream& istrm, size_type& m, size_type& n, size_type& nz
	, std::valarray<value_type>& a, std::valarray<indice_type>& ivect
	, std::valarray<indice_type>& jvect);

}	// end namespace AbstractLinAlgPack 

#endif // SPARSECOOREADMATRIX_H
