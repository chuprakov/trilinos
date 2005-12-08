// ///////////////////////////////////////////////////////////
// DenseLinAlgPack_random_vector.hpp
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

#ifndef RANDOM_VECTOR_H
#define RANDOM_VECTOR_H

#include "DenseLinAlgPack_Types.hpp"

namespace DenseLinAlgPack {

///
/* * Seed the random number generator
  */
void seed_random_vector_generator( unsigned int );

///
/* * Generate a random vector with elements uniformly
  * distrubuted elements.
  * 
  * The elements are randomly generated between
  * [l, u].
  */
void random_vector( value_type l, value_type u, DVectorSlice* v );

}	// end namespace DenseLinAlgPack

#endif    // RANDOM_VECTOR_H
