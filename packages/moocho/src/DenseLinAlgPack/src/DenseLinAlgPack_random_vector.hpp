// ///////////////////////////////////////////////////////////
// random_vector.h

#ifndef RANDOM_VECTOR_H
#define RANDOM_VECTOR_H

#include "LinAlgPackTypes.h"

namespace LinAlgPack {

///
/** Seed the random number generator
  */
void seed_random_vector_generator( unsigned int );

///
/** Generate a random vector with elements uniformly
  * distrubuted elements.
  * 
  * The elements are randomly generated between
  * [l, u].
  */
void random_vector( value_type l, value_type u, VectorSlice* v );

}	// end namespace LinAlgPack

#endif    // RANDOM_VECTOR_H