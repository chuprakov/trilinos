// ///////////////////////////////////////////////////////////
// random_vector.cpp

#include <stdlib.h>

#include <stdexcept>

#include "LinAlgPack/include/random_vector.h"
#include "LinAlgPack/include/VectorClass.h"

void LinAlgPack::seed_random_vector_generator( unsigned int s )
{
	srand(s);
}

void LinAlgPack::random_vector( value_type l, value_type u, VectorSlice* v )
{
	if(!v)
		throw std::invalid_argument( "random_vector(...) : Error, "
			"v can not be NULL" );
	if( l > u )
		throw std::invalid_argument( "random_vector(...) : Error, "
			"l can not be greater than u" );
	for( VectorSlice::iterator itr = v->begin(); itr != v->end(); )
		*itr++ = l + (double(rand())/RAND_MAX) * (u -l);
}
