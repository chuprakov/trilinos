// /////////////////////////////////////////////////////////////////////////////
// SparseCOOReadMatrix.cpp

#include <stdlib.h>

#include "../include/SparseCOOReadMatrix.h"

// Throw an exception if the char is not ':'
namespace {
inline void assert_sep_char(char c) {
	if(c != ':')
		throw SparseLinAlgPack::InputException("Sparse COO matrix input stream error:  The seperator between the element, row indice and column indice must be a \':\'");
}
inline void assert_eof(std::istream& istrm) {
	if(istrm.eof())
		throw SparseLinAlgPack::InputException("Sparse COO matrix input stream error:  Premature end to the input file.");
}
}

void SparseLinAlgPack::read_coo_into_valarrays(std::istream& istrm, size_type& m, size_type& n, size_type& nz
	, std::valarray<value_type>& a, std::valarray<indice_type>& ivect
	, std::valarray<indice_type>& jvect)
{
	// read in dimensions and resize
	istrm >> m;		assert_eof(istrm);
	istrm >> n;		assert_eof(istrm);
	istrm >> nz;
	a.resize(nz);
	ivect.resize(nz);
	jvect.resize(nz);

	// Read in the non-zero elements
	value_type	*p_a =			&a[0],
				*p_a_last =		p_a + nz;
	indice_type	*p_ivect =		&ivect[0],
				*p_jvect =		&jvect[0];

	for(; p_a != p_a_last; ++p_a, ++p_ivect, ++p_jvect) {
		const int bs = 50;
		char num[bs];
		char c;
		assert_eof(istrm);
		istrm.get(num, bs-1, ':');  assert_eof(istrm);  *p_a = ::atof(num);	// Read in ak 
		istrm.get(c); assert_eof(istrm); assert_sep_char(c);				// Read in ':'
		istrm.get(num, bs-1, ':'); assert_eof(istrm); *p_ivect = ::atoi(num);// Read in ik 
		istrm.get(c); assert_eof(istrm); assert_sep_char(c);				// Read in ':'
		istrm >> *p_jvect;										// Read in jk
	}
}
