// /////////////////////////////////////////////////////////////////////////////
// PermIn.cpp

#include "../include/PermIn.h"
#include "../include/IVector.h"

std::istream& LinAlgPack::operator>>(std::istream& istrm, IVector& perm) {
	size_type size;
	istrm >> size;
	perm.resize(size);

	IVector::iterator		itr_perm		= perm.begin(),
							itr_perm_end	= perm.end();
	for(;itr_perm != itr_perm_end;)
		istrm >> *itr_perm++;
	return istrm;
}
