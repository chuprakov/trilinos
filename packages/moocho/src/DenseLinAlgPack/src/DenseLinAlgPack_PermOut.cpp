// /////////////////////////////////////////////////////////////////////////////
// PermOut.cpp

#include <iomanip>

#include "../include/PermOut.h"
#include "../include/IVector.h"

std::ostream& LinAlgPack::operator<<(std::ostream& o, const IVector& perm) {
	int w = o.width(0) - 1; // get the set width
	o << perm.size() << "\n";
	IVector::const_iterator
		itr_perm		= perm.begin(),
		itr_perm_end	= perm.end();
	for(;itr_perm != itr_perm_end;)
		o << std::setw(w) << *itr_perm++ << ' ';
	o << std::endl;
	return o;
}
