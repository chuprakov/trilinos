// /////////////////////////////////////////////////////
// MatlabPack.cpp

#include <limits>
#include <iomanip>

#include "../include/MatlabPack.h"
#include "../include/GenMatrixClass.h"

namespace LinAlgPack {

std::ostream& MatlabPack::out( std::ostream& o, const char* name, const VectorSlice& vs
	, BLAS_Cpp::Transp trans )
{
	int p = o.precision();
	o.precision( std::numeric_limits<value_type>::digits10 + 3 );
	try {
		o << name << " =  [ ";
		for( VectorSlice::const_iterator itr = vs.begin(); itr != vs.end(); ++itr )
			o << ' ' << *itr << ";";
		o << "];" << (trans == BLAS_Cpp::no_trans ? ' ' : '\'' ) << std::endl;
	}
	catch(...) {
		o.precision(p);
		throw;
	}
	o.precision(p);
	return o;
}

std::ostream& MatlabPack::out( std::ostream& o, const char* name, const GenMatrixSlice& gms
	, BLAS_Cpp::Transp trans )
{
	int p = o.precision();
	o.precision( std::numeric_limits<value_type>::digits10 + 3 );
	try {
		o << name << " =  [\n";
		for( size_type i = 1; i <= gms.rows(); ++i ) {
			const VectorSlice vs = gms.row(i);
			for( VectorSlice::const_iterator itr = vs.begin(); itr != vs.end(); ++itr )
				o << *itr << ", ";
			o << ";\n";
		}
		o << "];" << (trans == BLAS_Cpp::no_trans ? ' ' : '\'' ) << std::endl;
	}
	catch(...) {
		o.precision(p);
		throw;
	}
	o.precision(p);
	return o;
}

}	// end namespace LinAlgPack