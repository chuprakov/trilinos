// /////////////////////////////////////////////////////////////////////////////////
// GenMatrixOutFunc.cpp

#include <ostream>
#include <iomanip>

#include "../include/GenMatrixOutFunc.h"
#include "../include/VectorOutFunc.h"
#include "../include/GenMatrixClass.h"

std::ostream& LinAlgPack::output(std::ostream& os, const GenMatrixSlice& gms
	, LinAlgPackIO::fmtflags extra_flags )
{
	int w = os.width(0) - 1; // get the set width (minus 1 since a space is inserted)

	if( !(extra_flags & LinAlgPackIO::ignore_dim_bit) ) {
		std::ios_base::fmtflags old = os.flags();
		os	<< std::setw(0) << std::left << gms.rows() << ' ' << gms.cols()
			<< std::endl;
		os.flags(old);
	}

	if( gms.rows() && gms.cols() ) {
		for(size_type i = 1; i <= gms.rows();++i) {
			const VectorSlice& vs =gms.row(i); 
			VectorSlice::const_iterator itr = vs.begin();
			for( size_type j = 1; itr != vs.end(); ++j, ++itr ) {
				os << " " << std::setw(w) << (*itr) << ":" << i << ":" << j;
			}
			os << std::endl;
		}
	}
	
	return os;
}
