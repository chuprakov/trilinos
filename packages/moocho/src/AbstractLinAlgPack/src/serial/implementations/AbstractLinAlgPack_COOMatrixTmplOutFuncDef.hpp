// /////////////////////////////////////////////////////////////////////////////////
// COOMatrixTmplOutFuncDef.h

#ifndef COO_MATRIX_TMPL_OUT_FUNC_DEF_H
#define COO_MATRIX_TMPL_OUT_FUNC_DEF_H

#include <ostream>
#include <iomanip>

#include "COOMatrixTmplOutFuncDecl.h"

namespace SparseLinAlgPack {

template <class T_COOM>
std::ostream& output_COOM(std::ostream& os, const T_COOM& coom
	, SparseLinAlgPackIO::fmtflags extra_flags)
{
	int w = os.width(0) - 1; // get the set width (minus 1 since a space is inserted)

	if(    !(extra_flags & SparseLinAlgPackIO::ignore_dim_bit)
		|| !(extra_flags & SparseLinAlgPackIO::ignore_nz_bit) )
	{

		os << std::setw(0) << std::left;

		if( !(extra_flags & SparseLinAlgPackIO::ignore_dim_bit) )
			os << coom.rows() << ' ' << coom.cols() << ' ';

		if( !(extra_flags & SparseLinAlgPackIO::ignore_nz_bit) )
			os << coom.nz();

		os << std::endl << std::right;
	}
	
	if(!coom.nz()) return os;	// no elements to output
	
	T_COOM::difference_type
		row_offset	= coom.row_offset(),
		col_offset	= coom.col_offset();
	for(T_COOM::const_iterator itr = coom.begin(); itr != coom.end();++itr) {
		os	<< " " << std::setw(w) << itr->value()
			<< ':' << itr->row_i() + row_offset
			<< ':' << itr->col_j() + col_offset;
	}

	if( !(extra_flags & SparseLinAlgPackIO::no_insert_newlines_bit) )
		os << std::endl;

	return os;
}

}	// end namespace SparseLinAlgPack 

#endif	// COO_MATRIX_TMPL_OUT_FUNC_DEF_H