// //////////////////////////////////////////////////////////////////////////////
// COOMatrixOutFunc.cpp

#include "../include/COOMatrixOutFunc.h"
#include "../include/COOMatrixClass.h"

std::ostream& SparseLinAlgPack::output(std::ostream& o, const COOMatrix& coom) {

	o	<< coom.rows() << " " << coom.cols() << " " << coom.nz() << "\n";

	const COOMatrix::value_type
		*itr_val		= coom.const_val(),
		*itr_val_end	= coom.const_val() + coom.nz();
	const COOMatrix::indice_type
		*itr_ivect	= coom.const_ivect(),
		*itr_jvect	= coom.const_jvect();

	for(; itr_val != itr_val_end; ++itr_val, ++itr_ivect, ++itr_jvect)
		o << *itr_val << ":" << *itr_ivect << ":" << *itr_jvect << " ";
	o << "\n";

	return o;
}
