// //////////////////////////////////////////////////////////////////////////////////
// MatrixWithOpConcreteEncap.h

#ifndef MATRIX_WITH_OP_CONCRETE_ENCAP_H
#define MATRIX_WITH_OP_CONCRETE_ENCAP_H

#include "MatrixWithOp.h"

namespace SparseLinAlgPack {

///
/** This template class defines the storage for a concrete matrix
  * class that operations are based on.
  *
  * The default copy constructor and assignment operator are allowed.
  */
template<class M>
class MatrixWithOpConcreteEncap : public virtual MatrixWithOp
{
public:

	// /////////////////////////////////////////////////////
	/** @name Representation access */
	//@{

	/// The compiler did not generate this default constructor
	MatrixWithOpConcreteEncap()
	{}

	/// This constructor will have to be overridden.
	MatrixWithOpConcreteEncap(const M& m) : m_(m)
	{}

	/// Get the underlying M object
	M& m() {
		return m_;
	}

	///
	const M& m() const {
		return m_;
	}

	//@}	// end Representation access

	// /////////////////////////////////////////////////////
	// Overridden from Matrix

	///
	size_type rows() const;

	///
	size_type cols() const;

	// /////////////////////////////////////////////////////
	// Overridden from MatrixWithOp

	///
	MatrixWithOp& operator=(const MatrixWithOp& m);

private:
	M m_;

};	// end class MatrixWithOpConcreteEncap<M>

// Template definitions

template<class M>
size_type MatrixWithOpConcreteEncap<M>::rows() const {
	return m().rows();
}

template<class M>
size_type MatrixWithOpConcreteEncap<M>::cols() const {
	return m().cols();
}

template<class M>
MatrixWithOp& MatrixWithOpConcreteEncap<M>::operator=(const MatrixWithOp& m) {
	if(&m == this) return *this;	// assignment to self
	const MatrixWithOpConcreteEncap<M> *p_m = dynamic_cast<const MatrixWithOpConcreteEncap<M>*>(&m);
	if(p_m) {
		m_ = p_m->m_;
	}
	else {
		throw std::invalid_argument("MatrixWithOpConcreteEncap<M>::operator=(const MatrixWithOp& m)"
			" : The concrete type of m is not a subclass of MatrixWithOpConcreteEncap<M> as expected" );
	}
	return *this;
}

}	// end namespace SparseLinAlgPack 

#endif	// MATRIX_WITH_OP_CONCRETE_ENCAP_H
