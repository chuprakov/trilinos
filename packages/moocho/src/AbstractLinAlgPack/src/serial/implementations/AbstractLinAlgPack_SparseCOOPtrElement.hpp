// //////////////////////////////////////////////////////////////////////
// SparseCOOPtrElement.h

#ifndef SPARSE_COO_PTR_ELEMENT_H
#define SPARSE_COO_PTR_ELEMENT_H

#include "SparseLinAlgPackTypes.h"

namespace SparseLinAlgPack {

///
/** Sparse pointer element type for a COO matrix (val, ivect, jvect).
  *
  * This class abstracts a sparse element of a templated
  * type from a coordinate matrix. It
  * has a pointer to the value of the element.
  *
  * The default assignment operator and copy constructor
  * are allowed.
  */
template <class T_Indice, class T_Value>
class SparseCOOPtrElement {
public:
	/** @name Public Typedefs. */
	//@{

	///
	typedef T_Value							value_type;
	///
	typedef T_Indice						indice_type;

	//@}

	/** @name Constructors */
	//@{

	/// Construct uninitialized (poiner to value set to zero) (#indice() == 0#).
	SparseCOOPtrElement() : pvalue_(0), row_i_(0), col_j_(0)
	{}

	/// Construct with a pointer to the value and indice set
	SparseCOOPtrElement(value_type* pvalue, indice_type row_i, indice_type col_j)
		: pvalue_(pvalue), row_i_(row_i), col_j_(col_j)
	{}

	/// Initialize
	void initialize(value_type* pvalue, indice_type row_i, indice_type col_j) {
		pvalue_	= pvalue;
		row_i_	= row_i;
		col_j_	= col_j;
	}
	
	//@}

	/** @name Value and indice access */
	//@{ 

	///
	value_type& value()
	{
		return *pvalue_;
	}
	///
	value_type value() const
	{
		return *pvalue_;
	}
	///
	indice_type row_i() const
	{
		return row_i_;
	}
	///
	indice_type col_j() const
	{
		return col_j_;
	}
	/// Change the indices
	void change_indices(indice_type row_i, indice_type col_j)
	{
		row_i_ = row_i;
		col_j_ = col_j;
	}
	/// Change the element pointer
	void change_value_ptr(value_type* pvalue)
	{
		pvalue_ = pvalue;
	}

	//@}
private:
	value_type*				pvalue_;
	indice_type				row_i_, col_j_;

};	// end class SparseCOOPtrElement

} // end namespace SparseLinAlgPack 

#endif // SPARSE_COO_PTR_ELEMENT_H