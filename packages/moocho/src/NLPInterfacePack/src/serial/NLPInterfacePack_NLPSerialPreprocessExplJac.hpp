// ///////////////////////////////////////////////////////////////////////
// NLPSerialPreprocessExplJac.h
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#ifndef NLP_SERIAL_PREPROCESS_EXPL_JAC_H
#define NLP_SERIAL_PREPROCESS_EXPL_JAC_H

#include <valarray>

#include "NLPSerialPreprocess.h"
#include "NLPFirstOrderInfo.h"
#include "LinAlgPack/include/VectorClass.h"
#include "AbstractFactory.h"

namespace NLPInterfacePack {

///
/** NLP node subclass for explicit Jacobians.
 *
 * ToDo: Finish documentation!
 */
class NLPSerialPreprocessExplJac
	: virtual public NLPSerialPreprocess
	, virtual public NLPFirstOrderInfo
{
public:

	/** @name Public types */
	//@{
	
	///
	typedef MemMngPack::ref_count_ptr<
		const MemMngPack::AbstractFactory<MatrixWithOp> >    factory_mat_ptr_t;

	//@}

	/** @name Constructors / initializers */
	//@{

	///
	/** Calls <tt>this->set_mat_factories()</tt>.
	 */
	NLPSerialPreprocessExplJac(
		const factory_mat_ptr_t     &factory_Gc_orig = MemMngPack::null
		,const factory_mat_ptr_t    &factory_Gh_orig = MemMngPack::null
		);

	///
	/** Initialize with matrix factories for \c Gc and \c Gh.
	 *
	 * @param  factory_Gc_orig
	 *                [in] Smart pointer to matrix factory for \c Gc_orig.  If
	 *                <tt>factory_Gc_orig.get() == NULL</tt> then the concrete matrix
	 *                type ??? will be used as the default.
	 * @param  factory_Gh_orig
	 *                [in] Smart pointer to matrix factory for \c Gh_orig.  If
	 *                <tt>factory_Gh_orig.get() == NULL</tt> then the concrete matrix
	 *                type ??? will be used as the default. */
	void set_mat_factories(
		const factory_mat_ptr_t     &factory_Gc_orig
		,const factory_mat_ptr_t    &factory_Gh_orig
		);

	//@}

	/** @name Overridden public members from NLP */
	//@{

	///
	void initialize();	
	///
	bool is_initialized() const;

	//@}

	/** @name Overridden public members from NLPFirstOrderInfo */
	//@{
	
	///
	const mat_fcty_ptr_t factory_Gc() const;
	///
	const mat_fcty_ptr_t factory_Gh() const;
	/// Validates the type of Gc is correct
	void set_Gc(MatrixWithOp* Gc);
	/// Validates the type of Gh is correct
	void set_Gh(MatrixWithOp* Gh);

	//@}

	/** @name Overridden public members from NLPVarReductPerm */
	//@{

	///
	bool get_next_basis(
		Permutation*  P_var,   Range1D* var_dep
		,Permutation* P_equ,   Range1D* equ_decomp
		,Permutation* P_inequ, Range1D* inequ_decomp
		);
	///
	void set_basis(
		const Permutation   &P_var,   const Range1D  &var_dep
		,const Permutation  *P_equ,   const Range1D  *equ_decomp
		,const Permutation  *P_inequ, const Range1D  *inequ_decomp
		);

	//@}


protected:

	/** @name Overridden protected members from NLPFirstOrderInfo */
	//@{

	///
	void imp_calc_Gc(
		const VectorWithOp& x, bool newx
		,const FirstOrderInfo& first_order_info
		) const;
	///
	void imp_calc_Gh(
		const VectorWithOp& x, bool newx
		,const FirstOrderInfo& first_order_info
		) const;
	
	//@}

	/** @name Protected types */
	//@{

	///
	/** Struct for zero and explicit first order quantities that subclass must fill in.
	 */
	struct FirstOrderExplInfo {
		///
		typedef std::valarray<value_type>    val_full_t;
		///
		typedef std::valarray<index_type>    ivect_full_t;
		//
		typedef std::valarray<index_type>    jvect_full_t;
		///
		FirstOrderExplInfo()
			:Gc_val(NULL), Gc_ivect(NULL), Gc_jvect(NULL)
			,Gh_val(NULL), Gh_ivect(NULL), Gh_jvect(NULL)
			,f(NULL)
		{}
		///
		FirstOrderExplInfo(
			index_type* Gc_nz_in, val_full_t* Gc_val_in, ivect_full_t* Gc_ivect_in, jvect_full_t* Gc_jvect_in
			,index_type* Gh_nz_in, val_full_t* Gh_val_in, ivect_full_t* Gh_ivect_in, jvect_full_t* Gh_jvect_in
			,const ObjGradInfoSerial& obj_grad
			)
			:Gc_nz(Gc_nz_in), Gc_val(Gc_val_in), Gc_ivect(Gc_ivect_in), Gc_jvect(Gc_jvect_in)
			,Gh_nz(Gh_nz_in), Gh_val(Gh_val_in), Gh_ivect(Gh_ivect_in), Gh_jvect(Gh_jvect_in)
			,Gf(obj_grad.Gf), f(obj_grad.f), c(obj_grad.c), h(obj_grad.h)
		{}
		///
		size_type*    Gc_nz;
		///
		val_full_t*   Gc_val;
		///
		ivect_full_t* Gc_ivect;
		///
		jvect_full_t* Gc_jvect;
		///
		size_type*    Gh_nz;
		///
		val_full_t*   Gh_val;
		///
		ivect_full_t* Gh_ivect;
		///
		jvect_full_t* Gh_jvect;
		///
		VectorSlice   Gf;
		///
		value_type*   f;
		///
		VectorSlice   c;
		///
		VectorSlice   h;
	}; // end struct FirstOrderExplInfo

	//@}

	/** @name Pure virtual template methods to be defined by subclasses */
	//@{

	///
	/** Return the number of nonzero elements in \c Gc before elements are removed for fixed variables.
	  *
	  * The value returned from this method before the first time \c imp_calc_Gc() is called
	  * is an upper estimate of the number of nonzeros.  To get the actual number
	  * of nonzeros, call this function again after \c imp_calc_Gc() has been called.
	  */
	virtual size_type imp_Gc_nz_full() const = 0;

	///
	/** Return the number of nonzero elements in \c Gh before elements are removed for fixed variables.
	  *
	  * The value returned from this method before the first time \c imp_calc_Gh() is called
	  * is an upper estimate of the number of nonzeros.  To get the actual number
	  * of nonzeros, call this function again after \c imp_calc_Gh() has been called.
	  */
	virtual size_type imp_Gh_nz_full() const = 0;

	///
	/** Calculate the COOR matrix for the gradient for all of the \a c(x) constaints in the full %NLP.
	 *
	 * ToDo: Finish documentation!
	 */
	virtual void imp_calc_Gc_full(
		const VectorSlice& x_full, bool newx
		, const FirstOrderExplInfo& first_order_expl_info
		) const = 0;

	///
	/** Calculate the COOR matrix for the gradient for all of the \c h(x) constaints in the full %NLP.
	 *
	 * ToDo: Finish documentation!
	 */
	virtual void imp_calc_Gh_full(
		const VectorSlice& x_full, bool newx
		, const FirstOrderExplInfo& first_order_expl_info
		) const = 0;

	//@}

	/** @name Protected member functions for subclasses to use */
	//@{

	/// Assert if we have been initizlized (throws UnInitialized)
	void assert_initialized() const;

	///
	const FirstOrderExplInfo first_order_expl_info() const;

	//@}

private:

	// ////////////////////////////////////////
	// Private data members
	
	bool initialized_;              // Flag for if the NLP has has been properly initialized

	factory_mat_ptr_t   factory_Gc_orig_;
	factory_mat_ptr_t   factory_Gh_orig_;
	mat_fcty_ptr_t      factory_Gc_;
	mat_fcty_ptr_t      factory_Gh_;

//	mutable size_type	Gc_nz_;        // Number of nonzeros in the shrunken NLP Gc
//	mutable size_type	Gh_nz_;        // Number of nonzeros in the shrunken NLP Gh
	mutable size_type                         Gc_nz_full_;    // Number of nonzeros in the full NLP Gc
	mutable FirstOrderExplInfo::val_full_t    Gc_val_full_;   // Storage for explicit nonzeros of full Gc
	mutable FirstOrderExplInfo::ivect_full_t  Gc_ivect_full_;
	mutable FirstOrderExplInfo::jvect_full_t  Gc_jvect_full_;
	mutable size_type                         Gh_nz_full_;    // Number of nonzeros in the full NLP Gc
	mutable FirstOrderExplInfo::val_full_t    Gh_val_full_;   // Storage for explicit nonzeros of full Gh
	mutable FirstOrderExplInfo::ivect_full_t  Gh_ivect_full_;
	mutable FirstOrderExplInfo::jvect_full_t  Gh_jvect_full_;

	mutable bool                              Gc_perm_new_basis_updated_;      // Flag for if a new basis was set!
	mutable bool                              Gh_perm_new_basis_updated_;      // Flag for if a new basis was set!

	// ////////////////////////////
	// Private member functions

	//
	void imp_calc_Gc_or_Gh(
		bool calc_Gc
		,const VectorWithOp& x, bool newx
		,const FirstOrderInfo& first_order_info
		) const;

};	// end class NLPSerialPreprocessExplJac

// ///////////////////////////
// inline members

inline
const NLPSerialPreprocessExplJac::FirstOrderExplInfo
NLPSerialPreprocessExplJac::first_order_expl_info() const
{
	return FirstOrderExplInfo(
		&Gc_nz_full_,&Gc_val_full_,&Gc_ivect_full_,&Gc_jvect_full_
		,&Gh_nz_full_,&Gh_val_full_,&Gh_ivect_full_,&Gh_jvect_full_
		,obj_grad_full_info()
		);
}

}	// end namespace NLPInterfacePack 

#endif // NLP_SERIAL_PREPROCESS_EXPL_JAC_H
