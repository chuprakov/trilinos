// //////////////////////////////////////////////////////////////
// QPInitFixedFreeStd.h
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

#ifndef QP_INIT_FIXED_FREE_STD_H
#define QP_INIT_FIXED_FREE_STD_H

#include "QPSchur.h"

namespace ConstrainedOptimizationPack {
namespace QPSchurPack {

///
/** General (and flexible) implementation class for a QPSchur QP
  * problem.
  *
  * The basic idea of this class is to just build the QP from its
  * various components in a way that is easy and flexible for the
  * client.  The class will also do consistency testing if asked
  * to.
  */
class QPInitFixedFreeStd : public QP {
public:

	/// Construct uninitialized
	QPInitFixedFreeStd();

	///
	/** Initialize.
	  *
	  * The pointers and references to the objects pointed to by the
	  * arguments to this function must not be modified by the caller.
	  * Copies of these objects are not made internally so these
	  * objects must remain valid while this object is in use.
	  *
	  * If the sizes of the arguments do not match up or some consistency
	  * test fails then exceptions may be thrown with (hopefully) helpful
	  * messages.
	  *
	  *	@param	g 	[in] vector (size #n#): objective gradient
	  *	@param	G	[in] matrix (size #n x n#): objective Hessian
	  *	@param	A	[in] matrix (size #n x m#): full rank equality constraints
	  *					in #Ko#.  If #A==NULL# then there are no equality constraints
	  *					in #Ko# and m will be zero.
	  *	@param	n_R	[in] number of initially free variables
	  *	@param	i_x_free
	  *				[in] array (size #n_R#): #i_x_free[l-1], l = 1...n_R# defines
	  *					the matrix #Q_R# as:\\
	  *					#Q_R(:,l) = e(i_x_free[l-1]), l = 1...n_R#\\
	  *					The ordering of these indices is significant.  It is allowed
	  *                 for #i_x_free == NULL# in which case it will be
	  *                 considered to be identity.
	  *	@param	i_x_fixed
	  *				[in] array (size #n_X = n - n_R#):
	  *					#i_x_fixed[l-1], l = 1...n_X# defines the matrix #Q_X# as:\\
	  *					#Q_X(:,l) = e(i_x_fixed[l-1]), l = 1...n_X#\\
	  *					The ordering of these indices is significant.
	  *	@param	bnd_fixed
	  *				[in] array (size #n_X = n - n_R#):
	  *					#bnd_fixed[l-1], l = 1...n_X# defines the initial active set as:\\
	  \begin{verbatim}
	                    / LOWER : b_X(l) = xL(i_x_fixed[l-1])
	  bnd_fixed[l-1] = |  UPPER : b_X(l) = xU(i_x_fixed[l-1])
	                    \ EQUALITY : b_X(l) = xL(i) = xU(i) (i = i_x_fixed[l-1])
	  \end{verbatim}
	  *	@param	b_X	[in] vector (size #n_X = n - n_R#):
	  *				Initial varaible bounds (see #bnd_fixed#)
	  *	@param	Ko	[in] matrix (size #(n_R+m) x (n_R+m)#):  Initial KKT matrix
	  *	@param	fo 	[in] vector (size #n_R + m#): Initial KKT system rhs vector
	  *	@param	constraints
	  *				[in] Constraints object for the extra constraints
	  *					#cL_bar <= A_bar'*x <= cU_bar#
	  *	@param	out	[out] If #out!=NULL#, then any warning or error messages will
	  *					be printed here.
	  *	@param	test_setup
	  *				[in] If set to true, then consistency checks will be
	  *					made on all the input arguments.  The cost of the
	  *					tests will not be too excessive in runtime or
	  *					storge costs and do not completly validate everything
	  *	@param	waring_tol
	  *				[in] Warning tolerance for tests.
	  *	@param	error_tol
	  *				[in] Error tolerance for tests.  If the relative error
	  *					of any test exceeds this limit, then an error
	  *					message will be printed to out (if #out!=NULL#) and then
	  *					a runtime exception will be thrown.
	  *	@param	print_all_warnings
	  *				[in] If set to #true#, then any relative errors for tests
	  *					that are above #warning_tol# will be printed to
	  *					#out# (if #out!= NULL#) (O(#n#) output).
	  *					Otherwise, if #false#, then
	  *					only the number of violations and the maximum
	  *					violation will be printed (O(1) output).
	  */
	void initialize(
		  const VectorSlice						&g
		, const MatrixSymWithOp					&G
		, const MatrixWithOp					*A
		, size_type								n_R
		, const size_type						i_x_free[]
		, const size_type						i_x_fixed[]
		, const EBounds							bnd_fixed[]
		, const VectorSlice						&b_X
		, const MatrixSymWithOpFactorized		&Ko
		, const VectorSlice						&fo
		, Constraints							*constraints
		, std::ostream							*out				= NULL
		, bool									test_setup			= false
		, value_type							warning_tol			= 1e-10
		, value_type							error_tol			= 1e-5
		, bool									print_all_warnings	= false
		);

	// /////////////////////
	// Overridden from QP 

	///
	size_type n() const;
	///
	size_type m() const;
	///
	const VectorSlice g() const;
	///
	const MatrixSymWithOp& G() const;
	///
	const MatrixWithOp& A() const;
	///
	size_type n_R() const;
	///
	const x_init_t& x_init() const;
	///
	const l_x_X_map_t& l_x_X_map() const;
	///
	const i_x_X_map_t& i_x_X_map() const;
	///
	const VectorSlice b_X() const;
	///
	const GenPermMatrixSlice& Q_R() const;
	///
	const GenPermMatrixSlice& Q_X() const;
	///
	const MatrixSymWithOpFactorized& Ko() const;
	///
	const VectorSlice fo() const;
	///
	Constraints& constraints();
	///
	const Constraints& constraints() const;

private:

	// ///////////////////////////////////
	// Private types

	typedef std::vector<size_type>		row_i_t;
	typedef std::vector<size_type>		col_j_t;
	
	// ///////////////////////////////////
	// Private data members


	size_type				n_;
	size_type				n_R_;
	size_type				m_;
	VectorSlice				g_;	// will not be modified!
	const MatrixSymWithOp	*G_;
	const MatrixWithOp		*A_;	// If NULL not no equalities in Ko
	x_init_t				x_init_;
	l_x_X_map_t				l_x_X_map_;
	i_x_X_map_t				i_x_X_map_;
	VectorSlice				b_X_;	// will not be modified!
	GenPermMatrixSlice		Q_R_;
	row_i_t					Q_R_row_i_;
	col_j_t					Q_R_col_j_;
	GenPermMatrixSlice		Q_X_;
	row_i_t					Q_X_row_i_;
	col_j_t					Q_X_col_j_;
	const MatrixSymWithOpFactorized
							*Ko_;
	VectorSlice				fo_;	// will not be modified
	Constraints				*constraints_;

	// Private member function
	void assert_initialized() const;

};	// end class QPInitFixedFreeStd

}	// end namespace QPSchurPack
}	// end namespace ConstrainedOptimizationPack 

#endif	// QP_INIT_FIXED_FREE_STD_H
