// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncCalc1DQuadratic.h

#ifndef MERIT_FUNC_CALC_1D_QUADRATIC_H
#define MERIT_FUNC_CALC_1D_QUADRATIC_H

#include "MeritFuncCalc1D.h"
#include "LinAlgPack/include/VectorClass.h"

namespace ConstrainedOptimizationPack {

///
/** Adds the ability to compute phi(alpha) at alpha.
  *
  * Computes phi( x = d[0] + alpha * d[1] + alpha^2 * d[2] ).
  */
class MeritFuncCalc1DQuadratic : public MeritFuncCalc1D {
public:

	///
	/** The only constructor.
	  *
	  * Note that *x get updated as value(alpha) is called.
	  *
	  * The client must ensure that the memory pointed to by the
	  * vector slices in d must be desturbed while this object
	  * is in use.  To do so may have bad side effects.
	  */
	MeritFuncCalc1DQuadratic( const MeritFuncCalc& phi, size_type p, const VectorSlice d[]
		, VectorSlice* x );

	// ///////////////////////////////////////
	// Overridden from MeritFuncCalc1D

	/// Return phi( d[0] + alpha * d[1] + alpha^2 * d[2] + ... + alpha^p * d[p] ).
	value_type operator()(value_type alpha) const;

	/// Returns phi.deriv()
	value_type deriv() const;

	///
	void print_merit_func(std::ostream& out
		, const std::string& leading_str) const;

private:
	const MeritFuncCalc&		phi_;
	size_type					p_;
	const VectorSlice			d_[3];
	VectorSlice					*x_;

	// not defined and not to be called
	MeritFuncCalc1DQuadratic();
	MeritFuncCalc1DQuadratic& operator=( const MeritFuncCalc1DQuadratic& );


};	// end class MeritFuncCalc1DQuadratic

}	// end namespace ConstrainedOptimizationPack

#endif	// MERIT_FUNC_CALC_1D_QUADRATIC_H
