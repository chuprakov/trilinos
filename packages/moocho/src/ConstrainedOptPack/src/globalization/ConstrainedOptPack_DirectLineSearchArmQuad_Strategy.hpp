// //////////////////////////////////////////////////////////////////////////////////
// DirectLineSearchArmQuad_Strategy.h
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

#ifndef DIRECT_LINE_SEARCH_ARM_QUAD_STRATEGY_H
#define DIRECT_LINE_SEARCH_ARM_QUAD_STRATEGY_H

#include "DirectLineSearch_Strategy.h"

namespace ConstrainedOptimizationPack {

///
/** Performs a line search using the Armijo condition and
  * uses quadratic interpolation to select each new alpha.
  */
class DirectLineSearchArmQuad_Strategy : public DirectLineSearch_Strategy {
public:

	/// Constructs with default settings.
	DirectLineSearchArmQuad_Strategy(int max_iter = 20, value_type eta = 1.0e-4
		, value_type min_frac = 0.1, value_type max_frac = 0.5 );

	/// Set the Armijo cord test fractional reduction parameter.
	void eta(value_type eta);
	///
	value_type eta() const;

	/// The minimum fraction that alpha is reduced for each line search iteration.
	void min_frac(value_type min_frac);
	///
	value_type min_frac() const;

	/// The maximum fraction that alpha is reduced for each line search iteration.
	void max_frac(value_type max_frac);
	///
	value_type max_frac() const;

	// ///////////////////////////////////////////////////
	// Overridden from DirectLineSearch_Strategy

	///
	void set_max_iter(int max_iter);
	///
	int max_iter() const;
	///
	int num_iterations() const;
	///
	/** Performs the following line search:\\
	  *
	  \begin{verbatim}
	  num_iter = 0;
	  while( phi.value(alpha_k) > phi_k + eta * alpha_k * phi.deriv() ) 
	  {
	     if(num_iter >= max_iter) return true;
	     num_iter = num_iter + 1;
	     alpha_k = [ min_frac * alpha_k <= quadradic interpolation for alpha
						<= max_frac * alpha_k ];
	  }
	  return true;\\
	  \end{verbatim}
	  *
	  * If the maximum number of iterations is exceeded then false will be returned.
	  *
	  * The default values for the adjustable parameters (from D&S A6.3.1)
	  * are:\\
	  * max_iter = 20\\
	  * eta = 1.0e-4\\
	  * min_frac = 0.1\\
	  * max_frac = 0.5\\
	  */
	bool do_line_search( const MeritFuncCalc1D& phi, value_type phi_k
		, value_type* alpha_k, value_type* phi_kp1
		, std::ostream* out);

	///
	void print_algorithm(std::ostream& out, const std::string& leading_str) const;

private:
	int	max_iter_;
	int	num_iter_;	// stores the number of interations
	value_type eta_, min_frac_, max_frac_;

	// Throw an exception if the parameters are not in a proper range.
	void validate_parameters() const;

};	// end class DirectLineSearchArmQuad_Strategy

// ///////////////////////////////////////////
// Inline member functions

inline
void DirectLineSearchArmQuad_Strategy::eta(value_type eta) {
	eta_ = eta;
}

inline
value_type DirectLineSearchArmQuad_Strategy::eta() const {
	return eta_;
}

inline
void DirectLineSearchArmQuad_Strategy::min_frac(value_type min_frac) {
	min_frac_ = min_frac;
}

inline
value_type DirectLineSearchArmQuad_Strategy::min_frac() const {
	return min_frac_;
}

inline
void DirectLineSearchArmQuad_Strategy::max_frac(value_type max_frac) {
	max_frac_ = max_frac;
}

inline
value_type DirectLineSearchArmQuad_Strategy::max_frac() const {
	return max_frac_;
}

}	// end namespace ConstrainedOptimizationPack

#endif	// DIRECT_LINE_SEARCH_ARM_QUAD_STRATEGY_H
