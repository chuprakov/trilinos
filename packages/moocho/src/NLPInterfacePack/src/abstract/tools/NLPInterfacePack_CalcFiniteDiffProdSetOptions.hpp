// ////////////////////////////////////////////////////////////////
// NLPInterfacePack_CalcFiniteDiffProdSetOptions.hpp
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

#ifndef CALC_FINITE_DIFF_PROD_SET_OPTIONS_H
#define CALC_FINITE_DIFF_PROD_SET_OPTIONS_H

#include "NLPInterfacePack_CalcFiniteDiffProd.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace NLPInterfacePack {

///
/** Set options for \c CalcFiniteDiffProd from an
 * \c OptionsFromStream object.
 *
 * The default options group name is CalcFiniteDiffProd.
 *
 * The options group is:
 *
 \verbatim
 options_group CalcFiniteDiffProdSetOptions {
 *    fd_method_order = FD_ORDER_ONE;
 *    fd_method_order = FD_ORDER_TWO;
 *    fd_method_order = FD_ORDER_TWO_CENTRAL;
 *    fd_method_order = FD_ORDER_TWO_AUTO;
 *    fd_method_order = FD_ORDER_FOUR;
 *    fd_method_order = FD_ORDER_FOUR_CENTRAL;
 *    fd_method_order = FD_ORDER_FOUR_AUTO; *** (Default)
 *    fd_step_select = FD_STEP_ABSOLUTE; *** (Default)
 *    fd_step_select = FD_STEP_RELATIVE;
 *    fd_step_size = -1.0; *** (default)
 *    fd_step_size_min = -1.0; *** (default)
 *    fd_step_size_f = -1.0; *** (default)
 *    fd_step_size_c = -1.0; *** (default)
 }
 \endverbatim
 */
class CalcFiniteDiffProdSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			CalcFiniteDiffProd >
{
public:

	///
	CalcFiniteDiffProdSetOptions(
		CalcFiniteDiffProd* target = 0
		,const char opt_grp_name[] = "CalcFiniteDiffProd" );
	
protected:

	/// Overridden from SetOptionsFromStreamNode
	void setOption( int option_num, const std::string& option_value );

};	// end class CalcFiniteDiffProdSetOptions

}	// end namespace NLPInterfacePack

#endif	// CALC_FINITE_DIFF_PROD_SET_OPTIONS_H
