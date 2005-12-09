// ////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_DirectSparseSolverMA28SetOptions.hpp
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

#ifdef SPARSE_SOLVER_PACK_USE_MA28

#ifndef DIRECT_SPARSE_SOLVER_MA28_SET_OPTIONS_H
#define DIRECT_SPARSE_SOLVER_MA28_SET_OPTIONS_H

#include "AbstractLinAlgPack_DirectSparseSolverMA28.hpp"
#include "OptionsFromStreamPack_SetOptionsFromStreamNode.hpp"
#include "OptionsFromStreamPack_SetOptionsToTargetBase.hpp"

namespace AbstractLinAlgPack {

///
/** Set options for DirectSparseSolverMA28 from
  * OptionsFromStream object.
  *
  * The options group is:
  *
  \begin{verbatim}
  options_group DirectSparseSolverMA28 {
      estimated_fillin_ratio = 10.0;
      u = 0.1;
      grow = true;
      tol = 0.0;
      nsrch = 4;
      lbig = false;
      print_ma28_outputs = false;
      output_file_name = NONE;
  }
  \end{verbatim}
  *
  * See MA28 documentation for a description of these options.
  */
class DirectSparseSolverMA28SetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			DirectSparseSolverMA28 >
{
public:

	///
	DirectSparseSolverMA28SetOptions(
		DirectSparseSolverMA28* qp_solver = 0 );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void setOption( int option_num, const std::string& option_value );

};	// end class DirectSparseSolverMA28SetOptions

}	// end namespace AbstractLinAlgPack 

#endif	// DIRECT_SPARSE_SOLVER_MA28_SET_OPTIONS_H

#endif // SPARSE_SOLVER_PACK_USE_MA28
