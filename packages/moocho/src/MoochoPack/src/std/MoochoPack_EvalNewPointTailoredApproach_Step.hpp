// ////////////////////////////////////////////////////////////////////////////
// EvalNewPointTailoredApproach_Step.hpp
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

#ifndef EVAL_NEW_POINT_TAILORED_APPROACH_STEP_H
#define EVAL_NEW_POINT_TAILORED_APPROACH_STEP_H

#include "MoochoPack/src/MoochoPackTypes.hpp"
#include "IterationPack/src/AlgorithmStep.hpp"
#include "ConstrainedOptPack/src/misc/VariableBoundsTester.hpp"
#include "NLPInterfacePack/src/abstract/test/NLPDirectTester.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

///
/** Base class for evaluating a new point for the "Tailored Approach".
 *
 * Uses the \c NLPDirect interface to compute <tt>Z = [ -inv(C)*N; I ]</tt>
 * and <tt>py = -inv(C)*c(decomp)</tt> explicitly.  Subclasses determine how
 * <tt>py</tt> and <tt>Y</tt> are updated.
 */
class EvalNewPointTailoredApproach_Step
	: public IterationPack::AlgorithmStep // doxygen needs full path
{
public:

	/** @name Public types */
	//@{
	///
	enum EFDDerivTesting { FD_DEFAULT, FD_TEST, FD_NO_TEST };
	//@}

	/** @name Constructors / initializers */
	//@{

	///
	typedef Teuchos::RefCountPtr<const MatrixOp> D_ptr_t;
	/// <<std comp>> members for testing object for NLPDirect
	STANDARD_COMPOSITION_MEMBERS( NLPDirectTester, deriv_tester )
	/// <<std comp>> Members for variable bounds tester object
	STANDARD_COMPOSITION_MEMBERS( VariableBoundsTester, bounds_tester )
	///
	/** Set how and if finite derivatives are tested.
	 *
	 * ToDo: Finish documentation.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EFDDerivTesting, fd_deriv_testing )

	///
	EvalNewPointTailoredApproach_Step(
		  const deriv_tester_ptr_t& 	deriv_tester
		, const bounds_tester_ptr_t&	bounds_tester
		, EFDDerivTesting				fd_deriv_testing = FD_DEFAULT
		);

	//@}

	/** @name Overridden from AlgorithmStep */
	//@{
	///
	bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss);
	///
	void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
	//@}

	/** @name To be overridden by subclasses */
	//@{

	///
	/** Call to uninitialize the matrices.
	 *
	 * ToDo: Finish documentation!
	 */
	virtual void uninitialize_Y_Uy(
		MatrixOp         *Y
		,MatrixOp        *Uy
		) = 0;

	///
	/** Overridden by subclass to compute \c py, \c Y and \c Uy.
	 *
	 * @param  D    [in/out] Smart pointer to matrix <tt>D = -inv(C)*N</tt>.
	 *              On output, D->count() may be incremented in order to
	 *              initialize <tt>Y</tt>.
	 * @param  py   [in/out] On input <tt>py = -inv(C)*c(decomp)</tt>.
	 *              On output <tt>py = -inv((Gc(decomp)'*Y)*c(decomp)</tt>
	 * @param  Y    [in/out] On ouput <tt>Y</tt> is initialized properly.
	 * @param  Uy   [in/out] On ouput <tt>Uy</tt> is initialized properly.
	 * @param  olevel
	 *              [in] Determines output level.
	 * @param  out  [out] Journal outputting.
	 */
	virtual void calc_py_Y_Uy(
		const NLPDirect       &nlp
		,const D_ptr_t        &D
		,VectorMutable        *py
		,MatrixOp             *Y
		,MatrixOp             *Uy
		,EJournalOutputLevel  olevel
		,std::ostream         &out
		) = 0;

	///
	/** Overridden by subclass to recompute \c py and \c Ypy.
	 *
	 * @param  D    [in] matrix <tt>D = -inv(C)*N</tt>
	 * @param  py   [in/out] On input <tt>py = -inv(C)*c(decomp)</tt>.
	 *              On output <tt>py = -inv((Gc(decomp)'*Y)*c(decomp)</tt>
	 */
	virtual void recalc_py(
		const MatrixOp          &D
		,VectorMutable          *py
		,EJournalOutputLevel    olevel
		,std::ostream           &out
		) = 0;
	
	///
	/** Overridden by subclass to print how \c py and \c Y are computed.
	 */
	virtual void print_calc_py_Y_Uy(
		std::ostream& out, const std::string& leading_str
		) const = 0;

	//@}

private:

	// Not defined and not to be called
	EvalNewPointTailoredApproach_Step();

};	// end class EvalNewPointTailoredApproach_Step

}	// end namespace MoochoPack 

#endif	// EVAL_NEW_POINT_TAILORED_APPROACH_STEP_H
