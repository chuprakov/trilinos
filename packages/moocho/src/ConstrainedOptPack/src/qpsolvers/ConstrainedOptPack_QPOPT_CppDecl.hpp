// ////////////////////////////////////////////////////////////////////////////////
// QPOPT_CppDecl.h
//
// C declarations for QPSOL functions.  These declarations should not have to change
// for different platforms.  As long as the fortran object code uses capitalized
// names for its identifers then the declarations in fortran_types.h should be
// sufficent for portability.

#ifndef QPOPT_CPP_DECL_H
#define QPOPT_CPP_DECL_H

#include "Misc/include/fortran_types.h"

namespace QPOPT_CppDecl {

// Declarations that will link to the fortran object file.
// These may change for different platforms

using FortranTypes::f_int;			// INTEGER
using FortranTypes::f_real;			// REAL
using FortranTypes::f_dbl_prec;		// DOUBLE PRECISION
using FortranTypes::f_logical;		// LOGICAL

// /////////////////////////////////////////////////////////////////////////
/** @name QPOPT interface functions.
  */

//@{

typedef FORTRAN_FUNC_PTR_DECL(void,qphess_func) ( const f_int& N, const f_int& LDH
	, const f_int& JTHCOL, const f_dbl_prec* HESS, const f_dbl_prec* X, f_dbl_prec* HX
	, f_int* IW, const f_int& LENIW, f_dbl_prec* W, const f_int& LENW );

/// Call QPOPT through C++ declaration.
void qpopt( const f_int& N, const f_int& NCLIN
	, const f_int& LDA, const f_int& LDH, const f_dbl_prec* A
	, const f_dbl_prec* BL, const f_dbl_prec* BU, const f_dbl_prec* CVEC
	, const f_dbl_prec* H, qphess_func QPHESS, f_int* ISTATE, f_dbl_prec* X
	, f_int& INFORM, f_int& ITER, f_dbl_prec& OBJ, f_dbl_prec* AX
	, f_dbl_prec* CLAMDA, f_int* IW, const f_int& LENIW, f_dbl_prec* W
	, const f_int& LENW );

// //////////////////////////////////////////////////////////
// Enumerations for QPOPT options

/// Use to specify problem types
enum EQPOPT_problem_type {
		FP								= 1,
		LP								= 2,
		QP1								= 3,
		QP2								= 4,
		QP3								= 5,
		QP4								= 6
};

/// Integer options
enum EQPOPT_int_option {
		CHECK_FREQUENCY					= 1,
		EXPAND_FREQUENCY				= 2,
		FEASIBILITY_PHASE_ITER_LIMIT	= 3,
		OPTIMALITY_PHASE_ITER_LIMIT		= 4,
		HESSIAN_ROWS					= 5,
		ITERATION_LIMIT					= 6,
		MAXIMUM_DEGREES_OF_FREEDOM		= 7,
		PRINT_FILE						= 8,
		PRINT_LEVEL						= 9,
		PROBLEM_TYPE					= 10,	// Use EQPOPT_problem_type
		SUMMARY_FILE					= 11
};


/// Logical options
enum EQPOPT_logical_option {
		WARM_START						= 1,
		LIST							= 2,
		MIN_SUM							= 3
};

/// Real options
enum EQPOPT_real_option {
		CRASH_TOLERANCE					= 1,
		FEASIBILITY_TOLERANCE			= 2,
		INFINITE_BOUND_SIZE				= 3,
		INFINITE_STEP_SIZE				= 4,
		OPTIMALITY_TOLERANCE			= 5,
		RANK_TOLERANCE					= 6
};

// ///////////////////////////////////////////////////////////
// C++ functions for setting QPOPT options.

/// Reset all of QPOPT's options to the defaults.
void reset_defaults();

/// Set an integer valued option
void set_int_option(EQPOPT_int_option option, const f_int&);

/// Set a logical valued option
void set_logical_option(EQPOPT_logical_option option, const f_logical&);

/// Set a real valued option
void set_real_option(EQPOPT_real_option option, const f_dbl_prec&);

// ////////////////////////////////////////////////////
// Declarations to link with Fortran QPSOL procedures

extern "C" {

FORTRAN_FUNC_DECL(void,QPOPT) ( const f_int& N, const f_int& NCLIN
	, const f_int& LDA, const f_int& LDH, const f_dbl_prec* A
	, const f_dbl_prec* BL, const f_dbl_prec* BU, const f_dbl_prec* CVEC
	, const f_dbl_prec* H, qphess_func QPHESS, f_int* ISTATE, f_dbl_prec* X
	, f_int& INFORM, f_int& ITER, f_dbl_prec& OBJ, f_dbl_prec* AX
	, f_dbl_prec* CLAMDA, f_int* IW, const f_int& LENIW, f_dbl_prec* W
	, const f_int& LENW );

FORTRAN_FUNC_DECL(void,QPHESS) ( const f_int& N, const f_int& LDH
	, const f_int& JTHCOL, const f_dbl_prec* H, const f_dbl_prec* X, f_dbl_prec* HX
	, f_int* IW, const f_int& LENIW, f_dbl_prec* W, const f_int& LENW );

FORTRAN_FUNC_DECL(void,QPOPT_SET_DEFAULTS) ();

FORTRAN_FUNC_DECL(void,QPOPT_INT_OPT) (const f_int& option, const f_int& );

FORTRAN_FUNC_DECL(void,QPOPT_LOG_OPT) (const f_int& option, const f_logical& );

FORTRAN_FUNC_DECL(void,QPOPT_REAL_OPT) (const f_int& option, const f_dbl_prec& );

} // end extern "C"

// ///////////////////////////////////////////////////////////////////////////////
// Inline definitions.

inline
void qpopt( const f_int& N, const f_int& NCLIN
	, const f_int& LDA, const f_int& LDH, const f_dbl_prec* A
	, const f_dbl_prec* BL, const f_dbl_prec* BU, const f_dbl_prec* CVEC
	, const f_dbl_prec* H, qphess_func QPHESS, f_int* ISTATE, f_dbl_prec* X
	, f_int& INFORM, f_int& ITER, f_dbl_prec& OBJ, f_dbl_prec* AX
	, f_dbl_prec* CLAMDA, f_int* IW, const f_int& LENIW, f_dbl_prec* W
	, const f_int& LENW )
{
	FORTRAN_FUNC_CALL(QPOPT) ( N, NCLIN, LDA, LDH, A, BL, BU, CVEC, H, QPHESS
		, ISTATE, X, INFORM, ITER, OBJ, AX, CLAMDA, IW, LENIW, W, LENW );
}

inline
void reset_defaults()
{	FORTRAN_FUNC_CALL(QPOPT_SET_DEFAULTS) (); }

inline
void set_int_option(EQPOPT_int_option option, const f_int& val)
{	FORTRAN_FUNC_CALL(QPOPT_INT_OPT) ( option, val ); }

inline
void set_logical_option(EQPOPT_logical_option option, const f_logical& val)
{	FORTRAN_FUNC_CALL(QPOPT_LOG_OPT) ( option, val ); }

inline
void set_real_option(EQPOPT_real_option option, const f_dbl_prec& val)
{	FORTRAN_FUNC_CALL(QPOPT_REAL_OPT) ( option, val ); }

//@}

} // end namespace QPOPT_CppDecl

#endif // QPOPT_CPP_DECL_H