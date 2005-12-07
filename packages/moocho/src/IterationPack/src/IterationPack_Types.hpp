// ///////////////////////////////////////////////////////////////////////////////////////////////////
// IterationPack_Types.hpp
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

#ifndef GENERAL_ITERATION_PACK_TYPES_H
#define GENERAL_ITERATION_PACK_TYPES_H

#include <stdexcept>

namespace IterationPack {

///
enum EAssocStepType {
	PRE_STEP   = 0
	,POST_STEP = 1
};
///
enum EDoStepType {
	DO_MAIN_STEP  = 0
	,DO_PRE_STEP  = 1
	,DO_POST_STEP = 2
};
///
enum EAlgoReturn {
	TERMINATE_TRUE
	,TERMINATE_FALSE
	,MAX_ITER_EXCEEDED
	,MAX_RUN_TIME_EXCEEDED
	,INTERRUPTED_TERMINATE_TRUE
	,INTERRUPTED_TERMINATE_FALSE
};
///
class InvalidTypeCastException : public std::logic_error
{public: InvalidTypeCastException(const std::string& what_arg) : std::logic_error(what_arg) {}};

class IterQuantity;
template<class T> class IterQuantityAccess;
template<class T> class IterQuantityAccessContiguous;
template<class T_Base, class T_Derived> class IterQuantityAccessDerivedToBase;
class Algorithm;
class AlgorithmStep;
class AlgorithmState;
class AlgorithmTracker;
class AlgorithmTrackerComposite;
template<class T> class CastIQMember;

}	// end namespace IterationPack

#endif	// GENERAL_ITERATION_PACK_TYPES_H
