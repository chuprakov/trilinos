// ///////////////////////////////////////////////////////////////////////////////////////////////////
// GeneralIterationPackTypes.h

#ifndef GENERAL_ITERATION_PACK_TYPES_H
#define GENERAL_ITERATION_PACK_TYPES_H

#include <stdexcept>

namespace GeneralIterationPack {

///
enum EAssocStepType { PRE_STEP = 0, POST_STEP = 1 };	// Algorithm dependends on the values of these.
///
enum EDoStepType { DO_MAIN_STEP = 0, DO_PRE_STEP = 1, DO_POST_STEP = 2 };
///
enum EAlgoReturn { TERMINATE_TRUE ,TERMINATE_FALSE ,MAX_ITER_EXCEEDED
	, MAX_RUN_TIME_EXCEEDED };

///
class InvalidTypeCastException : public std::logic_error
{public: InvalidTypeCastException(const std::string& what_arg) : std::logic_error(what_arg) {}};


class IterQuantity;
template<class T> class IterQuantityAccess;
template<class T> class IterQuantityAccessContinuous;
template<class T_Base, class T_Derived> class IterQuantityAccessDerivedToBase;
class Algorithm;
class AlgorithmStep;
class AlgorithmAssocStep;
class AlgorithmState;
class AlgorithmTrack;
template<class T> class CastIQMember;

}	// end namespace GeneralIterationPack

#endif	// GENERAL_ITERATION_PACK_TYPES_H
