// //////////////////////////////////////////////////////////////////////////
// active_set_change.h

#ifndef ACTIVE_SET_CHANGE_H
#define ACTIVE_SET_CHANGE_H

#include <iosfwd>

#include "../ReducedSpaceSQPPackTypes.h"

namespace ReducedSpaceSQPPack {

///
/** Calculate the change in the active set and output change
  * if asked to.
  *
  * ToDo: Add more description of the output you get.
  *
  * @param	nu_k	[in] Multipliers for variable bounds for iteration k.
  * @param	num_km1	[in] Multipliers for variable bounds for iteration k-1
  * @param	olevel	[in] Specifies the output level.  We have:\\
  *				PRINT_NOTHING : No output is sent to out\\
  *				PRINT_ALGORITHM_STEPS :
  *					Just the number of additions
  *					and deletions to the active set and the total number
  *					of active constraints is output.\\ 
  *				PRINT_ACTIVE_SET : Enumerates
  *					which variable were added and dropped from the active set.\\
  * @param	num_adds
  *                 [out] Gives the total number of variables fixed at a bound
  *					added to the active set.
  * @param	num_drops
  *                 [out] Gives the total number of variables freed from a 
  *					bound and dropped from the active set.
  * @param	num_adds_indep
  *                 [out] Gives the number of independent variables fixed at a bound
  *					added to the active set.
  * @param	num_drops_indep
  *                 [out] Gives the number of independent variables freed from a 
  *					bound and dropped from the active set.
  * @param	out	[O] Target for output.
  */
void active_set_change(
	const SpVectorSlice& nu_k, const SpVectorSlice& nu_km1, Range1D indep
	,EJournalOutputLevel olevel, std::ostream* out
	,size_type* num_adds, size_type* num_drops
	,size_type* num_active_indep, size_type* num_adds_indep, size_type* num_drops_indep
	); 

}	// end namespace ReducedSpaceSQPPack

#endif	// ACTIVE_SET_CHANGE_H
