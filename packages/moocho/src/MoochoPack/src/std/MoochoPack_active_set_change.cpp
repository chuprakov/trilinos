// //////////////////////////////////////////////////////////////////////////
// active_set_change.cpp
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

#include <assert.h>

#include <iomanip>
#include <ostream>

#include "../std/MoochoPack_active_set_change.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_SpVectorClass.hpp"

void MoochoPack::active_set_change(
	const SpVectorSlice& nu_k, const SpVectorSlice& nu_km1, Range1D indep
	,EJournalOutputLevel olevel, std::ostream* out
	,size_type* num_adds, size_type* num_drops
	,size_type* num_active_indep, size_type* num_adds_indep, size_type* num_drops_indep
	)
{
	using std::setw;
	using std::endl;
	
	const int w = 12;

	*num_adds = *num_drops = *num_adds_indep = *num_drops_indep = 0;
	*num_active_indep = nu_k(indep).nz();

	if( !nu_k.nz() && !nu_km1.nz() )
		return;

	assert( nu_k.is_sorted() && nu_km1.is_sorted() );

	bool dump_change = (int)olevel >= (int)PRINT_ACTIVE_SET;

	if( dump_change ) {
		*out
			<< "\n*** Changes in active set\n\n"
			<< setw(w) << "i"
			<< setw(w) << "uplo"
			<< setw(w) << "dep/indep"
			<< setw(w) << "change\n"
			<< setw(w) << "----------"
			<< setw(w) << "----------"
			<< setw(w) << "----------"
			<< setw(w) << "----------\n";
	}

	SpVectorSlice::const_iterator
		nu_k_itr	= nu_k.begin(),
		nu_k_end	= nu_k.end(),
		nu_km1_itr	= nu_km1.begin(),
		nu_km1_end	= nu_km1.end();

	while( nu_k_itr != nu_k_end || nu_km1_itr != nu_km1_end ) {
		if( nu_k_itr != nu_k_end && ( nu_km1_itr == nu_km1_end || ( nu_k_itr != nu_k_end
				&& nu_k_itr->indice()+nu_k.offset() < nu_km1_itr->indice()+nu_km1.offset() ) ) )
		{
			// *nu_k_itr was added to active set.
			const size_type i = nu_k_itr->indice() + nu_k.offset();
			const bool is_indep = indep.in_range(i);
			if(is_indep)
				(*num_adds_indep)++;
			(*num_adds)++;
			if(dump_change)
				*out
					<< setw(w) << i
					<< setw(w) << ( nu_k_itr->value() >= 0.0 ? "upper" : "lower" )
					<< setw(w) << ( is_indep ? "indep" : "dep" )
					<< setw(w) << "added" << endl;
			nu_k_itr++;
		}
		else if( nu_km1_itr != nu_km1_end && ( nu_k_itr == nu_k_end || ( nu_km1_itr != nu_km1_end
				&& nu_k_itr->indice()+nu_k.offset() > nu_km1_itr->indice()+nu_km1.offset() ) ) )
		{
			// *nu_km1_itr was removed from the active set.
			const size_type i = nu_km1_itr->indice() + nu_km1.offset();
			const bool is_indep = indep.in_range(i);
			if(is_indep)
				(*num_drops_indep)++;
			(*num_drops)++;
			if(dump_change)
				*out
					<< setw(w) << i
					<< setw(w) << ( nu_km1_itr->value() >= 0.0 ? "upper" : "lower" )
					<< setw(w) << ( is_indep ? "indep" : "dep" )
					<< setw(w) << "dropped" << endl;
			nu_km1_itr++;
		}
		else {
			// same variable (but the bound may have changed)
			const size_type i = nu_k_itr->indice() + nu_k.offset();
			const bool is_indep = indep.in_range(i);
			if( nu_k_itr->value() * nu_km1_itr->value() < 0.0 ) {
				// Switched bounds.
				if(is_indep) {
					(*num_adds_indep)++;
					(*num_drops_indep)++;
				}
				(*num_adds)++;
				(*num_drops)++;
			if(dump_change)
				*out
					<< setw(w) << i
					<< setw(w) << ( nu_k_itr->value() >= 0.0 ? "upper" : "lower" )
					<< setw(w) << ( is_indep ? "indep" : "dep" )
					<< setw(w) << "switch bnd" << endl;
			}
			nu_k_itr++;
			nu_km1_itr++;
		}
	}

	// Output summary
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		*out
			<< "\n*** Active set change summary\n"
			<< "nact_old  = "			<< nu_km1.nz()			<< endl
			<< "nact_new  = "			<< nu_k.nz()			<< endl
			<< "num_adds  = "			<< *num_adds			<< endl
			<< "num_drops = "			<< *num_drops			<< endl
			<< "nact_indep_old  = "		<< nu_km1(indep).nz()	<< endl
			<< "nact_indep_new  = "		<< *num_active_indep	<< endl
			<< "num_indep_adds  = "		<< *num_adds_indep		<< endl
			<< "num_indep_drops = "		<< *num_drops_indep		<< endl;
	}
}
