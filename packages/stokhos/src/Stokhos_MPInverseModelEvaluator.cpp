// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Stokhos_MPInverseModelEvaluator.hpp"

#include "Teuchos_TestForException.hpp"
#include "Stokhos_ProductEpetraVector.hpp"
#include "Stokhos_ProductEpetraMultiVector.hpp"
#include "Epetra_Map.h"

Stokhos::MPInverseModelEvaluator::MPInverseModelEvaluator(
  const Teuchos::RCP<EpetraExt::ModelEvaluator>& me_,
  const Teuchos::Array<int>& mp_p_index_,
  const Teuchos::Array<int>& non_mp_p_index_,
  const Teuchos::Array<int>& mp_g_index_,
  const Teuchos::Array<int>& non_mp_g_index_,
  const Teuchos::Array< Teuchos::RCP<const Epetra_Map> >& base_p_maps_,
  const Teuchos::Array< Teuchos::RCP<const Epetra_Map> >& base_g_maps_) 
  : me(me_),
    mp_p_index(mp_p_index_),
    non_mp_p_index(non_mp_p_index_),
    mp_g_index(mp_g_index_),
    non_mp_g_index(non_mp_g_index_),
    base_p_maps(base_p_maps_),
    base_g_maps(base_g_maps_),
    num_p_mp(mp_p_index.size()),
    num_p(non_mp_p_index.size()),
    num_g_mp(mp_g_index.size()),
    num_g(non_mp_g_index.size()),
    block_p(num_p_mp),
    block_g(num_g_mp),
    block_dgdp(num_g_mp)
{
  TEST_FOR_EXCEPTION(
    base_p_maps.size() != num_p_mp, std::logic_error,
    std::endl 
    << "Error!  Stokhos::MPInverseModelEvaluator::MPInverseModelEvaluator():"
    << "  Base parameter map array size does not match size of index array!");
  TEST_FOR_EXCEPTION(
    base_g_maps.size() != num_g_mp, std::logic_error,
    std::endl 
    << "Error!  Stokhos::MPInverseModelEvaluator::MPInverseModelEvaluator():"
    << "  Base response map array size does not match size of index array!");

  InArgs me_inargs = me->createInArgs();
  OutArgs me_outargs = me->createOutArgs();
  
  // Create parameter MP blocks
  for (int i=0; i<num_p_mp; i++) {
    Teuchos::RCP<const Epetra_Map> mp_p_map = me->get_p_map(mp_p_index[i]);
    block_p[i] = Teuchos::rcp(new Epetra_Vector(*mp_p_map));
  }

  // Create response MP blocks
  for (int i=0; i<num_g_mp; i++) {
    Teuchos::RCP<const Epetra_Map> mp_g_map = me->get_g_map(mp_g_index[i]);
    
    // Create g MP blocks
    block_g[i] = Teuchos::rcp(new Epetra_Vector(*mp_g_map));
    
    // Create dg/dp MP blocks
    block_dgdp[i].resize(num_p);
    for (int j=0; j<num_p; j++) {
      Teuchos::RCP<const Epetra_Map> p_map = me->get_p_map(non_mp_p_index[j]);
      block_dgdp[i][j] = 
	Teuchos::rcp(new Epetra_MultiVector(*mp_g_map, p_map->NumMyElements()));
    }
  }
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
Stokhos::MPInverseModelEvaluator::
get_x_map() const
{
  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::MPInverseModelEvaluator::
get_f_map() const
{
  return Teuchos::null;
}

Teuchos::RCP<const Epetra_Map>
Stokhos::MPInverseModelEvaluator::
get_p_map(int l) const
{
  TEST_FOR_EXCEPTION(
    l >= num_p || l < 0, std::logic_error,
    std::endl << "Error!  Stokhos::MPInverseModelEvaluator::get_p_map():"
    << "  Invalid parameter index l = " << l << std::endl);
  return me->get_p_map(non_mp_p_index[l]);
}

Teuchos::RCP<const Epetra_Map>
Stokhos::MPInverseModelEvaluator::
get_p_mp_map(int l) const
{
  TEST_FOR_EXCEPTION(
    l >= num_p_mp || l < 0, std::logic_error,
    std::endl << "Error!  Stokhos::MPInverseModelEvaluator::get_p_mp_map():"
    << "  Invalid parameter index l = " << l << std::endl);
  return base_p_maps[l];
}

Teuchos::RCP<const Epetra_Map>
Stokhos::MPInverseModelEvaluator::
get_g_map(int l) const
{
  TEST_FOR_EXCEPTION(
    l >= num_g || l < 0, std::logic_error,
    std::endl << "Error!  Stokhos::MPInverseModelEvaluator::get_g_map():"
    << "  Invalid response index l = " << l << std::endl);
  return me->get_g_map(non_mp_g_index[l]);
}

Teuchos::RCP<const Epetra_Map>
Stokhos::MPInverseModelEvaluator::
get_g_mp_map(int l) const
{
  TEST_FOR_EXCEPTION(
    l >= num_g_mp || l < 0, std::logic_error,
    std::endl << "Error!  Stokhos::MPInverseModelEvaluator::get_g_mp_map():"
    << "  Invalid response index l = " << l << std::endl);
  return base_g_maps[l];
}

Teuchos::RCP<const Teuchos::Array<std::string> >
Stokhos::MPInverseModelEvaluator::
get_p_names(int l) const
{
  TEST_FOR_EXCEPTION(
    l >= num_p || l < 0, std::logic_error,
    std::endl << "Error!  Stokhos::MPInverseModelEvaluator::get_p_names():"
    << "  Invalid parameter index l = " << l << std::endl);
  return me->get_p_names(non_mp_p_index[l]);
}

Teuchos::RCP<const Teuchos::Array<std::string> >
Stokhos::MPInverseModelEvaluator::
get_p_mp_names(int l) const
{
  TEST_FOR_EXCEPTION(
    l >= num_p_mp || l < 0, std::logic_error,
    std::endl << "Error!  Stokhos::MPInverseModelEvaluator::get_p_mp_names():"
    << "  Invalid parameter index l = " << l << std::endl);
  return me->get_p_names(mp_p_index[l]);
}

Teuchos::RCP<const Epetra_Vector> 
Stokhos::MPInverseModelEvaluator::
get_p_init(int l) const
{
  TEST_FOR_EXCEPTION(
    l >= num_p || l < 0, std::logic_error,
    std::endl << "Error!  Stokhos::MPInverseModelEvaluator::get_p_init():"
    << "  Invalid parameter index l = " << l << std::endl);
  return me->get_p_init(non_mp_p_index[l]);
}

EpetraExt::ModelEvaluator::InArgs
Stokhos::MPInverseModelEvaluator::createInArgs() const
{
  InArgsSetup inArgs;
  InArgs me_inargs = me->createInArgs();

  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);
  inArgs.set_Np_mp(num_p_mp); 
  
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
Stokhos::MPInverseModelEvaluator::createOutArgs() const
{
  OutArgsSetup outArgs;
  OutArgs me_outargs = me->createOutArgs();

  outArgs.setModelEvalDescription(this->description());

  outArgs.set_Np_Ng(num_p, num_g);
  for (int i=0; i<num_g; i++)
    for (int j=0; j<num_p; j++)
      outArgs.setSupports(OUT_ARG_DgDp, i, j, 
			  me_outargs.supports(OUT_ARG_DgDp, 
					      non_mp_g_index[i], 
					      non_mp_p_index[j]));

  outArgs.set_Np_Ng_mp(num_p, num_g_mp);
  for (int i=0; i<num_g_mp; i++)
    for (int j=0; j<num_p; j++)
      outArgs.setSupports(OUT_ARG_DgDp_mp, i, j, 
			  me_outargs.supports(OUT_ARG_DgDp, 
					      mp_g_index[i], 
					      non_mp_p_index[j]));
  
  return outArgs;
}

void 
Stokhos::MPInverseModelEvaluator::evalModel(const InArgs& inArgs,
					    const OutArgs& outArgs) const
{
  // Create underlying inargs
  InArgs me_inargs = me->createInArgs();

  // Pass parameters
  for (int i=0; i<num_p; i++)
    me_inargs.set_p(non_mp_p_index[i], inArgs.get_p(i));

  // Pass MP parameters
  for (int i=0; i<num_p_mp; i++) {
    mp_const_vector_t p_mp = inArgs.get_p_mp(i);
    if (p_mp != Teuchos::null) {
      p_mp->assignToBlockVector(*(block_p[i]));
      me_inargs.set_p(mp_p_index[i], block_p[i]);
    }
  }

  // Create underlying outargs
  OutArgs me_outargs = me->createOutArgs();

  // Responses
  for (int i=0; i<num_g; i++) {
    // g
    me_outargs.set_g(non_mp_g_index[i], outArgs.get_g(i));

    // dg/dp
    for (int j=0; j<num_p; j++)
      if (!outArgs.supports(OUT_ARG_DgDp, i, j).none())
	me_outargs.set_DgDp(non_mp_g_index[i], non_mp_p_index[j], 
			    outArgs.get_DgDp(i,j));
  }


  // MP Responses
  for (int i=0; i<num_g_mp; i++) {
    // g
    mp_vector_t g_mp = outArgs.get_g_mp(i);
    if (g_mp != Teuchos::null) {
      g_mp->assignToBlockVector(*(block_g[i]));
      me_outargs.set_g(mp_g_index[i], block_g[i]);
    }

    // dg/dp
    for (int j=0; j<num_p; j++) {
      if (!outArgs.supports(OUT_ARG_DgDp_mp, i, j).none()) {
	MPDerivative dgdp_mp = outArgs.get_DgDp_mp(i,j);
	Teuchos::RCP<Stokhos::ProductEpetraMultiVector> dgdp_mp_mv
	  = dgdp_mp.getMultiVector();
	if (dgdp_mp_mv != Teuchos::null) {
	  dgdp_mp_mv->assignToBlockMultiVector(*(block_dgdp[i][j]));
	  me_outargs.set_DgDp(mp_g_index[i], non_mp_p_index[j], 
			      Derivative(block_dgdp[i][j],
					 dgdp_mp.getMultiVectorOrientation()));
	}
	TEST_FOR_EXCEPTION(dgdp_mp.getLinearOp() != Teuchos::null, 
			   std::logic_error,
			   "Error!  Stokhos::MPInverseModelEvaluator::evalModel"
			   << " cannot handle operator form of dg/dp!");
      }
    }

  }

  // Compute the functions
  me->evalModel(me_inargs, me_outargs);

  // Copy block MP components
  for (int i=0; i<num_g_mp; i++) {
    // g
    mp_vector_t g_mp = outArgs.get_g_mp(i);
    if (g_mp != Teuchos::null) {
      g_mp->assignFromBlockVector(*(block_g[i]));
    }

    // dg/dp
    for (int j=0; j<num_p; j++) {
      if (!outArgs.supports(OUT_ARG_DgDp_mp, i, j).none()) {
	MPDerivative dgdp_mp = outArgs.get_DgDp_mp(i,j);
	Teuchos::RCP<Stokhos::ProductEpetraMultiVector> dgdp_mp_mv
	  = dgdp_mp.getMultiVector();
	if (dgdp_mp_mv != Teuchos::null)
	  dgdp_mp_mv->assignFromBlockMultiVector(*(block_dgdp[i][j]));
      }
    }
  }

}
