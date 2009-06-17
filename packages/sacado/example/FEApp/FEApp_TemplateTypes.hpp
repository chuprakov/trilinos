// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef FEAPP_TEMPLATETYPES_HPP
#define FEAPP_TEMPLATETYPES_HPP

#include "Sacado_ConfigDefs.h"

// Conditionally enable Stokhos stochastic Galerkin expansions
#ifdef HAVE_SACADO_STOKHOS
#define SG_ACTIVE 1
#else
#define SG_ACTIVE 0
#endif

// Include the MPL vector
#include "Sacado_mpl_vector.hpp"

// Include all of our AD types
#include "Sacado_MathFunctions.hpp"
#include "Sacado.hpp"
//#include "Sacado_CacheFad_DFad.hpp"
#if SG_ACTIVE
#include "Sacado_PCE_OrthogPoly.hpp"
#endif

// Include ScalarParameterLibrary to specialize traits
#include "Sacado_ScalarParameterLibrary.hpp"

// Typedef AD types to standard names
typedef double RealType;
typedef Sacado::Fad::DFad<double> FadType;
#if SG_ACTIVE
typedef Sacado::PCE::OrthogPoly<double> SGType;
typedef Sacado::Fad::DFad< Sacado::PCE::OrthogPoly<double> > SGFadType;
#endif

// Set up evaluation traits
namespace FEApp {

  // Evaluation types
  struct ResidualType {};
  struct JacobianType {};
  struct TangentType {};
#if SG_ACTIVE
  struct SGResidualType {};
  struct SGJacobianType {};
#endif

  // Evaluation traits
  struct EvaluationTraits {
    template <class EvalT> struct apply {};
  };

  template <> struct EvaluationTraits::apply<ResidualType> { 
    typedef RealType type; };
  template <> struct EvaluationTraits::apply<JacobianType> { 
    typedef FadType type; };
  template <> struct EvaluationTraits::apply<TangentType> { 
    typedef FadType type; };
#if SG_ACTIVE
  template <> struct EvaluationTraits::apply<SGResidualType> { 
    typedef SGType type; };
  template <> struct EvaluationTraits::apply<SGJacobianType> { 
    typedef SGFadType  type; };
#endif
}

// Synonym for the ScalarParameterLibrary/Vector on our traits
typedef Sacado::ScalarParameterLibrary<FEApp::EvaluationTraits> ParamLib;
typedef Sacado::ScalarParameterVector<FEApp::EvaluationTraits> ParamVec;

// Build the MPL vector containing our valid evaluation types
typedef Sacado::mpl::vector<> EvalTypes0;
typedef Sacado::mpl::push_back<EvalTypes0, FEApp::ResidualType>::type EvalTypes1;
typedef Sacado::mpl::push_back<EvalTypes1, FEApp::JacobianType>::type EvalTypes2;
typedef Sacado::mpl::push_back<EvalTypes2, FEApp::TangentType>::type EvalTypes3;
#if SG_ACTIVE
typedef Sacado::mpl::push_back<EvalTypes3, FEApp::SGResidualType>::type EvalTypes4;
typedef Sacado::mpl::push_back<EvalTypes4, FEApp::SGJacobianType>::type EvalTypes5;
#else
typedef EvalTypes3 EvalTypes5;
#endif
typedef EvalTypes5 EvalTypes;

// Turn on/off explicit template instantiation
#define SACADO_ETI

// Define macro for explicit template instantiation
#define INSTANTIATE_TEMPLATE_CLASS_RESIDUAL(name) template class name<FEApp::ResidualType>;
#define INSTANTIATE_TEMPLATE_CLASS_JACOBIAN(name) template class name<FEApp::JacobianType>;
#define INSTANTIATE_TEMPLATE_CLASS_TANGENT(name) template class name<FEApp::TangentType>;
#if SG_ACTIVE
#define INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL(name) template class name<FEApp::SGResidualType>;
#define INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN(name) template class name<FEApp::SGJacobianType>;
#else
#define INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL(name)
#define INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN(name)
#endif

#define INSTANTIATE_TEMPLATE_CLASS(name)      \
  INSTANTIATE_TEMPLATE_CLASS_RESIDUAL(name)   \
  INSTANTIATE_TEMPLATE_CLASS_JACOBIAN(name)   \
  INSTANTIATE_TEMPLATE_CLASS_TANGENT(name)    \
  INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL(name) \
  INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN(name)

#endif // FEAPP_TEMPLATETYPES_HPP
