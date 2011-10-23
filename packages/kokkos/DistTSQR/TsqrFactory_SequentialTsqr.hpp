//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef __TSQR_Trilinos_TsqrFactory_SequentialTsqr_hpp
#define __TSQR_Trilinos_TsqrFactory_SequentialTsqr_hpp

/// \file TsqrFactory_SequentialTsqr.hpp
/// \brief Declaration and definition of SequentialTsqrFactory.
///

#include "Tsqr_SequentialTsqr.hpp"
#include "Tsqr.hpp"
#include "Teuchos_ParameterListExceptions.hpp"


namespace TSQR {
  namespace Trilinos {

    /// \class SequentialTsqrFactory
    ///
    /// Subclass of \c TsqrFactory that knows how to instantiate
    /// \c SequentialTsqr as the intranode TSQR implementation.
    ///
    /// \tparam LO "LocalOrdinal": the type of indices into the
    ///   node-local part of the matrix.
    ///
    /// \tparam S "Scalar": the type of entries in the node-local part
    ///   of the matrix.
    ///
    /// All of this class' public methods, other than the constructor
    /// and destructor, are implemented in the parent class.
    template<class LO, class S>
    class SequentialTsqrFactory :
      public TsqrFactory<LO, S, SequentialTsqr<LO, S>, DistTsqr<LO, S> > {
    public:
      // This class' parent class.
      typedef TsqrFactory<LO, S, SequentialTsqr<LO, S>, DistTsqr<LO, S> > base_type;

      // Pull in the typedefs from the base class.  C++ doesn't do
      // this when both the base and the derived classes are
      // templated.
      typedef typename base_type::node_tsqr_type node_tsqr_type;
      typedef typename base_type::dist_tsqr_type dist_tsqr_type;
      typedef typename base_type::tsqr_type tsqr_type;
      typedef typename base_type::scalar_messenger_type scalar_messenger_type;

      SequentialTsqrFactory () {}
      virtual ~SequentialTsqrFactory () {}
    };

  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrFactory_SequentialTsqr_hpp
