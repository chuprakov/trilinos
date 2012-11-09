// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_PATTERNFACTORY_DECL_HPP
#define MUELU_PATTERNFACTORY_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Level_fwd.hpp"

namespace MueLu {

  /*!
    @class PatternFactory class.
    @brief Factory for building nonzero patterns for energy minimization.
    @ingroup MueLuTransferClasses
    */

  template<class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::sparseOps>
  class PatternFactory : public SingleLevelFactoryBase {
#undef MUELU_PATTERNFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors.
    //@{

    /*! @brief Constructor.
      User can supply a factory for generating the prolongator. The nonzero pattern of the
      constructed prolongator will be used for the nonzero pattern.
      */
    PatternFactory(RCP<const FactoryBase> PFact = Teuchos::null);

    //! Destructor.
    virtual ~PatternFactory();

    //@}

    //! @name Input
    //@{

    void DeclareInput(Level& currentLevel) const;

    //@}

    //! @name Build methods.
    //@{

    /*!
      @brief Build method.

      Builds nonzero pattern (graph) and returns it in <tt>coarseLevel</tt>.
      */
    void Build(Level & currentLevel) const;

    //@}

  private:
    RCP<const FactoryBase> PFact_;              //!< Prolongator factory
  }; // class PatternFactory


} // namespace MueLu

#define MUELU_PATTERNFACTORY_SHORT
#endif // MUELU_PATTERNFACTORY_DECL_HPP
