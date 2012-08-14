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
/*
 * MueLu_BraessSarazinSmoother_decl.hpp
 *
 *  Created on: Apr 16, 2012
 *      Author: wiesner
 */

#ifndef MUELU_BRAESSSARAZINSMOOTHER_DECL_HPP_
#define MUELU_BRAESSSARAZINSMOOTHER_DECL_HPP_
#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_ParameterList.hpp>

//Xpetra
#include <Xpetra_MapExtractor_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>
#include <Xpetra_CrsOperator_fwd.hpp>
#include <Xpetra_Operator_fwd.hpp>

//MueLu
#include "MueLu_BraessSarazinSmoother_fwd.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_FactoryManagerBase_fwd.hpp"
#include "MueLu_SmootherBase_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

#include "MueLu_SchurComplementFactory_fwd.hpp"
#include "MueLu_DirectSolver_fwd.hpp"
#include "MueLu_SmootherFactory_fwd.hpp"
#include "MueLu_FactoryManager_fwd.hpp"

namespace MueLu {

  /*!
    @class BraessSarazinSmoother
    @brief BraessSarazin smoother for 2x2 block matrices

  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class BraessSarazinSmoother : public SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>
  {
    typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractorClass;

#undef MUELU_BRAESSSARAZINSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors / destructors
    //@{

    /*! @brief Constructor
    */
    BraessSarazinSmoother(const LocalOrdinal sweeps = 1, const Scalar omega = 1.0, RCP<const FactoryBase> AFact = Teuchos::null);

    //! Destructor
    virtual ~BraessSarazinSmoother();
    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;


    //! Set factory manager for BraessSarazin internal SchurComplement handling
    void SetFactoryManager(RCP<FactoryManager> FactManager);

    //@}

    //! @name Setup and Apply methods.
    //@{

    /*! @brief Setup routine
     * In the Setup method the whole block operator and each of the blocks are stored. Also the Schur complement
     * is calculated and a smoother is set to (partially) solve it in the Apply method. Note that the smoother set
     * here is based on the Schur complement and not in A. Without the smoo_ smoother set we cannot call
     * BraessSarazinSmoother::Apply.
     */
    void Setup(Level &currentLevel);

    /*! @brief Apply the Braess Sarazin smoother.
    Smooth the matrix from the linear system <tt>AX=B</tt> using the constructed smoother.
    @param X initial guess
    @param B right-hand side
    @param InitialGuessIsZero This option has no effect.
    */
    void Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero=false) const;
    //@}

    RCP<SmootherPrototype> Copy() const;

    //! @name Overridden from Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const;

    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

    //@}

  private:

    //! smoother type
    std::string type_;

    const LocalOrdinal                    nSweeps_;         //!< number of Braess Sarazin sweeps
    const Scalar                          omega_;           //!< damping/scaling factor



    RCP<const FactoryBase>                AFact_;           //!< A Factory
    RCP<FactoryManager>                   FactManager_;     //!< Factory manager for creating the Schur Complement

    //! block operator
    RCP<Operator>                         A_;               // < ! internal blocked operator "A" generated by AFact_

    RCP<const MapExtractorClass>          rangeMapExtractor_;  //!< range  map extractor (from A_ generated by AFact)
    RCP<const MapExtractorClass>          domainMapExtractor_; //!< domain map extractor (from A_ generated by AFact)

    //Not used yet
    std::map<size_t,size_t>               bgsOrderingIndex2blockRowIndex_; //!< map: block GaussSeidel ordering 2 block row ordering

    //! matrices
    Teuchos::RCP<Vector>                  diagFinv_;                      //!< inverse diagonal of fluid operator (vector). A vector is used to make use of Multivector::elementWiseMultiply
    Teuchos::RCP<Operator>                F_;                             //!< fluid operator
    Teuchos::RCP<Operator>                G_;                             //!< pressure gradient operator
    Teuchos::RCP<Operator>                D_;                             //!< divergence operator
    Teuchos::RCP<Operator>                Z_;                             //!< pressure stabilization terms or null block

    Teuchos::RCP<SmootherBase>                  smoo_;                                  //!< Smoother to solve the schur complement


  }; // class Amesos2Smoother

} // namespace MueLu

#define MUELU_BRAESSSARAZINSMOOTHER_SHORT
#endif /* MUELU_BRAESSSARAZINSMOOTHER_DECL_HPP_ */
