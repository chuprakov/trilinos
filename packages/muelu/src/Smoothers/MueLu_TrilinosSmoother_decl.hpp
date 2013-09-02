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
#ifndef MUELU_TRILINOSSMOOTHER_DECL_HPP
#define MUELU_TRILINOSSMOOTHER_DECL_HPP

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TrilinosSmoother_fwd.hpp"
#include "MueLu_SmootherPrototype.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_IfpackSmoother_fwd.hpp"

#if defined(HAVE_MUELU_IFPACK2)
#include "MueLu_Ifpack2Smoother_fwd.hpp"
#include "MueLu_Ifpack2Smoother.hpp"
#endif

// Note: TrilinosSmoother is a SmootherPrototype that cannot be turned into a smoother using Setup().
//       When this prototype is cloned using Copy(), the clone is an Ifpack or an Ifpack2 smoother.
//       The clone can be used as a smoother after calling Setup().

namespace MueLu {

  /*!
    @class TrilinosSmoother
    @brief Class that encapsulates external library smoothers.

    Autoselection of Ifpack or Ifpack2 according to the underlying linear algebra library.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType, class LocalMatOps = typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class TrilinosSmoother : public SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>
  {
#undef MUELU_TRILINOSSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors / destructors
    //@{

    /*! @brief Constructor

         @param[in] type Smoother type.  Can currently be
           - RELAXATION for (symmetric) Gauss-Seidel or Jacobi smoothing
           - CHEBYSHEV for Chebyshev polynomial smoothing
           - ILU for incomplete LU smoothing

         @param[in] paramList  A list holding parameters understood by either Ifpack or Ifpack2 to guide the smoother
         setup.  See MueLu::Ifpack2Smoother for summary of the most commonly used parameters. See the Ifpack/Ifpack2
         documentation for the full set.
    */

    template<class Scalar2, class LocalOrdinal2, class GlobalOrdinal2, class Node2, class LocalMatOps2>
    friend class TrilinosSmoother;


    TrilinosSmoother(std::string const & type = "", Teuchos::ParameterList const & paramList = Teuchos::ParameterList(), LO const &overlap=0, RCP<FactoryBase> AFact = Teuchos::null);

    //! Destructor
    virtual ~TrilinosSmoother() { }

    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    //@}

    //! @name Setup and Apply methods.
    //@{

    //! TrilinosSmoother cannot be turned into a smoother using Setup(). Setup() always returns a RuntimeError exception.
    void Setup(Level &currentLevel);

    //! TrilinosSmoother cannot be applied. Apply() always returns a RuntimeError exception.
    void Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero=false) const;

    //@}

    //! When this prototype is cloned using Copy(), the clone is an Ifpack or an Ifpack2 smoother.
    RCP<SmootherPrototype> Copy() const;

    template<typename Node2, typename LocalMatOps2>
    RCP<MueLu::TrilinosSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node2,LocalMatOps2> > clone(const RCP<Node2>& node2, const Teuchos::RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node2,LocalMatOps2> >& A_newnode) const;



    //! Convert an Ifpack2 preconditioner name to Ifpack
    // As a temporary solution.
    // See https://software.sandia.gov/bugzilla/show_bug.cgi?id=5283#c5 for what I proposed to do
    static std::string Ifpack2ToIfpack1Type(std::string const & type);

    //! Convert an Ifpack2 parameter list to Ifpack
    // As a temporary solution.
    // See https://software.sandia.gov/bugzilla/show_bug.cgi?id=5283#c5 for what I proposed to do
    /* TODO:
      ifpackList.set("type", "Chebyshev");
      ifpackList.set("chebyshev: degree", (int) 1);
      ifpackList.set("chebyshev: max eigenvalue", (double) 2.0);
      ifpackList.set("chebyshev: min eigenvalue", (double) 1.0);
      ifpackList.set("chebyshev: zero starting solution", false);
    */
    static Teuchos::ParameterList Ifpack2ToIfpack1Param(Teuchos::ParameterList const & ifpack2List);

    //! @name Overridden from Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const;

    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    //void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const {
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

    //@}

  private:

    //! ifpack1/2-specific key phrase that denote smoother type
    std::string type_;

    //! parameter list that is used by Ifpack/Ifpack2 internally
    Teuchos::ParameterList paramList_;

    //! overlap when using the smoother in additive Schwarz mode
    LO overlap_;

    //! A Factory
    RCP<FactoryBase> AFact_;

    //
    // Underlying Smoother
    //

    //! Smoother
    RCP<SmootherPrototype> s_;

  }; // class TrilinosSmoother

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  template<typename Node2, typename LocalMatOps2>
  Teuchos::RCP<MueLu::TrilinosSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node2,LocalMatOps2> >
  TrilinosSmoother<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::clone(const RCP<Node2>& node2, const Teuchos::RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node2, LocalMatOps2> >& A_newnode) const {
#if defined(HAVE_MUELU_IFPACK2)
    RCP<TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node2> > cloneSmoother = rcp(new TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node2>(type_, paramList_, overlap_, AFact_));
    Teuchos::RCP<MueLu::SmootherBase<Scalar, LocalOrdinal, GlobalOrdinal, Node2, LocalMatOps2> >  cloneSB;

    Teuchos::RCP<MueLu::Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > ifpack2Smoother = Teuchos::rcp_dynamic_cast<MueLu::Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >(this->s_);
    if (ifpack2Smoother != Teuchos::null){
      cloneSB = ifpack2Smoother->template clone<Node2, LocalMatOps2>(node2, A_newnode);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::invalid_argument, "MueLu::TrilinosSmoother: "
          "Invalid smoother type to clone (not type Ifpack2)\"");
    }
    RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node2, LocalMatOps2> > clonedProto = Teuchos::rcp_dynamic_cast<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node2, LocalMatOps2>  >(cloneSB);
    cloneSmoother->s_ = clonedProto;
    cloneSmoother->IsSetup(true);
    return cloneSmoother;
#else
  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::invalid_argument, "MueLu::SmootherClone: "
      "clone() only available with IFPACK2 enabled.");
#endif
  }


} // namespace MueLu

#define MUELU_TRILINOSSMOOTHER_SHORT
#endif // MUELU_TRILINOSSMOOTHER_DECL_HPP
