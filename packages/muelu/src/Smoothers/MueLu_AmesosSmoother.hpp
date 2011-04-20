#ifndef MUELU_AMESOS_SMOOTHER_HPP
#define MUELU_AMESOS_SMOOTHER_HPP

#include "MueLu_SmootherBase.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Utilities.hpp"

#ifdef HAVE_MUELU_AMESOS
#include "Amesos_BaseSolver.h"
#include "Amesos.h"
#include "Epetra_LinearProblem.h"

namespace MueLu {

template <class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node, class LocalMatOps>
class Level;

/*!
  @class AmesosSmoother
  @brief Class that encapsulates Amesos direct solvers.

  This class creates an Amesos preconditioner factory.  The factory is capable of generating direct solvers
  based on the type and ParameterList passed into the constructor.  See the constructor for more information.
*/

  template<class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class AmesosSmoother : public SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>
  {

#include "MueLu_UseShortNames.hpp"

  private:

    //! amesos-specific key phrase that denote smoother type (same as SmootherBase::Type_)
    std::string amesosType_;
    //! pointer to Amesos solver object
    RCP<Amesos_BaseSolver> prec_;
    //! matrix operator 
    Teuchos::RCP<Operator> A_;
    //! parameter list that is used by Amesos internally
    Teuchos::ParameterList list_;
    //! Problem that Amesos uses internally.
    RCP<Epetra_LinearProblem> AmesosLinearProblem;

  protected:
    Teuchos::RCP<Teuchos::FancyOStream> out_;

  public:

    //! @name Constructors / destructors
    //@{

    /*! @brief Constructor

        Creates a MueLu interface to the direct solvers in the Amesos package.  The options are those specified in
        the Amesos user's manual.

        @param type solver type
        @param list options for the particular solver type

        Here is how to select the more commonly used direct solvers:

        - KLU (serial sparse direct solver)
            - <tt>type</tt> = <tt>Amesos-KLU</tt>
            - parameter list options
                - none required

        - SuperLU (serial sparse super-nodal direct solver)
            - <tt>type</tt> = <tt>Amesos-SuperLU</tt>
            - parameter list options
                - none required

        See also Amesos_Klu and Amesos_Superlu.

    */
    AmesosSmoother(std::string const & type, Teuchos::ParameterList const & list)
      : amesosType_(type), list_(list), out_(this->getOStream())
    {
      SmootherBase::SetType(type);
    }

    //! Destructor
    virtual ~AmesosSmoother() {}
    //@}

    //! @name Set/Get methods
    //@{

    //! @brief This has no effect and will throw an error.
    void SetNIts(LO const &nIts) {
      throw(Exceptions::RuntimeError("Only one iteration of Amesos solve is supported."));
    }

    //! @brief Returns 1.
    LO GetNIts() {
      return 1;
    }
    //@}

    //! @name Setup and Apply methods.
    //@{


    /*! @brief Set up the direct solver.

       This creates the underlying Amesos solver object according to the parameter list options passed into the
       AmesosSmoother constructor.  This includes doing a numeric factorization of the matrix.
    */
    void Setup(Level &level) {
      Teuchos::OSTab tab(out_);
      //MueLu_cout(Teuchos::VERB_HIGH) << "AmesosSmoother::Setup()" << std::endl;
      SmootherPrototype::IsSetup(true);
      A_ = level.GetA();
      RCP<Epetra_CrsMatrix> epA = Utils::Op2NonConstEpetraCrs(A_);
      AmesosLinearProblem = rcp(new Epetra_LinearProblem());
      AmesosLinearProblem->SetOperator(&*epA); //FIXME RCP probably has a safer way to do this
      Amesos factory;
      prec_ = rcp(factory.Create(amesosType_, *AmesosLinearProblem));
      if (prec_ == Teuchos::null) {
        std::string msg = "Amesos::Create: factorization type '" + amesosType_ + "' is not supported";
        throw(Exceptions::RuntimeError(msg));
      }
      prec_->SetParameters(list_);
      int rv = prec_->NumericFactorization();
      if (rv != 0) {
        std::ostringstream buf;
        buf << rv;
        std::string msg = "Amesos_BaseSolver::NumericFactorization return value of " + buf.str(); 
        throw(Exceptions::RuntimeError(msg));
      }
    }

    /*! @brief Apply the direct solver.

        Solves the linear system <tt>AX=B</tt> using the constructed solver.

        @param X initial guess
        @param B right-hand side
        @param InitialGuessIsZero This option has no effect.
    */
    void Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero=false)
    {
      if (!SmootherPrototype::IsSetup())
        throw(Exceptions::RuntimeError("Setup has not been called"));

      Epetra_MultiVector &epX = Utils::MV2NonConstEpetraMV(X);
      Epetra_MultiVector const &epB = Utils::MV2EpetraMV(B);
      //Epetra_LinearProblem takes the right-hand side as a non-const pointer.
      //I think this const_cast is safe because Amesos won't modify the rhs.
      Epetra_MultiVector &nonconstB = const_cast<Epetra_MultiVector&>(epB);
      AmesosLinearProblem->SetLHS(&epX);
      AmesosLinearProblem->SetRHS(&nonconstB);

      prec_->Solve();

      // Don't keep pointers to our vectors in the Epetra_LinearProblem.
      AmesosLinearProblem->SetLHS(0);
      AmesosLinearProblem->SetRHS(0);
    }
    //@}

    //! @name Utilities.
    //@{

    void Print(std::string prefix) {
      throw(Exceptions::NotImplemented("AmesosSmoother::Print is not implemented"));
    }

    RCP<SmootherPrototype> Copy()
    {
      return rcp(new AmesosSmoother(*this) );
    }

    void CopyParameters(RCP<SmootherPrototype> source)
    {
      RCP<AmesosSmoother> amesosSmoo = rcp_dynamic_cast<AmesosSmoother>(source);
      //TODO check if dynamic cast fails
      amesosType_ = amesosSmoo->amesosType_;
      prec_ = amesosSmoo->prec_;
      A_ = amesosSmoo->A_;
      list_ = amesosSmoo->list_;
    }
    //@}

  }; //class AmesosSmoother

} //namespace MueLu

#define MUELU_AMESOS_SMOOTHER_SHORT

#endif //ifdef HAVE_MUELU_AMESOS

#endif //ifndef MUELU_AMESOS_SMOOTHER_HPP
