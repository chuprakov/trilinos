#ifndef MUELU_AMESOS_SMOOTHER_HPP
#define MUELU_AMESOS_SMOOTHER_HPP

#ifdef HAVE_MUELU_AMESOS

#include <Amesos_config.h>
#include <Amesos.h>
#include <Amesos_BaseSolver.h>
#include <Epetra_LinearProblem.h>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SmootherBase.hpp"
#include "MueLu_SmootherPrototype.hpp"

namespace MueLu {

  /*!
    @class AmesosSmoother
    @brief Class that encapsulates Amesos direct solvers.
    
    This class creates an Amesos preconditioner factory.  The factory is capable of generating direct solvers
    based on the type and ParameterList passed into the constructor.  See the constructor for more information.
  */

  class AmesosSmoother : public SmootherPrototype<double, int, int>
  {

    typedef double Scalar;
    typedef int    LocalOrdinal;
    typedef int    GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;
    typedef Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps LocalMatOps;
#include "MueLu_UseShortNames.hpp"

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

        If you are using type=="", then either SuperLU or KLU are used by default.

        See also Amesos_Klu and Amesos_Superlu.

    */

    AmesosSmoother(std::string const & type = "", Teuchos::ParameterList const & paramList = Teuchos::ParameterList())
      : type_(type), paramList_(paramList)
    {

#if defined(HAVE_AMESOS_SUPERLU)
      type_ = "Superlu";
#elif defined(HAVE_AMESOS_KLU)
      type_ = "Klu";
#endif
      TEST_FOR_EXCEPTION(type_ == "", Exceptions::RuntimeError, "MueLu::AmesosSmoother::AmesosSmoother(): Amesos compiled without KLU and SuperLU. Cannot define a solver by default for this AmesosSmoother object");

      TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == true, Exceptions::RuntimeError, "TO BE REMOVED");
    }

    //! Destructor
    virtual ~AmesosSmoother() {}

    //@}

    //! @name Setup and Apply methods.
    //@{
    
    /*! @brief Set up the direct solver.
      This creates the underlying Amesos solver object according to the parameter list options passed into the
      AmesosSmoother constructor.  This includes doing a numeric factorization of the matrix.
    */
    void Setup(Level &level) {
      TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == true, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Setup(): Setup() has already been called"); //TODO: Valid. To be replace by a warning.

      A_ = level.Get< RCP<Operator> >("A",NULL); //TODO NULL!!

      RCP<Epetra_CrsMatrix> epA = Utils::Op2NonConstEpetraCrs(A_);
      linearProblem_ = rcp( new Epetra_LinearProblem() );
      linearProblem_->SetOperator(epA.get());

      Amesos factory;
      prec_ = rcp(factory.Create(type_, *linearProblem_));
      TEST_FOR_EXCEPTION(prec_ == Teuchos::null, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Setup(): Solver '" + type_ + "' not supported by Amesos");

      prec_->SetParameters(paramList_);

      int r = prec_->NumericFactorization();
      TEST_FOR_EXCEPTION(r != 0, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Setup(): Amesos solver returns value of " + Teuchos::Utils::toString(r) + " during NumericFactorization()");

      SmootherPrototype::IsSetup(true);
    }

    /*! @brief Apply the direct solver.

        Solves the linear system <tt>AX=B</tt> using the constructed solver.

        @param X initial guess
        @param B right-hand side
        @param InitialGuessIsZero This option has no effect with this smoother
    */
    void Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero = false) const
    {
      TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Apply(): Setup() has not been called");

      Epetra_MultiVector &epX = Utils::MV2NonConstEpetraMV(X);
      Epetra_MultiVector const &epB = Utils::MV2EpetraMV(B);
      //Epetra_LinearProblem takes the right-hand side as a non-const pointer.
      //I think this const_cast is safe because Amesos won't modify the rhs.
      Epetra_MultiVector &nonconstB = const_cast<Epetra_MultiVector&>(epB);

      linearProblem_->SetLHS(&epX);
      linearProblem_->SetRHS(&nonconstB);

      prec_->Solve();

      // Don't keep pointers to our vectors in the Epetra_LinearProblem.
      linearProblem_->SetLHS(0);
      linearProblem_->SetRHS(0);
    }

    //@}

    RCP<SmootherPrototype> Copy() const {
      return rcp( new AmesosSmoother(*this) );
    }
    
    //! @name Overridden from Teuchos::Describable 
    //@{
    
    //! Return a simple one-line description of this object.
    std::string description() const {
      std::ostringstream out;
      out << SmootherPrototype::description();
      out << "{type = " << type_ << "} ";
      return out.str();
    }
    
    //! Print the object with some verbosity level to an FancyOStream object.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const {
      using std::endl;
      int vl = (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;
      if (vl == VERB_NONE) return;
      
      if (vl == VERB_LOW) { out << description() << endl; } else { out << SmootherPrototype::description() << endl; }
      
      Teuchos::OSTab tab1(out);

      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        out << "Prec. type: " << type_ << endl;
        out << "Parameter list: " << endl; { Teuchos::OSTab tab2(out); out << paramList_; }
      }
      
      if (vl == VERB_HIGH || vl == VERB_EXTREME) {
        if (prec_ != Teuchos::null) { prec_->PrintStatus(); prec_->PrintTiming(); }
      }

      if (vl == VERB_EXTREME) {
        out << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << endl;
        out << "-" << endl;
        out << "RCP<A_>: " << A_ << std::endl;
        out << "RCP<linearProblem__>: " << linearProblem_ << std::endl;
        out << "RCP<prec_>: " << prec_ << std::endl; 
      }
    }

    //@}

  private:

    //! amesos-specific key phrase that denote smoother type
    std::string type_;

    //! parameter list that is used by Amesos internally
    Teuchos::ParameterList paramList_;

    //! pointer to Amesos solver object
    RCP<Amesos_BaseSolver> prec_;

    //! Problem that Amesos uses internally.
    RCP<Epetra_LinearProblem> linearProblem_;

    //! Operator. Not used directly, but held inside of linearProblem_. So we have to keep an RCP pointer to it!
    RCP<Operator> A_;

  }; // class AmesosSmoother

  //! Non-member templated function GetAmesosSmoother() returns a new AmesosSmoother object when <Scalar, LocalOrdinal, GlobalOrdinal> == <double, int, int>. Otherwise, an exception is thrown.
  //! This function simplifies the usage of AmesosSmoother objects inside of templates as templates do not have to be specialized for <double, int, int> (see DirectSolver for an example).
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > GetAmesosSmoother(std::string const & type = "", Teuchos::ParameterList const & paramList = Teuchos::ParameterList()) { 
    TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "AmesosSmoother cannot be used with Scalar != double, LocalOrdinal != int, GlobalOrdinal != int");
    return Teuchos::null;
  }

} // namespace MueLu

#define MUELU_AMESOS_SMOOTHER_SHORT

#endif // HAVE_MUELU_AMESOS
#endif // MUELU_AMESOS_SMOOTHER_HPP
