#ifndef MUELU_SMOOTHERFACTORY_HPP
#define MUELU_SMOOTHERFACTORY_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SmootherFactoryBase.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

  /*!
    @class SmootherFactory
    @brief Generic Smoother Factory for generating the smoothers of the MG hierarchy

    This factory is generic and can produce any kind of Smoother.
 
    The design of the Smoother factory is based on the prototype design
    pattern.  The smoother factory uses prototypical instances of
    smoothers, and prototypes are cloned to produce the smoothers of the MG
    hierarchy.
 
    See also:
    - http://en.wikipedia.org/wiki/Factory_method_pattern
    - http://en.wikipedia.org/wiki/Prototype_pattern
    - http://www.oodesign.com/prototype-pattern.html
 
    There is one prototype for the pre-smoother and one for the
    post-smoother. Thus, this factory can produce two different smoothers
    for pre and post-smoothing. Prototypes are stored in this factory.
 
    The type of smoother to create is determined by the prototypical
    instances. Prototypes store their own parameters (nits, omega), and the
    factory doesn't have any knowledge about that.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class SmootherFactory : public SmootherFactoryBase {

#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    /*!
      @brief Constructor

      TODO update the documentation from Matlab syntax

      The constructor takes as arguments pre-configured smoother
      object(s). They are used as prototypes by the factory.
      
      SYNTAX

      SmootherFactory(preSmootherPrototype_, postSmootherPrototype_);
      
      preSmootherPrototype_  - prototype for pre-smoothers  (SmootherPrototype)
      postSmootherPrototype_ - prototype for post-smoothers
      (SmootherPrototype,optional,default=preSmootherPrototype_)

      EXAMPLES:

      SmootherFactory = SmootherFactory(ChebySmoother(5,1/30))

      or 

      nIts        = 5;
      lambdaRatio = 1/30;
      Smoother        = ChebySmoother()
      Smoother        = Smoother.SetIts(nIts);
      Smoother        = Smoother.SetLambdaRatio(1/30); 
      SmootherFactory = SmootherFactory(Smoother);

      To use different smoothers for pre and post smoothing, two prototypes can be passed in as argument:
      PreSmoother     = ChebySmoother(2, 1/30)
      PostSmoother    = ChebySmoother(10,1/30)
      SmootherFactory = SmootherFactory(PreSmoother, PostSmoother);

      [] is also a valid smoother which do nothing:
      PostSmoother    = ChebySmoother(10,1/30)
      SmootherFactory = SmootherFactory([], PostSmoother);
    */
    //Note: Teuchos::null as parameter allowed (= no smoother) but must be explicitly defined (no parameter default value)
    //Note: precondition: input smoother must not be Setup()
    SmootherFactory(RCP<SmootherPrototype> preAndPostSmootherPrototype)
      : preSmootherPrototype_(preAndPostSmootherPrototype), postSmootherPrototype_(preAndPostSmootherPrototype)
    { 
      TEST_FOR_EXCEPTION(preAndPostSmootherPrototype != Teuchos::null && preAndPostSmootherPrototype->IsSetup() == true, Exceptions::RuntimeError, "preAndPostSmootherPrototype is not a smoother prototype (IsSetup() == true)");
    }

    SmootherFactory(RCP<SmootherPrototype> preSmootherPrototype, RCP<SmootherPrototype> postSmootherPrototype)
      : preSmootherPrototype_(preSmootherPrototype), postSmootherPrototype_(postSmootherPrototype)
    { 
      TEST_FOR_EXCEPTION(preSmootherPrototype  != Teuchos::null && preSmootherPrototype->IsSetup()  == true, Exceptions::RuntimeError, "preSmootherPrototype is not a smoother prototype (IsSetup() == true)");
      TEST_FOR_EXCEPTION(postSmootherPrototype != Teuchos::null && postSmootherPrototype->IsSetup() == true, Exceptions::RuntimeError, "postSmootherPrototype is not a smoother prototype (IsSetup() == true)");
    }

    virtual ~SmootherFactory() {}
    //@}

    //! @name Set/Get methods.
    //@{

    //! Set smoother prototypes.
    void SetSmootherPrototypes(RCP<SmootherPrototype> preAndPostSmootherPrototype) {
      TEST_FOR_EXCEPTION(preAndPostSmootherPrototype != Teuchos::null && preAndPostSmootherPrototype->IsSetup() == true, Exceptions::RuntimeError, "preAndPostSmootherPrototype is not a smoother prototype (IsSetup() == true)");

      preSmootherPrototype_ = preAndPostSmootherPrototype;
      postSmootherPrototype_ = preAndPostSmootherPrototype;
    }

    //! Set smoother prototypes.
    void SetSmootherPrototypes(RCP<SmootherPrototype> preSmootherPrototype, RCP<SmootherPrototype> postSmootherPrototype) {
      TEST_FOR_EXCEPTION(preSmootherPrototype  != Teuchos::null && preSmootherPrototype->IsSetup()  == true, Exceptions::RuntimeError, "preSmootherPrototype is not a smoother prototype (IsSetup() == true)");
      TEST_FOR_EXCEPTION(postSmootherPrototype != Teuchos::null && postSmootherPrototype->IsSetup() == true, Exceptions::RuntimeError, "postSmootherPrototype is not a smoother prototype (IsSetup() == true)");
      preSmootherPrototype_ = preSmootherPrototype;
      postSmootherPrototype_ = postSmootherPrototype;
    }
    
    //! Get smoother prototypes.
    void GetSmootherPrototypes(RCP<SmootherPrototype> & preSmootherPrototype, RCP<SmootherPrototype> & postSmootherPrototype) const {
      preSmootherPrototype = preSmootherPrototype_;
      postSmootherPrototype = postSmootherPrototype_;
    }

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const { }

    //@}

    //! @name Build methods.
    //@{

    /*! @brief Creates pre and post smoothers.

    Factory.Build() clones its prototypes and calls Setup() on the
    new objects to create fully functional smoothers for the current
    hierarchy level. Note that smoothers do their own setup
    (ie: ILUSmoother.Setup() computes the ILU factorization). 
       
    If pre and post smoother are identical, the Setup() phase
    is done only once: to create the post-smoother, the
    pre-smoother is duplicated and its parameters are changed
    according to the parameters of the post-smoother prototype.
       
    If the parameters of pre and post smoothers are not the same,
    the Setup() phase is also done only once when parameters
    don't change the result of the setup computation.
    */
    bool Build(Level& currentLevel) const {
      return BuildSmoother(currentLevel, BOTH);
    }
    
    bool BuildSmoother(Level & currentLevel, PreOrPost const preOrPost = BOTH) const {
      RCP<SmootherPrototype> preSmoother;
      RCP<SmootherPrototype> postSmoother;
      
      if ((preOrPost == BOTH || preOrPost == PRE) && (preSmootherPrototype_ != Teuchos::null)) {
        preSmoother = preSmootherPrototype_->Copy();
        //preSmoother = rcp( new SmootherPrototype(preSmootherPrototype_) );
        //TODO if outputlevel high enough
        //TODO preSmoother.Print();
        preSmoother->Setup(currentLevel);
        
        // Level Set
        currentLevel.Set<RCP<SmootherBase> >("PreSmoother", preSmoother, this);
        currentLevel.Set<RCP<SmootherBase> >("PreSmoother", preSmoother);
      }
      
      if ((preOrPost == BOTH || preOrPost == POST) && (postSmootherPrototype_ != Teuchos::null))
          {
            if (preOrPost == BOTH && preSmootherPrototype_ == postSmootherPrototype_) {
              
              // Very simple reuse. TODO: should be done in MueMat too
              postSmoother = preSmoother;
              
//            }  else if (preOrPost == BOTH &&
//                        preSmootherPrototype_ != Teuchos::null &&
//                        preSmootherPrototype_->GetType() == postSmootherPrototype_->GetType()) {
              
//               // More complex reuse case: need implementation of CopyParameters() and a smoothers smart enough to know when parameters affect the setup phase.
              
//               // YES: post-smoother == pre-smoother 
//               // => copy the pre-smoother to avoid the setup phase of the post-smoother.
//               postSmoother = preSmoother->Copy();
//               // If the post-smoother parameters are different from
//               // pre-smoother, the parameters stored in the post-smoother
//               // prototype are copied in the new post-smoother object.
//               postSmoother->CopyParameters(postSmootherPrototype_);
//               // If parameters don't influence the Setup phase (it is the case
//               // for Jacobi, Chebyshev...), PostSmoother is already setup. Nothing
//               // more to do. In the case of ILU, parameters of the smoother
//               // are in fact the parameters of the Setup phase. The call to
//               // CopyParameters resets the smoother (only if parameters are
//               // different) and we must call Setup() again.
//               postSmoother->Setup(currentLevel);

//               // TODO: if CopyParameters do not exist, do setup twice.

            } else {
              
              // NO reuse: preOrPost==POST or post-smoother != pre-smoother
              // Copy the prototype and run the setup phase.
              postSmoother = postSmootherPrototype_->Copy();
              postSmoother->Setup(currentLevel);
              
            }
            
            // Level Set
            currentLevel.Set<RCP<SmootherBase> >("PostSmoother", postSmoother, this);
            currentLevel.Set<RCP<SmootherBase> >("PostSmoother", postSmoother);
          }
        
        return true; //?
        
    } // Build()
    
    //@}


    //! @name Overridden from Teuchos::Describable 
    //@{
    
    //! Return a simple one-line description of this object.
    std::string description() const {
      std::ostringstream out;
      out << SmootherFactoryBase::description();
      std::string preStr  = (preSmootherPrototype_ == Teuchos::null) ? "null" : preSmootherPrototype_->description();
      std::string postStr = (preSmootherPrototype_ == postSmootherPrototype_) ? "pre" : ( (postSmootherPrototype_ == Teuchos::null) ? "null" : preSmootherPrototype_->description() );
      out << "{pre = "  << preStr << ", post = "<< postStr << "} ";
      return out.str();
    }
    
    //! Print the object with some verbosity level to an FancyOStream object.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const { //TODO: remove Teuchos::Describable::
      using std::endl;
      int vl = (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;
      if (vl == VERB_NONE) return;
      
      if (vl == VERB_LOW) { out << description() << endl; } else { out << SmootherFactoryBase::description() << endl; }

      Teuchos::OSTab tab1(out);
        
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        out << "PreSmoother: "; if (preSmootherPrototype_  == Teuchos::null) { out << "null" << endl; } else { Teuchos::OSTab tab2(out); preSmootherPrototype_->describe(out, verbLevel); }

        out << "PostSmoother: ";
        if      (postSmootherPrototype_ == preSmootherPrototype_) { out << "same as PreSmoother" << endl; } 
        else if (postSmootherPrototype_ == Teuchos::null)         { out << "null" << endl; }
        else { 
          { Teuchos::OSTab tab2(out); postSmootherPrototype_->describe(out, verbLevel); }
          out << "PostSmoother is different than PreSmoother (not the same object)" << endl;
        }
      }

      if (vl == VERB_EXTREME) {
        if (preSmootherPrototype_  != Teuchos::null || postSmootherPrototype_ != Teuchos::null) { out << "-" << endl; }
        if (preSmootherPrototype_  != Teuchos::null) { out << "RCP<preSmootherPrototype_> : " << preSmootherPrototype_  << std::endl; }
        if (postSmootherPrototype_ != Teuchos::null) { out << "RCP<postSmootherPrototype_>: " << postSmootherPrototype_ << std::endl; }
      }
    }

    //@}

  private:
    RCP<SmootherPrototype> preSmootherPrototype_;
    RCP<SmootherPrototype> postSmootherPrototype_;

    //@}

  }; // class SmootherFactory

} // namespace MueLu

#define MUELU_SMOOTHERFACTORY_SHORT

#endif // ifndef MUELU_SMOOTHERFACTORY_HPP

//TODO: doc: setup done twice if PostSmoother object != PreSmoother object and no adv. reused capability
