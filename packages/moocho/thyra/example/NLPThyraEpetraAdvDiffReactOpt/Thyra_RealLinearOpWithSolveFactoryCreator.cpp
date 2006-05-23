



#include "Thyra_RealLinearOpWithSolveFactoryCreator.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#ifdef HAVE_AZTECOO_THYRA
#  include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"
#endif
#ifdef HAVE_BELOS_THYRA
#  include "Thyra_BelosLinearOpWithSolveFactory.hpp"
#endif
#ifdef HAVE_IFPACK_THYRA
#  include "Thyra_IfpackPreconditionerFactory.hpp"
#endif
#if defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
#  include "Teuchos_XMLParameterListHelpers.hpp"
#endif



namespace Thyra {

const RealLinearOpWithSolveFactoryCreator::ELOWSFactoryType
RealLinearOpWithSolveFactoryCreator::LOWSFactoryTypeValues_[
  RealLinearOpWithSolveFactoryCreator::numLOWSFactoryTypes_
  ] = {
  LOWSF_AMESOS
#ifdef HAVE_AZTECOO_THYRA
  ,LOWSF_AZTECOO
#endif
#ifdef HAVE_BELOS_THYRA
  ,LOWSF_BELOS
#endif
};

const char*
RealLinearOpWithSolveFactoryCreator:: LOWSFactoryTypeNames_[
  RealLinearOpWithSolveFactoryCreator::numLOWSFactoryTypes_
  ] = {
  "amesos"
#ifdef HAVE_AZTECOO_THYRA
  ,"aztecoo"
#endif
#ifdef HAVE_BELOS_THYRA
  ,"belos"
#endif
};

RealLinearOpWithSolveFactoryCreator::RealLinearOpWithSolveFactoryCreator()
  :lowsFactoryType_(LOWSF_AMESOS)
#if defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
  ,lowsfParamsFile_("")
  ,lowsfExtraParams_("")
  ,lowsfParamsUsedFile_("")
#endif
  ,usePrec_(true)
{}

void RealLinearOpWithSolveFactoryCreator::setupCLP(
  Teuchos::CommandLineProcessor *clp
  )
{
  clp->setOption(
    "lowsf", &lowsFactoryType_
    ,numLOWSFactoryTypes_,LOWSFactoryTypeValues_,LOWSFactoryTypeNames_
    ,"The implementation for the LinearOpWithSolveFactory object used to solve the state linear systems"
    );
#if defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
  clp->setOption(
    "lowsf-params-file", &lowsfParamsFile_
    ,"LOWSF parameters XML file (must be compatible with --lowsf=???)"
    );
  clp->setOption(
    "lowsf-extra-params", &lowsfExtraParams_
    ,"Extra LOWSF parameters specified as a string in XML format (must be compatible with --lowsf=???)"
    );
  clp->setOption(
    "lowsf-params-used-file", &lowsfParamsUsedFile_
    ,"File to write the LOWSF parameters that where actually used to."
    );
#endif
  clp->setOption( "use-prec", "no-use-prec",  &usePrec_, "Flag for if preconditioning is used or not" );
}

Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> >
RealLinearOpWithSolveFactoryCreator::createLOWSF( std::ostream *out_arg ) const
{

  using Teuchos::OSTab;
  typedef double Scalar;

  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::getFancyOStream(Teuchos::rcp(out_arg,false));
    
  Teuchos::RefCountPtr<Thyra::LinearOpWithSolveFactoryBase<Scalar> > lowsFactory;
  switch(lowsFactoryType_) {
    case LOWSF_AMESOS:
      *out << "\nCreating a Thyra::AmesosLinearOpWithSolveFactory object ...\n";
      lowsFactory = Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory());
      break;
#ifdef HAVE_AZTECOO_THYRA
    case LOWSF_AZTECOO:
      *out << "\nCreating a Thyra::AztecOOLinearOpWithSolveFactory object ...\n";
      lowsFactory = Teuchos::rcp(new Thyra::AztecOOLinearOpWithSolveFactory());
      break;
#endif // HAVE_AZTECOO_THYRA
#ifdef HAVE_BELOS_THYRA
    case LOWSF_BELOS:
      *out << "\nCreating a Thyra::BelosLinearOpWithSolveFactory object ...\n";
      lowsFactory = Teuchos::rcp(new Thyra::BelosLinearOpWithSolveFactory<Scalar>());
      break;
#endif // HAVE_BELOS_THYRA
    default:
      TEST_FOR_EXCEPT(true); // should never get here!
  }
#ifdef HAVE_IFPACK_THYRA
  if( usePrec_ && lowsFactory->acceptsPreconditionerFactory() ) {
    *out << "\nCreating a Thyra::IfpackPreconditionerFactory object ...\n";
    lowsFactory->setPreconditionerFactory(Teuchos::rcp(new Thyra::IfpackPreconditionerFactory()),"");
  }
#endif // HAVE_IFPACK_THYRA
  Teuchos::RefCountPtr<Teuchos::ParameterList>
    lowsfPL = Teuchos::rcp(new Teuchos::ParameterList("LOWSF"));
  if(1) {
#if defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
    if(lowsfParamsFile_.length()) {
      Teuchos::updateParametersFromXmlFile(lowsfParamsFile_,&*lowsfPL);
      if(out.get()) {
        *out << "\nLOWSF parameters read from the file \""<<lowsfParamsFile_<<"\":\n";
        lowsfPL->print(*OSTab(out).getOStream(),0,true);
      }
    }
    if(lowsfExtraParams_.length()) {
      Teuchos::updateParametersFromXmlString(lowsfExtraParams_,&*lowsfPL);
      if(out.get()) {
        *out << "\nUpdated with extra LOWSF parameters taken from the command-line:\n";
        lowsfPL->print(*OSTab(out).getOStream(),0,true);
      }
    }
#endif // defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
    lowsFactory->setParameterList(lowsfPL);
    if(out.get()) {
      *out << "\nList of all valid LOWSF parameters:\n";
      OSTab tab(out);
      *out << lowsFactory->getValidParameters()->name() << " ->\n";
      tab.incrTab();
      lowsFactory->getValidParameters()->print(*out,0,true);
    }
  }
  
  return lowsFactory;
  
}

void RealLinearOpWithSolveFactoryCreator::writeParamsUsedFile(
  const LinearOpWithSolveFactoryBase<double> &lowsFactory
  ) const
{
  // Write the LOWSF parameters that were used:
#if defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
  if(lowsfParamsUsedFile_ != "" ) {
    Teuchos::writeParameterListToXmlFile(*lowsFactory.getParameterList(),lowsfParamsUsedFile_);
  }
#endif // defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
}


} // namespace Thyra
