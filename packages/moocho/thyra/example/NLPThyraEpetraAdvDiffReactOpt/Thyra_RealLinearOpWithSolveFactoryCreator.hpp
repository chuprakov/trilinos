#include "Moocho_ConfigDefs.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

namespace Thyra {

/** \brief A utility class for creating a <tt>LinearOpWithSolveFactory</tt>
 * from the command-line and an XML file.
 */
class RealLinearOpWithSolveFactoryCreator {
public:

  /** \brief . */
  RealLinearOpWithSolveFactoryCreator();

  /** \brief . */
  void setupCLP( Teuchos::CommandLineProcessor *clp );

  /** \brief . */
  Teuchos::RefCountPtr<LinearOpWithSolveFactoryBase<double> >
  createLOWSF( std::ostream *out = NULL ) const;

  /** \brief . */
  void writeParamsUsedFile(
    const LinearOpWithSolveFactoryBase<double> &lowsFactory
    ) const;

private:

  enum ELOWSFactoryType {
    LOWSF_AMESOS
#ifdef HAVE_AZTECOO_THYRA
    ,LOWSF_AZTECOO
#endif
#ifdef HAVE_BELOS_THYRA
    ,LOWSF_BELOS
#endif
  };

  const static int numLOWSFactoryTypes_ =
    1
#ifdef HAVE_AZTECOO_THYRA
    +1
#endif
#ifdef HAVE_BELOS_THYRA
    +1
#endif
    ;
  const static ELOWSFactoryType LOWSFactoryTypeValues_[numLOWSFactoryTypes_];
  const static char* LOWSFactoryTypeNames_[numLOWSFactoryTypes_];
  
  ELOWSFactoryType    lowsFactoryType_;
#if defined(HAVE_TEUCHOS_EXTENDED) && defined(HAVE_TEUCHOS_EXPAT)
  std::string         lowsfParamsFile_;
  std::string         lowsfExtraParams_;
  std::string         lowsfParamsUsedFile_;
#endif
  bool                usePrec_;
  
};

} // namespace Thyra
