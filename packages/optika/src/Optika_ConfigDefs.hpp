#ifndef OPTIKA_CONFIG_DEFS_HPP
#define OPTIKA_CONFIG_DEFS_HPP

/*! \file Optika_ConfigDefs.hpp
 * \brief A Header file that includes some of the commonly used includes
 * throughtout Optika as well as some overall information used throughout
 * Optika.
 */
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_StandardDependencies.hpp"
#include "Teuchos_DependencySheet.hpp"
#include "Teuchos_TwoDArray.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"


namespace Optika{

  //Common Teuchos classes that are used.
  using Teuchos::ParameterList;
  using Teuchos::ParameterEntry;
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::null;
  using Teuchos::is_null;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_static_cast;
  using Teuchos::ParameterEntryValidator;
  using Teuchos::EnhancedNumberValidator;
  using Teuchos::EnhancedNumberTraits;
  using Teuchos::FileNameValidator;
  using Teuchos::ArrayValidator;
  using Teuchos::DependencySheet;
  using Teuchos::Dependency;
  using Teuchos::VisualDependency;
  using Teuchos::any;
  using Teuchos::any_cast;
  using Teuchos::XMLParameterListWriter;
  using Teuchos::XMLObject;
  using Teuchos::getValue;
  using Teuchos::getParametersFromXmlFile;
  using Teuchos::TwoDArray;


  
} //namespace Optika

//Declarations of the supported Array types
//as QMEATTYPES so that we can encapsulate them
//in QVariant objects.
Q_DECLARE_METATYPE(Teuchos::Array<int>)
Q_DECLARE_METATYPE(Teuchos::Array<short>)
Q_DECLARE_METATYPE(Teuchos::Array<float>)
Q_DECLARE_METATYPE(Teuchos::Array<double>)
Q_DECLARE_METATYPE(Teuchos::Array<std::string>)
Q_DECLARE_METATYPE(Teuchos::TwoDArray<int>)
Q_DECLARE_METATYPE(Teuchos::TwoDArray<short>)
Q_DECLARE_METATYPE(Teuchos::TwoDArray<float>)
Q_DECLARE_METATYPE(Teuchos::TwoDArray<double>)
Q_DECLARE_METATYPE(Teuchos::TwoDArray<std::string>)


#endif //OPTIKA_CONFIG_DEFS_HPP
