#ifndef OPTIKA_CONFIG_DEFS_HPP
#define OPTIKA_CONFIG_DEFS_HPP
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_StandardDependencies.hpp"
#include "Teuchos_DependencySheet.hpp"
#include "Teuchos_Array.hpp"
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
  using Teuchos::ArrayNumberValidator;
  using Teuchos::ArrayFileNameValidator;
  using Teuchos::DependencySheet;
  using Teuchos::Dependency;
  using Teuchos::VisualDependency;
  using Teuchos::any;
  using Teuchos::any_cast;
  using Teuchos::XMLParameterListWriter;
  using Teuchos::XMLObject;
  using Teuchos::getValue;
  using Teuchos::getParametersFromXmlFile;


} //namespace Optika

#endif //OPTIKA_CONFIG_DEFS_HPP
