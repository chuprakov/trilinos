// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER


#ifndef TEUCHOS_XMLPARAMETERLISTEXCEPTIONS_HPP_
#define TEUCHOS_XMLPARAMETERLISTEXCEPTIONS_HPP_

/*! \file Teuchos_XMLParameterListExceptions.hpp
 * \brief A collection of Exceptions that can be potentially
 * thrown when converting a ParameterList to and from XML
 */
#include <stdexcept>

namespace Teuchos {

/**
 * \brief Thrown when an appropriate ParameterEntryXMLConverter
 * can't be found.
 */
class CantFindParameterEntryConverterException : public std::logic_error{

public: 

  /**
   * \brief Constructs an CantFindParameterEntryConverterException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  CantFindParameterEntryConverterException(const std::string& what_arg):std::logic_error(what_arg){}

};

/** \brief Constructs a CantFindParameterEntryConverterException */

/**
 * \brief Thrown when two validators in an XML file have the same ID.
 */
class DuplicateValidatorIDsException : public std::logic_error{

public: 

  /**
   * \brief Constructs an DuplicateValidatorIDsException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  DuplicateValidatorIDsException(const std::string& what_arg):std::logic_error(what_arg){}

};

/**
 * \brief Thrown when two parameters in an XML file have the same ID.
 */
class DuplicateParameterIDsException : public std::logic_error{

public: 

  /**
   * \brief Constructs an DuplicateParameterIDsException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  DuplicateParameterIDsException(const std::string& what_arg):std::logic_error(what_arg){}

};

/**
 * \brief Thrown when a bad validator xml converter is used.
 */
class BadValidatorXMLConverterException : public std::logic_error{

public: 
  /**
   * \brief Constructs an BadValidatorXMLConverterException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  BadValidatorXMLConverterException(const std::string& what_arg):std::logic_error(what_arg){}

};


/**
 * \brief Thrown when the ValidatorXMLConverterDB can't find an 
 * appropriate converter.
 */
class CantFindValidatorConverterException : public std::logic_error{

public: 
  /**
   * \brief Constructs a CantFindValidatorConverterException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  CantFindValidatorConverterException(const std::string& what_arg):std::logic_error(what_arg){}

};


/**
 * \brief Thrown when a converter is being used to convert either and XML tag or
 * ParameterEntry with an innappropriate type.
 */
class BadParameterEntryXMLConverterTypeException : public std::logic_error{

public: 

  /**
   * \brief Constructs a BadParmaeterEntryXMLConverterTypeException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  BadParameterEntryXMLConverterTypeException(const std::string& what_arg):std::logic_error(what_arg){}

};



/**
 * \brief Thrown when a parameter entry tag is missing it's value attribute.
 */
class NoValueAttributeExecption : public std::logic_error{
public: 
  /**
   * \brief Constructs a NoValueAttributeExecption.
   *
   * @param what_arg The error message to be associated with this error.
   */
  NoValueAttributeExecption(const std::string& what_arg):std::logic_error(what_arg){}
};


/**
 * \brief Thrown when a parameter entry tag is missing it's type attribute.
 */
class NoTypeAttributeExecption : public std::logic_error{
public: 
  /**
   * \brief Constructs a NoTypeAttributeExecption.
   *
   * @param what_arg The error message to be associated with this error.
   */
  NoTypeAttributeExecption(const std::string& what_arg):std::logic_error(what_arg){}
};

/**
 * \brief Thrown when a parameter entry tag is missing it's name attribute.
 */
class NoNameAttributeExecption : public std::logic_error{
public: 
  /**
   * \brief Constructs a NoNameAttributeExecption.
   *
   * @param what_arg The error message to be associated with this error.
   */
  NoNameAttributeExecption(const std::string& what_arg):std::logic_error(what_arg){}
};

/**
 * \brief Thrown when an element inside a parameter list is bad.
 */
class BadParameterListElementException : public std::logic_error{
public: 
  /**
   * \brief Constructs a BadParameterListElementException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  BadParameterListElementException(const std::string& what_arg):std::logic_error(what_arg){}
};

/**
 * \brief Thrown when the root xml tag for a parameter list is incorrect.
 */
class BadXMLParameterListRootElementException : public std::logic_error{
public: 
  /**
   * \brief Constructs a BadXMLParameterListRootElementException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  BadXMLParameterListRootElementException(const std::string& what_arg):std::logic_error(what_arg){}
};

/**
 * \brief Thrown when the xml tag for a parameter entry has a reference to
 * a non-existant validator.
 */
class MissingValidatorDefinitionException : public std::logic_error{
public: 
  /**
   * \brief Constructs a MissingValidatorDefinitionException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  MissingValidatorDefinitionException(const std::string& what_arg):std::logic_error(what_arg){}
};

/**
 * \brief Thrown when xml tag is encountered that is either unrecognized or 
 * inappropriate for a given context.
 */
class BadTagException : public std::logic_error{
public: 
  /**
   * \brief Constructs a MissingValidatorDefinitionException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  BadTagException(const std::string& what_arg):std::logic_error(what_arg){}
};


} // namespace Teuchos
#endif //TEUCHOS_XMLPARAMETERLISTEXCEPTIONS_HPP_
