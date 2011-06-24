// @HEADER
// ***********************************************************************
// 
//         Optika: A Tool For Developing Parameter Obtaining GUIs
//                Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, with Sandia Corporation, the 
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Kurtis Nusbaum (klnusbaum@gmail.com) 
// 
// ***********************************************************************
// @HEADER
#ifndef OPTIKA_ARRAYHELPERFUNCTIONS_HPP_
#define OPTIKA_ARRAYHELPERFUNCTIONS_HPP_
#include <QStringList>
#include <QVariant>
#include "Optika_Types.hpp"
#include "Optika_ConfigDefs.hpp"
/*! \file Opitka_ArrayHelperFunctions.hpp
    \brief Helper functions used for manipulating
    and querying Arrays
*/
namespace Optika{

/**
 * \brief Determines the type of array stored in a parameter.
 *
 * @param parameter The parameter whose array type is in question.
 * @return A QString containing the type of array in the parameter.
 */
QString determineArrayType(RCP<const ParameterEntry> parameter);

/**
 * \brief Creates a QVariant containing the array that is in
 * arrayEntry.
 *
 * @param arrayEntry The parameter entry containing the array.
 * @param type The array's template type.
 * @return A QVariant containing an Array that is equal to the Array in array entry.
 */
QVariant arrayEntryToVariant(
  RCP<const ParameterEntry> arrayEntry, QString type);

/**
 * \brief Given a type string, determines the template type of the Array.
 *
 * @param itemType The type string describing the array.
 * @return The template type of the Array.
 */
QString getArrayType(QString itemType);

/**
 * \brief Determines wether or no the array inside a ParameterEntry is empty.
 *
 * @param arrayEntry The parameter entry containging the array.
 * @param type The template type of the array.
 * @return True if the array is empty, false otherwise.
 */
bool isArrayEmpty(RCP<const ParameterEntry> arrayEntry, QString type);

}
#endif /* OPTIKA_ARRAYHELPERFUNCTIONS_HPP_ */
