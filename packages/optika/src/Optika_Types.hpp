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
#ifndef OPTIKA_TYPES_HPP
#define OPTIKA_TYPES_HPP
#include <QString>

/*! \file Optika_Types.hpp
    \brief A collection of types
    that Optika uses to to keep track
    of ParameterEntry types.
*/
namespace Optika{
/**
 * The id's for the supported Optika types.
 */

/** * \brief .  */
extern QString listId;
/** * \brief .  */
extern QString intId;
/** * \brief .  */
extern QString shortId;
//extern QString longlongId = QString("longlong");
/** * \brief .  */
extern QString doubleId;
/** * \brief .  */
extern QString floatId;
/** * \brief .  */
extern QString boolId;
/** * \brief .  */
extern QString stringId;
/** * \brief .  */
extern QString arrayId;
/** * \brief .  */
extern QString unrecognizedId;
}
#endif //OPTIKA_TYPES_HPP
