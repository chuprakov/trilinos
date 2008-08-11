// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_PARAMETERLIBRARYBASE_HPP
#define SACADO_PARAMETERLIBRARYBASE_HPP

#include "Teuchos_Array.hpp"

#include "Sacado_ParameterFamilyBase.hpp"
#include "Sacado_ParameterVectorBase.hpp"

namespace Sacado {

  using std::string;

  /*! 
   * \brief Class to provide a centralized library for setting/retrieving 
   * numerical parameter values.
   */
  template <typename FamilyType, template<typename> class EntryType>
  class ParameterLibraryBase  {

  protected:

    //! Map of all parameter families
    typedef std::map<string, Teuchos::RCP<FamilyType> > FamilyMap;

  public:

    //! Iterator typename
    typedef typename FamilyMap::iterator iterator;

    //! Const iterator typename
    typedef typename FamilyMap::const_iterator const_iterator;
  
    //! Default constructor
    ParameterLibraryBase();

    //! Destructor
    virtual ~ParameterLibraryBase();

    //! Determine if parameter of name \em name is in the library
    bool isParameter(const std::string& name) const;

    //! Determine if parameter of name \em name has type \em type
    template <typename ValueType>
    bool isParameterForType(const std::string& name) const;

    //! Create a new scalar parameter family
    /*!
     * Returns true if successful in adding family to library, false 
     * otherwise.
     */
    bool addParameterFamily(const std::string& name, bool supports_ad, 
			    bool supports_analytic);
    

    //! Add a new scalar parameter using custom entry
    /*!
     * Returns true if successful in adding entry to library, false 
     * otherwise.
     */
    template <typename ValueType>
    bool addEntry(const std::string& name, 
		  const Teuchos::RCP< EntryType<ValueType> >& entry);

    //! Return parameter entry
    template <typename ValueType>
    Teuchos::RCP< EntryType<ValueType> >
    getEntry(const std::string& name);

    //! Return number of parameters in library
    unsigned int size() const { return library.size(); }

    //! Iterator pointing at beginning of library
    iterator begin() { return library.begin(); }

    //! Iterator pointing at beginning of library
    const_iterator begin() const { return library.begin(); }

    //! Iterator pointing at end of library
    iterator end() { return library.end(); }

    //! Iterator pointing at end of library
    const_iterator end() const { return library.end(); }

    //! Fill a vector with the supplied parameter names and values
    template <typename BaseValueType>
    void
    fillVector(const Teuchos::Array<std::string>& names,
	       const Teuchos::Array<BaseValueType>& values,
	       ParameterVectorBase<FamilyType,BaseValueType>& pv);

  private:

    //! Private to prohibit copying
    ParameterLibraryBase(const ParameterLibraryBase&);

    //! Private to prohibit copying
    ParameterLibraryBase& operator = (const ParameterLibraryBase&);

  protected:

    //! Scalar parameter library
    FamilyMap library;
  };
}

// Include template definitions
#include "Sacado_ParameterLibraryBaseImp.hpp"

#endif
