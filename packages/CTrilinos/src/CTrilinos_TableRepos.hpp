
/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */


/*! @file CTrilinos_TableRepos.hpp
 * @brief Central table repository for CTrilinos. */


#ifndef CTRILINOS_TABLEREPOS_HPP
#define CTRILINOS_TABLEREPOS_HPP


#include "CTrilinos_config.h"


#include <string>

#include "CTrilinos_enums.h"
#include "CTrilinos_TableRepos_Utils.hpp"
#include "CTrilinos_utils_templ.hpp"
#include "Teuchos_RCP.hpp"


namespace CTrilinos {


/*! Central table repository for storing Trilinos objects. */
class TableRepos
{
  public:

    /*! constructor */
    TableRepos()
      : tabPoly(),
        tabNonPoly()
    { }

    /*! destructor */
    ~TableRepos()
    { }

    /*! retrieve the object */
    template <class T, typename S>
    const Teuchos::RCP<T> get(S id);

    /*! retrieve the object */
    template <class T, typename S>
    const Teuchos::RCP<const T> getConst(S id);

    /*! store a non-const RCP to object of type T */
    template <class T, typename S>
    S store(T* pobj, bool owned = true);

    /*! store a const RCP to object of type T */
    template <class T, typename S>
    S store(const T* pobj, bool owned = false);

    /*! remove an object from the table and invalidate the id struct */
    template <typename S>
    void remove(S * id);

    /*! dump the table's contents but keep it's properties */
    void purge()
    {
      tabPoly.purge();
      tabNonPoly.purge();
    }

    /*! create an alias for the object in another table */
    CTrilinos_Universal_ID_t alias(CTrilinos_Universal_ID_t aid,
        CTrilinos_Table_ID_t tab, bool keepold = true)
    {
	if (poly_checker_runtime(aid)) {
	    return tabPoly.alias(aid, tab, keepold);
	} else {
	    throw std::string("cannot alias a non-polymorphic type");
	}
    }

  private:

    TableReposPoly<true> tabPoly;
    TableReposPoly<false> tabNonPoly;

    template <bool POLY>
    TableReposPoly<POLY> &getPolyTable();

    template <class T>
    TableReposPoly<poly_checker<T>::is_poly> &getTable()
    { return getPolyTable<poly_checker<T>::is_poly>(); }

};

template <bool POLY>
TableReposPoly<POLY> & TableRepos::getPolyTable()
{ return tabPoly; }

template <>
TableReposPoly<false> & TableRepos::getPolyTable()
{ return tabNonPoly; }

template <class T, typename S>
S TableRepos::store(T* pobj, bool owned)
{
    return concreteType<S>(getTable<T>().store(pobj, owned));
}

template <class T, typename S>
S TableRepos::store(const T* pobj, bool owned)
{
    return concreteType<S>(getTable<T>().store(pobj, owned));
}

template <class T, typename S>
const Teuchos::RCP<T> TableRepos::get(S id)
{
    if (id.is_const)
        throw std::string("const violation");

    CTrilinos_Universal_ID_t aid = abstractType<S>(id);

    return getTable<T>().get<T>(aid);
}

template <class T, typename S>
const Teuchos::RCP<const T> TableRepos::getConst(S id)
{
    if (!id.is_const)
        return get<T,S>(id);

    CTrilinos_Universal_ID_t aid = abstractType<S>(id);

    return getTable<T>().getConst<T>(aid);
}

template <typename S>
void TableRepos::remove(S * id)
{
    CTrilinos_Universal_ID_t aid = abstractType<S>(*id);

    if (poly_checker_runtime(aid)) {
        tabPoly.remove(&aid);
    } else {
        tabNonPoly.remove(&aid);
    }
    *id = concreteType<S>(aid);
}


} // namespace CTrilinos


#endif

