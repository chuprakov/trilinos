#ifndef CTRILINOS_TABLE_HPP
#define CTRILINOS_TABLE_HPP


#include <string>
#include <typeinfo>

#include "Teuchos_RCP.hpp"
#include "Teuchos_SimpleObjectTable.hpp"
#include "Teuchos_Exceptions.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"
#include "CTrilinos_utils.hpp"


namespace CTrilinos
{

using Teuchos::RCP;
using Teuchos::rcp;

/* stringify the enum name -- defined in CTrilinos_utils.cpp */
std::string enum2str( CTrilinos_Type_ID_t ty );

template <class T>
class Table
{
  public:

    /* constructor -- use is_const = true if table will store
     * objects of type const T instead of T */
    Table(CTrilinos_Type_ID_t type, std::string str, boolean is_const = FALSE);

    /* destructor */
    ~Table();

    /* retrieve the object */
    const RCP<T> get(CTrilinos_Object_ID_t id);

    /* store an object of type T */
    CTrilinos_Object_ID_t store(T* pobj);

    /* store an object whose base class is T */
    template <class Told>
    CTrilinos_Object_ID_t store(Told* pobj);

    /* store a non-owned object of type T */
    CTrilinos_Object_ID_t storeShared(T* pobj);

    /* store a non-owned object whose base class is T */
    template <class Told>
    CTrilinos_Object_ID_t storeShared(Told* pobj);

    /* remove an object from the table and invalidate the id struct */
    int remove(CTrilinos_Object_ID_t * id);

    /* dump the table's contents but keep it's properties */
    void purge();

    /* cast the object to type T and store a copy */
    template <class Told>
    CTrilinos_Object_ID_t cast(const RCP<Told> & rold);

  private:

    /* build full exception msg on the fly */
    std::string typeMismatchMsg(CTrilinos_Object_ID_t id, std::string act);

    /* build full exception msg on the fly */
    std::string badCastMsg(std::string type, std::string act);

    /* table for storing objects */
    Teuchos::SimpleObjectTable<T> sot;

    /* properties of the table */
    CTrilinos_Type_ID_t ttype;  /* enum value for stored objects */
    std::string tstr;           /* string for exception msgs */
    std::string tstr2;          /* string for exception msgs */
    boolean tconst;                /* if table holds const T */
};


/* constructor -- use is_const = true if table will store
 * objects of type const T instead of T */
template <class T>
Table<T>::Table(CTrilinos_Type_ID_t type, std::string str, boolean is_const)
  : ttype(type),
    tconst(is_const)
{
    /* assemble exception error message for future use */
    std::string tmp1 = "[CTrilinos::Table]: Expected type ";
    std::string tmp2 = " (const)";
    std::string tmp3 = " (non-const)";
    std::string tmp4 = "; found ";
    std::string tmp5 = typeid(T).name();
    if (is_const) {
      tstr = tmp1 + str + tmp2 + tmp4;
      tstr2 = tmp1 + tmp5 + tmp2 + tmp4;
    } else {
      tstr = tmp1 + str + tmp3 + tmp4;
      tstr2 = tmp1 + tmp5 + tmp3 + tmp4;
    }
}

/* destructor */
template <class T>
Table<T>::~Table()
{
  purge();
}

/* retrieve the object */
template <class T>
const RCP<T> Table<T>::get(CTrilinos_Object_ID_t id)
{
    if ((id.type == ttype) && (id.is_const == tconst)) {
        return sot.getRCP(id.index);
    } else {
        throw CTrilinosTypeMismatchError(typeMismatchMsg(id, std::string("get()")));
        return Teuchos::null;
    }
}

/* store an object of type T */
template <class T>
CTrilinos_Object_ID_t Table<T>::store(T* pobj)
{
    if (pobj == NULL)
        throw Teuchos::NullReferenceError("[CTrilinos::Table]: Cannot store NULL pointer");

    CTrilinos_Object_ID_t id;
    id.type = CT_Invalid_ID;
    id.index = -1;
    id.is_const = tconst;

    if ((id.index = sot.storeNew(pobj)) != -1)
        id.type = ttype;

    return id;
}

/* store an object whose base class is T */
template <class T>
template <class Told>
CTrilinos_Object_ID_t Table<T>::store(Told* pobj)
{ /* prevent adding wrong types */
    if (pobj == NULL)
        throw Teuchos::NullReferenceError("[CTrilinos::Table]: Cannot store NULL pointer");

    CTrilinos_Object_ID_t id;
    id.type = CT_Invalid_ID;
    id.index = -1;
    id.is_const = tconst;

    T* pnew = dynamic_cast<T*>(pobj);
    if (pnew != NULL) {
        id = store(pnew);
    } else {
        throw CTrilinosTypeMismatchError(badCastMsg(typeid(*pobj).name(), std::string("store()")));
    }

    return id;
}

/* store a non-owned object of type T */
template <class T>
CTrilinos_Object_ID_t Table<T>::storeShared(T* pobj)
{
    if (pobj == NULL)
        throw Teuchos::NullReferenceError("[CTrilinos::Table]: Cannot store NULL pointer");

    CTrilinos_Object_ID_t id;
    id.type = CT_Invalid_ID;
    id.index = -1;
    id.is_const = tconst;

    if ((id.index = sot.storeNew(pobj, false)) != -1)
        id.type = ttype;

    return id;
}

/* store a non-owned object whose base class is T */
template <class T>
template <class Told>
CTrilinos_Object_ID_t Table<T>::storeShared(Told* pobj)
{ /* prevent adding wrong types */
    if (pobj == NULL)
        throw Teuchos::NullReferenceError("[CTrilinos::Table]: Cannot store NULL pointer");

    CTrilinos_Object_ID_t id;
    id.type = CT_Invalid_ID;
    id.index = -1;
    id.is_const = tconst;

    T* pnew = dynamic_cast<T*>(pobj);
    if (pnew != NULL) {
        id = storeShared(pnew);
    } else {
        throw CTrilinosTypeMismatchError(badCastMsg(typeid(*pobj).name(), std::string("storeShared()")));
    }

    return id;
}

/* remove an object from the table and invalidate the id struct */
template <class T>
int Table<T>::remove(CTrilinos_Object_ID_t * id)
{
    int ret = -1;

    if ((id->type == ttype) && (id->is_const == tconst)) {
        ret = (sot.removeRCP(id->index) < 0 ? -1 : 0);
        if (ret == 0) id->type = CT_Invalid_ID;
    } else {
        throw CTrilinosTypeMismatchError(typeMismatchMsg(*id, std::string("remove()")));
    }

    return ret;
}

/* dump the table's contents but keep it's properties */
template <class T>
void Table<T>::purge()
{
    sot.purge();
}

/* cast the object to type T and store a copy */
template <class T>
template <class Told>
CTrilinos_Object_ID_t Table<T>::cast(const RCP<Told> & rold)
{
    CTrilinos_Object_ID_t newid;
    newid.type = CT_Invalid_ID;
    newid.index = -1;
    newid.is_const = tconst;

    newid.index = sot.storeCastedRCP(rold);

    if (newid.index != -1)
        newid.type = ttype;

    return newid;
}

/* build full exception msg on the fly */
template <class T>
std::string Table<T>::typeMismatchMsg(CTrilinos_Object_ID_t id, std::string act)
{
    std::string s = enum2str(id.type);

    if (id.is_const) {
      std::string tmp = " (const)";
      s += tmp;
    } else {
      std::string tmp = " (non-const)";
      s += tmp;
    }
    std::string tmp = " when attempting to " + act;
    s += tmp;

    return (tstr + s);
}

/* build full exception msg on the fly */
template <class T>
std::string Table<T>::badCastMsg(std::string type, std::string act)
{
    std::string s = type + " when attempting to " + act;

    return (tstr2 + s);
}


} // namespace CTrilinos


#endif // CTRILINOS_TABLE_HPP


