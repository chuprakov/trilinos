#include "Teuchos_RCP.hpp"
#include "Teuchos_SimpleObjectTable.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"


#ifndef CTRILINOS_TABLE_HPP
#define CTRILINOS_TABLE_HPP


namespace CTrilinos
{

using Teuchos::RCP;
using Teuchos::rcp;


template <class T>
class Table
{
  public:

    /* constructor -- use is_const = true if table will store
     * objects of type const T instead of T */
    Table(CTrilinos_Type_ID_t type, std::string str, bool is_const = false);

    /* destructor */
    ~Table();

    /* retrieve the object */
    const RCP<T> get(CTrilinos_Object_ID_t id);

    /* store an object of type T */
    CTrilinos_Object_ID_t store(T* pobj);

    /* store an object whose base class is T */
    template <class Told>
    CTrilinos_Object_ID_t store(Told* pobj);

    /* store a copy of an object of type T */
    CTrilinos_Object_ID_t storeCopy(T* pobj);

    /* store a copy of an object whose base class is T */
    template <class Told>
    CTrilinos_Object_ID_t storeCopy(Told* pobj);

    /* remove an object from the table and invalidate the id struct */
    int remove(CTrilinos_Object_ID_t * id);

    /* dump the table's contents but keep it's properties */
    void purge();

    /* cast the object to type T and store a copy */
    template <class Told>
    CTrilinos_Object_ID_t cast(const RCP<Told> & rold);

  private:

    /* table for storing objects */
    Teuchos::SimpleObjectTable<T> sot;

    /* properties of the table */
    CTrilinos_Type_ID_t ttype;  /* enum value for stored objects */
    std::string tstr;           /* string form of enum (for exception msgs) */
    bool tconst;                /* if table holds const T */
};


/* constructor -- use is_const = true if table will store
 * objects of type const T instead of T */
template <class T>
Table<T>::Table(CTrilinos_Type_ID_t type, std::string str, bool is_const)
  : ttype(type),
    tconst(is_const)
{
    /* assemble exception error message for future use */
    std::string tmp1 = "Expected type ";
    std::string tmp2 = "(const)";
    std::string tmp3 = "(non-const)";
    if (is_const) {
      tstr = tmp1 + str + tmp2;
    } else {
      tstr = tmp1 + str + tmp3;
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
        throw CTrilinosTypeMismatchError(tstr);
        return Teuchos::null;
    }
}

/* store an object of type T */
template <class T>
CTrilinos_Object_ID_t Table<T>::store(T* pobj)
{
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
    CTrilinos_Object_ID_t id;
    id.type = CT_Invalid_ID;
    id.index = -1;
    id.is_const = tconst;

    T* pnew = dynamic_cast<T*>(pobj);
    if (pnew != NULL) {
        id = store(pnew);
    } else {
        throw CTrilinosTypeMismatchError(tstr);
    }

    return id;
}

/* store a copy of an object of type T */
template <class T>
CTrilinos_Object_ID_t Table<T>::storeCopy(T* pobj)
{
    CTrilinos_Object_ID_t id;
    id.type = CT_Invalid_ID;
    id.index = -1;
    id.is_const = tconst;

    if ((id.index = sot.storeNew(pobj, false)) != -1)
        id.type = ttype;

    return id;
}

/* store a copy of an object whose base class is T */
template <class T>
template <class Told>
CTrilinos_Object_ID_t Table<T>::storeCopy(Told* pobj)
{ /* prevent adding wrong types */
    CTrilinos_Object_ID_t id;
    id.type = CT_Invalid_ID;
    id.index = -1;
    id.is_const = tconst;

    T* pnew = dynamic_cast<T*>(pobj);
    if (pnew != NULL) {
        id = storeCopy(pnew);
    } else {
        throw CTrilinosTypeMismatchError(tstr);
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
        throw CTrilinosTypeMismatchError(tstr);
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


} // namespace CTrilinos


#endif // CTRILINOS_TABLE_HPP


