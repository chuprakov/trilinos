#include "Teuchos_RCP.hpp"
#include "Teuchos_SimpleObjectTable.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_exceptions.hpp"


#ifndef CTRILINOS_TABLE_HPP
#define CTRILINOS_TABLE_HPP


namespace CTrilinos
{

using Teuchos::RCP;


template <class T>
class Table
{
  public:

    Table(CTrilinos_Type_ID_t type, std::string str);

    ~Table();

    const RCP<T> get(CTrilinos_Object_ID_t id);

    CTrilinos_Object_ID_t store(T* pobj);

    CTrilinos_Object_ID_t store(const T* pobj);

    template <class Told>
    CTrilinos_Object_ID_t store(Told* pobj);

    int remove(CTrilinos_Object_ID_t * id);

    void purge();

    template <class Told>
    CTrilinos_Object_ID_t cast(const RCP<Told> & rold);

  private:

    Teuchos::SimpleObjectTable<T> sot;

    CTrilinos_Type_ID_t ttype;
    std::string tstr;

};

template <class T>
Table<T>::Table(CTrilinos_Type_ID_t type, std::string str)
{
  ttype = type;
  tstr = str;
}

template <class T>
Table<T>::~Table()
{
  purge();
}

template <class T>
const RCP<T> Table<T>::get(CTrilinos_Object_ID_t id)
{
    if (id.type == ttype) {
        return sot.getRCP(id.index);
    } else {
        throw CTrilinosTypeMismatchError(tstr);
        return Teuchos::null;
    }
}

template <class T>
CTrilinos_Object_ID_t Table<T>::store(T* pobj)
{
    CTrilinos_Object_ID_t id;
    id.type = CT_Invalid_ID;
    id.index = -1;

    if ((id.index = sot.storeNew(pobj)) != -1)
        id.type = ttype;

    return id;
}

template <class T>
CTrilinos_Object_ID_t Table<T>::store(const T* pobj)
{
    return store(const_cast<T*>(pobj));
}

template <class T>
template <class Told>
CTrilinos_Object_ID_t Table<T>::store(Told* pobj)
{ /* prevent adding wrong types */
    T* pbase = dynamic_cast<T*>(pobj);
    if (pbase != NULL) {
        return store(pbase);
    } else {
        throw CTrilinosTypeMismatchError(tstr);
    }

    CTrilinos_Object_ID_t id;
    id.type = CT_Invalid_ID;
    id.index = -1;

    return id;
}

template <class T>
int Table<T>::remove(CTrilinos_Object_ID_t * id)
{
    int ret = -1;

    if (id->type == ttype) {
        ret = (sot.removeRCP(id->index) < 0 ? -1 : 0);
        if (ret == 0) id->type = CT_Invalid_ID;
    } else {
        throw CTrilinosTypeMismatchError(tstr);
    }

    return ret;
}

template <class T>
void Table<T>::purge()
{
    sot.purge();
}

template <class T>
template <class Told>
CTrilinos_Object_ID_t Table<T>::cast(const RCP<Told> & rold)
{
    CTrilinos_Object_ID_t newid;
    newid.type = CT_Invalid_ID;
    newid.index = -1;

    newid.index = sot.storeCastedRCP(rold);

    if (newid.index != -1)
        newid.type = ttype;

    return newid;
}

} // namespace CTrilinos


#endif // CTRILINOS_TABLE_HPP


