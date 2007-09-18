/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 *
 *  'TypeList' templates significantly enhanced from
 *  Alexandrescu's "Modern C++ Design" book.
 */

#ifndef util_TypeList_h
#define util_TypeList_h

#include <util/Basics.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

/** @class  TypeList
 *  @brief  Linked list of types.
 */
template<typename ValueType, typename ListType>
struct TypeList {
  /** Value for the current entry */
  typedef ValueType TypeListValue ;

  /** Remainder of the list */
  typedef ListType  TypeListTail ;
};

/** The end of a TypeList is the first encounter with an entry of TypeListEnd.
 *  Thus all entries after a TypeListEnd entry are ignored.  This allows
 *  The MakeTypeList<...> template to define TypeLists of the desired length.
 */
struct TypeListEnd {};

//----------------------------------------------------------------------

/** Length of a TypeList */
template< class ListType>
struct TypeListLength /* { enum { value = <> }; } */ ;

/** Location of the Ith occurance ValueType in the TypeList.  */
template< class ListType, typename ValueType, unsigned I = 0>
struct TypeListIndex /* { enum { value = <> }; } */ ;

/** Count of appearances of ValueType in the TypeList */
template< class ListType, typename ValueType>
struct TypeListCount /* { enum { value = <> }; } */ ;

/** Last ValueType in the TypeList */
template< class ListType >
struct TypeListLast /* { typedef <> type ; } */ ;

/** ValueType and sub-ListType at a location in a TypeList */
template< class ListType, unsigned I>
struct TypeListAt /* { typedef <> type ; typedef <> list_type ; } */ ;

/** TypeList member of a ValueType in a TypeList */
template< class ListType, typename ValueType>
struct TypeListMember /* { typedef <> list_type ; } */ ;

/** Erase type from TypeList at I */
template< class ListType, unsigned I >
struct TypeListEraseAt /* { typedef <> list_type ; } */ ;

/** Check for uniqueness of a TypeList */
template< class ListType >
struct TypeListUnique /* { enum { value = <> }; } */ ;

/** Check if SuperList contains SubList */
template< class SuperList , class SubList >
struct TypeListContains /* { enum { value = <> }; } */ ;

/** Check if ListA is disjoint from ListB */
template< class ListA , class ListB >
struct TypeListDisjoint /* { enum { value = <> }; } */ ;

/** Truncate a TypeList at the first appearance of TypeListEnd */
template<class ListType>
struct TypeListClean /* { typedef <> list_type ; } */ ;

/** Make a TypeList from a sequence of type entries.
 *  Implemented to support list of up to thirtytwo (32) types.
 */
template< typename T0 = TypeListEnd ,
	  typename T1 = TypeListEnd ,
	  typename T2 = TypeListEnd ,
	  typename T3 = TypeListEnd ,
	  typename T4 = TypeListEnd ,
	  typename T5 = TypeListEnd ,
	  typename T6 = TypeListEnd ,
	  typename T7 = TypeListEnd ,
	  typename T8 = TypeListEnd ,
	  typename T9 = TypeListEnd ,
	  typename T10 = TypeListEnd ,
	  typename T11 = TypeListEnd ,
	  typename T12 = TypeListEnd ,
	  typename T13 = TypeListEnd ,
	  typename T14 = TypeListEnd ,
	  typename T15 = TypeListEnd ,
	  typename T16 = TypeListEnd ,
	  typename T17 = TypeListEnd ,
	  typename T18 = TypeListEnd ,
	  typename T19 = TypeListEnd ,
	  typename T20 = TypeListEnd ,
	  typename T21 = TypeListEnd ,
	  typename T22 = TypeListEnd ,
	  typename T23 = TypeListEnd ,
	  typename T24 = TypeListEnd ,
	  typename T25 = TypeListEnd ,
	  typename T26 = TypeListEnd ,
	  typename T27 = TypeListEnd ,
	  typename T28 = TypeListEnd ,
	  typename T29 = TypeListEnd ,
	  typename T30 = TypeListEnd ,
	  typename T31 = TypeListEnd >
struct MakeTypeList
/* { typedef <> type ; enum { length = <> , unique = <> }; } */ ;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Clean: End -> End
// Clean: <End,  Tail> -> End
// Clean: <Value,Tail> -> < Value , Clean<Tail> >

template<>
struct TypeListClean<TypeListEnd> {
  typedef TypeListEnd list_type ;
};

template<typename Tail>
struct TypeListClean< TypeList< TypeListEnd , Tail > > {
  typedef TypeListEnd list_type ;
};

template<class ListType>
struct TypeListClean {
private:
  typedef typename ListType::TypeListValue                 ValueType ;
  typedef typename ListType::TypeListTail                  InputTailType ;
  typedef typename TypeListClean<InputTailType>::list_type TailType ;
public:
  typedef TypeList< ValueType , TailType > list_type ;
};

//----------
// Length: End          -> 0
// Length: <End,  Tail> -> 0
// Length: <Value,Tail> -> 1 + Length<Tail>

template<> struct TypeListLength<TypeListEnd>
{ enum { value = 0 }; };

template<class ListType>
struct TypeListLength< TypeList<TypeListEnd,ListType> >
{ enum { value = 0 }; };

template<class ListType>
struct TypeListLength
{
private: typedef typename ListType::TypeListTail TailType ;
public:  enum { value = 1 + TypeListLength<TailType>::value };
};

//----------
// Index:  < End ,                  ValueType >  ->  -1
// Index:  < List ,                 End >        ->  -1
// Index:  < < End ,       Tail > , ValueType >  ->  -1
// Index:  < < ValueType , Tail > , ValueType >  ->   0
// Index:  < < OtherType , Tail > , ValueType >  ->
//         ( I = Index<List::Tail,ValueType> , I == -1 ? -1 : I + 1 )

template<typename ValueType, unsigned I>
struct TypeListIndex< TypeListEnd , ValueType , I> {
  enum { value = -1 };
  typedef TypeListEnd tail_type ;
};

template<class ListType, unsigned I>
struct TypeListIndex< ListType , TypeListEnd , I > {
  enum { value = -1 };
  typedef TypeListEnd tail_type ;
};

template<class Tail, typename ValueType, unsigned I >
struct TypeListIndex< TypeList<TypeListEnd,Tail> , ValueType, I >
{
  enum { value = -1 };
  typedef TypeListEnd tail_type ;
};

template<typename ValueType, class Tail>
struct TypeListIndex< TypeList<ValueType,Tail> , ValueType , 0 >
{
  enum { value = 0 };
  typedef Tail tail_type ;
};

// Condition: ValueType == ListType::TypeListValue && I matches

template<class ListType, typename ValueType, unsigned I>
struct TypeListIndex
{
private:
  enum { same = SameType<ValueType,typename ListType::TypeListValue>::value };
  enum { ord = I == 0 ? 0 : ( same ? I - 1 : I ) };
  typedef typename ListType::TypeListTail TailType ;
  typedef TypeListIndex< TailType , ValueType , ord > type_list_index ;
  enum { temp = type_list_index::value };
public:
  enum { value = temp == -1 ? -1 : 1 + temp };
  typedef typename type_list_index::tail_type tail_type ;
};

//----------
// Count :  < End , ValueType >  ->  0
// Count :  < List , End >       ->  0
// Count :  < < End ,       Tail > , ValueType >  ->  0
// Count :  < < ValueType , Tail > , ValueType >  ->  Count<Tail,ValueType> + 1
// Count :  < < OtherType , Tail > , ValueType >  ->  Count<Tail,ValueType>

template<typename ValueType>
struct TypeListCount< TypeListEnd , ValueType > { enum { value = 0 }; };

template<class ListType>
struct TypeListCount< ListType , TypeListEnd > { enum { value = 0 }; };

template<class Tail, typename ValueType>
struct TypeListCount< TypeList<TypeListEnd,Tail>,ValueType>
{ enum { value = 0 }; };

template<typename ValueType, class Tail>
struct TypeListCount< TypeList<ValueType,Tail> , ValueType>
{ enum { value = 1 + TypeListCount< Tail , ValueType >::value }; };

template<class ListType, typename ValueType>
struct TypeListCount
{
private: typedef typename ListType::TypeListTail TailType ;
public:  enum { value = TypeListCount< TailType , ValueType >::value };
};

//----------
// At :  < End ,                  0 >  ->  { End , End }
// At :  < End ,                  I >  ->  { End , End }
// At :  < < End ,       Tail > , I >  ->  { End , End }
// At :  < < ValueType , Tail > , 0 >  ->  { ValueType , < ValueType , Tail > }
// At :  < < ValueType , Tail > , I >  ->  At< Tail , I - 1 >

template<>
struct TypeListAt< TypeListEnd, 0>
{
  typedef TypeListEnd type ;
  typedef TypeListEnd list_type ;
};

template<unsigned I>
struct TypeListAt< TypeListEnd, I>
{
  typedef TypeListEnd type ;
  typedef TypeListEnd list_type ;
};

template< class ListType >
struct TypeListAt< ListType , 0 >
{
private:
  typedef typename ListType::TypeListTail Tail ;
public:
  typedef typename ListType::TypeListValue type ;
  typedef TypeList< type , Tail >          list_type ;
};

template<class Tail, unsigned I>
struct TypeListAt< TypeList<TypeListEnd,Tail>, I>
{
  typedef TypeListEnd type ;
  typedef TypeListEnd list_type ;
};

template<class ListType, unsigned I>
struct TypeListAt
{
private:
  typedef typename ListType::TypeListTail Tail ;
  typedef TypeListAt<Tail,I-1>            AtType ;
public:
  typedef typename AtType::type      type ;
  typedef typename AtType::list_type list_type ;
};

//----------
// Last : End -> End
// Last : < ValueType , End >             ->  ValueType
// Last : < ValueType , < End , Tail > >  ->  ValueType
// Last : < ValueType , Tail >            ->  Last< Tail >

template<>
struct TypeListLast< TypeListEnd >
{ typedef TypeListEnd type ; };

template<class ValueType>
struct TypeListLast< TypeList<ValueType,TypeListEnd> >
{ typedef ValueType type ; };

template<class ValueType,class Tail>
struct TypeListLast< TypeList<ValueType,TypeList<TypeListEnd,Tail> > >
{ typedef ValueType type ; };

template<class ValueType, class Tail>
struct TypeListLast< TypeList<ValueType,Tail> >
{ typedef typename TypeListLast<Tail>::type type ; };

//----------
// Member :
// Member :
// Member :
//

template< typename ValueType >
struct TypeListMember< TypeListEnd , ValueType >
{ typedef TypeListEnd list_type ; };

template< class Tail , typename ValueType >
struct TypeListMember< TypeList<TypeListEnd,Tail> , ValueType >
{ typedef TypeListEnd list_type ; };

template< typename ValueType , class ListType>
struct TypeListMember< TypeList<ValueType,ListType> , ValueType >
{ typedef TypeList<ValueType,ListType> list_type ; };

template< class ListType, typename ValueType>
struct TypeListMember
{
  private: typedef typename ListType::TypeListTail Tail ;
  public:  typedef typename TypeListMember<Tail,ValueType>::list_type list_type;
};

//----------
// EraseAt :  < End ,                  0 >  ->  { End , End }
// EraseAt :  < End ,                  I >  ->  { End , End }
// EraseAt :  < < End ,       Tail > , I >  ->  { End , End }
// EraseAt :  < < ListType ,           0 >  ->  { TypeList < ValueType , Tail > }
// EraseAt :  < < ListType , Tail > ,  I >  ->  { EraseAt< Tail , I - 1 > }

template<>
struct TypeListEraseAt< TypeListEnd, 0>
{
  typedef TypeListEnd				list_type ;
};

template<unsigned I>
struct TypeListEraseAt< TypeListEnd, I>
{
  typedef TypeListEnd				list_type ;
};

template<class Tail, unsigned I>
struct TypeListEraseAt< TypeList<TypeListEnd,Tail>, I>
{
  typedef TypeListEnd				list_type ;
};

template< class ListType >
struct TypeListEraseAt< ListType , 0 >
{
private:
  typedef typename ListType::TypeListTail	Tail ;
public:
  typedef Tail					list_type ;
};

template<class ListType, unsigned I>
struct TypeListEraseAt
{
private:
  typedef typename ListType::TypeListTail	Tail ;
  typedef TypeListEraseAt<Tail, I - 1>		EraseAtType ;
public:
  typedef TypeList<typename ListType::TypeListValue, typename EraseAtType::list_type>	list_type ;
};

//----------
// Unique :  End  ->  true
// Unique :  < End , Tail >  -> true
// Unique :  < ValueType , Tail >  ->
//           Index<Tail,ValueType> == -1 && Unique<Tail>

template<>
struct TypeListUnique<TypeListEnd> { enum { value = true }; };

template<class Tail>
struct TypeListUnique< TypeList<TypeListEnd,Tail> >
{ enum { value = true }; };

template< class ListType >
struct TypeListUnique
{
private:
  typedef typename ListType::TypeListValue ValueType ;
  typedef typename ListType::TypeListTail  TailType ;
public:
  // This ValueType does not appear in the remainder of the TypeList and
  // the remainder of the TypeList is also unique.
  enum { value = ( TypeListIndex<TailType,ValueType>::value == -1 ) &&
		   TypeListUnique<TailType>::value };
};

//----------
// Contains : < SuperList , End >  -> true
// Contains : < SuperList , < End , Tail > >  -> true
// Contains : < SuperList , SubList >  ->
//            Index<   SuperList,SubList::Value> != -1 &&
//            Contains<SuperList,SubList::Tail>

template<class SuperList>
struct TypeListContains<SuperList,TypeListEnd>
{ enum { value = true }; };

template<class SuperList,typename Tail>
struct TypeListContains<SuperList,TypeList<TypeListEnd,Tail> >
{ enum { value = true }; };

template<class SuperList, class SubList >
struct TypeListContains
{
private:
  typedef typename SubList::TypeListValue ValueType ;
  typedef typename SubList::TypeListTail  TailType ;
public:
  // The SuperList contains this ValueType and the remainder of the SubList
  enum { value = ( TypeListIndex<SuperList,ValueType>::value != -1 ) &&
		   TypeListContains<SuperList,TailType>::value };
};

//----------
// Disjoint : < ListA , End >  ->  true
// Disjoint : < ListA , < End , Tail > >  ->  true
// Disjoint : < ListA , ListB >  ->
//            Index<   ListA,ListB::Value> == -1 &&
//            Disjoint<ListA,ListB::Tail>

template<class SuperList>
struct TypeListDisjoint<SuperList,TypeListEnd>
{ enum { value = true }; };

template<class SuperList,typename Tail>
struct TypeListDisjoint<SuperList,TypeList<TypeListEnd,Tail> >
{ enum { value = true }; };

template<class ListA, class ListB>
struct TypeListDisjoint
{
private:
  typedef typename ListB::TypeListValue ValueType ;
  typedef typename ListB::TypeListTail  TailType ;
public:
  // ListA does not contain this ValueType and does not contain the remainder
  enum { value = ( TypeListIndex<ListA,ValueType>::value == -1 ) &&
		   TypeListDisjoint<ListA,TailType>::value };
};

//----------------------------------------------------------------------

template< typename T1 ,
	  typename T2 ,
	  typename T3 ,
	  typename T4 ,
	  typename T5 ,
	  typename T6 ,
	  typename T7 ,
	  typename T8 ,
	  typename T9 ,
	  typename T10 ,
	  typename T11 ,
	  typename T12 ,
	  typename T13 ,
	  typename T14 ,
	  typename T15 ,
	  typename T16 ,
	  typename T17 ,
	  typename T18 ,
	  typename T19 ,
	  typename T20 ,
	  typename T21 ,
	  typename T22 ,
	  typename T23 ,
	  typename T24 ,
	  typename T25 ,
	  typename T26 ,
	  typename T27 ,
	  typename T28 ,
	  typename T29 ,
	  typename T30 ,
	  typename T31 >
struct MakeTypeList<TypeListEnd,
		    T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,
		    T11,T12,T13,T14,T15,T16,T17,T18,T19,T20,
		    T21,T22,T23,T24,T25,T26,T27,T28,T29,T30,T31> {
  typedef TypeListEnd type ;
  enum { length = 0 };
  enum { unique = true };
};

template< typename T0 ,
	  typename T1 ,
	  typename T2 ,
	  typename T3 ,
	  typename T4 ,
	  typename T5 ,
	  typename T6 ,
	  typename T7 ,
	  typename T8 ,
	  typename T9 ,
	  typename T10 ,
	  typename T11 ,
	  typename T12 ,
	  typename T13 ,
	  typename T14 ,
	  typename T15 ,
	  typename T16 ,
	  typename T17 ,
	  typename T18 ,
	  typename T19 ,
	  typename T20 ,
	  typename T21 ,
	  typename T22 ,
	  typename T23 ,
	  typename T24 ,
	  typename T25 ,
	  typename T26 ,
	  typename T27 ,
	  typename T28 ,
	  typename T29 ,
	  typename T30 ,
	  typename T31 >
struct MakeTypeList {
  typedef typename TypeListClean<
	  TypeList< T0 ,
	  TypeList< T1,
	  TypeList< T2 ,
	  TypeList< T3 ,
	  TypeList< T4 ,
	  TypeList< T5 ,
	  TypeList< T6 ,
	  TypeList< T7 ,
	  TypeList< T8 ,
	  TypeList< T9 ,
	  TypeList< T10 ,
	  TypeList< T11 ,
	  TypeList< T12 ,
	  TypeList< T13 ,
	  TypeList< T14 ,
	  TypeList< T15 ,
	  TypeList< T16 ,
	  TypeList< T17 ,
	  TypeList< T18 ,
	  TypeList< T19 ,
	  TypeList< T20 ,
	  TypeList< T21 ,
	  TypeList< T22 ,
	  TypeList< T23 ,
	  TypeList< T24 ,
	  TypeList< T25 ,
	  TypeList< T26 ,
	  TypeList< T27 ,
	  TypeList< T28 ,
	  TypeList< T29 ,
	  TypeList< T30 ,
	  TypeList< T31 ,
	  TypeListEnd > > > > > > > > > > > > > > > >
		      > > > > > > > > > > > > > > > >
  >::list_type type ;
  enum { length = TypeListLength<type>::value };
  enum { unique = TypeListUnique<type>::value };
};

}

#endif

