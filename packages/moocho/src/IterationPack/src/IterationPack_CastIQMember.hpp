// ///////////////////////////////////////////////////////////
// CastIQMember.h

#ifndef CAST_IQ_MEMBER_H
#define CAST_IQ_MEMBER_H

#include <limits.h>

#include <typeinfo>

#include "AlgorithmState.h"
#include "IterQuantityAccess.h"

namespace GeneralIterationPack {

///
/** Base class for some of the implementation features.
  *
  * This class is included to avoid code blot with the templates.
  */
class CastIQMemberBase {
public:
	///
	const std::string& iq_name() const;
protected:
	///
	CastIQMemberBase( const std::string iq_name );
	///
	const std::string					iq_name_;
	///
	mutable AlgorithmState::iq_id_type	iq_id_;
	///
	void cache_iq_id( const AlgorithmState& s ) const;
	///
	void throw_cast_error( const std::string& iqa_name ) const;
private:
	enum { NOT_SET_YET = INT_MAX };
	CastIQMemberBase(); // not defined and not to be called.
};	// end class CastIQMemberBase

///
/** Template class to be used to lookup an interation quantity
  * , cast it to an IterQuantityAccess<T> and cache the iq_id
  *  for fast access later.
  *
  * The idea is that a Step class can create a data member
  * of this class and then access and interation quantity
  * and have the iq_id looked up the first time.
  * 
  * The best way to use this class to access an iteration quantity
  * is for a header file to be created for the iteration quantity (or several
  * iteration quantityies) that contains a new class.  For example, suppose
  * we have an double object that we want to add as an iteration quantity
  * with the name "x_step".  For this we might create an header file like:
  * 
  \begin{verbose}
    
    // /////////////////////////////////////////////////////////////////
    // x_step_iter_quantity.h
    
    #include "CastIQMember.h"
    
    class x_step_iq_member : public CastIQMember<double> {
    public:
        x_step_iq_member() : CastIQMember<double>("x_step") {}
    }
    
  \end{verbose} 
  *
  * Now lets suppose we have two step classes that need to access this
  * iteration quantity.  These step classes would each include
  * a data member of this new class.  For example, these
  * step classes might be implemented as:
  * 
  \begin{verbose}
   
    // /////////////////////////////////////////////////////////////////
    // MyStep1.h
    
    #include "x_step_iter_quantity.h"
    
    class MyStep1 : public AlgorithmStep {
    public:
        bool do_step( algo, ... )
        {
        	AlgorithmState &s = algo.state();
        	x_step_(s).set_k(0) = 5.0;
        }
    private:
    	x_step_iq_member x_step_;
    }
   
    // /////////////////////////////////////////////////////////////////
    // MyStep2.h
    
    #include "x_step_iter_quantity.h"
    
    class MyStep2 : public AlgorithmStep {
    public:
        bool do_step( algo, ... )
        {
        	AlgorithmState &s = algo.state();
        	double x_step = x_step_(s).get_k(0);
        	cout << "\nx_step = " << x_step << std::endl;
        }
    private:
    	x_step_iq_member x_step_;
    }
    
  \end{verbose}
  *
  * In the above example, on O(s.num_iter_quantities()) search for the iq_id
  * would only be performed the first time x_step_(s) was called by each
  * step object.  At later iterations, the cached iq_id would be used to access
  * the iteration quantity and the only price you would pay above a few
  * O(1) function calls in is an O(1) dynamic cast.
  */
template < class T >
class CastIQMember : public CastIQMemberBase {
public:
	/// Construct with the name of an iteration quantity.
	CastIQMember( const std::string iq_name );
	///
	/** Get the iteration quantity from an AlgorithmState object.
	  *
	  * If the iteration quantity of the name iq_namt does not
	  * exist then a AlgorithmState::DoesNotExist exception
	  * will be thrown.  If the type of the iteration quantity
	  * is not of the type IterQuantityAcess<T> (as determined
	  * by dynamic_cast<T>) then the exception InvalidTypeCastException:
	  * will be thrown with a helpful error message.
	  */
	IterQuantityAccess<T>& operator()( AlgorithmState& s ) const;
	///
	const IterQuantityAccess<T>& operator()( const AlgorithmState& s ) const;
private:
	CastIQMember();	// not defined and not to be called
};	// end class CastIQMember<T>

// //////////////////////////////////////////
// Definition of template members

template < class T >
CastIQMember<T>::CastIQMember( const std::string iq_name )
	:  CastIQMemberBase(iq_name)
{}

template < class T >
IterQuantityAccess<T>&
CastIQMember<T>::operator()( AlgorithmState& s ) const
{
	cache_iq_id(s);
	IterQuantityAccess<T>
		*p = dynamic_cast<IterQuantityAccess<T>*>( &s.iter_quant( iq_id_ ) );
	if( !p )
		throw_cast_error(typeid(T).name());
	return *p;	
}

template < class T >
const IterQuantityAccess<T>&
CastIQMember<T>::operator()( const AlgorithmState& s ) const
{
	cache_iq_id(s);
	const IterQuantityAccess<T>
		*p = dynamic_cast<const IterQuantityAccess<T>*>( &s.iter_quant( iq_id_ ) );
	if( !p )
		throw_cast_error(typeid(T).name());
	return *p;	
}

}	// namespace GeneralIterationPack

#endif	// CAST_IQ_MEMBER_H
