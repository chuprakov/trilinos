// /////////////////////////////////////////////////////////////////////////////////////
// IterQuantity.h
//
// Change Log:
//	11/18/99:
//		* last_updated() Added
//		* set_not_updated(k) Added

#ifndef ITER_QUANTITY_H
#define ITER_QUANTITY_H

#include <stdexcept>
#include <string>
#include <sstream>
#include <iomanip>
#include <limits>

#include "GeneralIterationPackTypes.h"

namespace GeneralIterationPack {

///
/** Iterface for information about Iteration Quantities.
  *
  * This class provides and interface to all concrete types of iteration quantities
  * and provides all of the services except storage access.
  *
  * ToDo: 7/27/98: Finish the documentation and give examples.
  */
class IterQuantity {
public:

	/// Constant for value returned when no iteration quantity has been updated.
	enum { NONE_UPDATED = INT_MIN };

	/// Thrown memory if attempted to be set that storage can not be allocated to.
	class NoStorageAvailable : public std::logic_error
	{public: NoStorageAvailable(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Thrown when memory access is attempted when it has not yet been updated.
	class QuanityNotSet : public std::logic_error
	{public: QuanityNotSet(const std::string& what_arg) : std::logic_error(what_arg) {}};

	// Virtual destructor
	virtual ~IterQuantity() {}

	// Clone this iteration quantity
	virtual IterQuantity* clone() const = 0;

	/// Return the name (zero terminated string) of this quantity.
	virtual const char* name() const = 0; 
	
	///
	/** Determine if there is storage advailable for the k #offset# iteration quanity.
	  *
	  * If this function returns true then #set_f(offset)# can be called to set
	  * the quanity for the kth iteration or #get_k(offset)# can be called if
	  * #updated_k(offset)# is already true.  
	  */
	virtual bool has_storage_k(int offset) const = 0;
	
	///
	/** Determine if the quanity for the k #offset# iteration has been accessed
	  * by a call to #set_k()#.
	  *
	  * This function does not confirm that the k #offset# quanity has been
	  * set to a meaningfull value, only that #set_k()# was called
	  * to get a reference to that quanity and #get_k()# can be called to
	  * get that same reference.
	  */
	virtual bool updated_k(int offset) const = 0;

	///
	/** Return the highest k such #updated_k(k)# returns true.
	  *
	  * If no #updated_k(k)# equals false for all k then this function
	  * will return NONE_UPDATED.
	  */
	virtual int last_updated() const = 0;

	///
	/** Causes #updated_k(k)# to return false.
	  *
 	  * Preconditions:\begin{itemize}
	  * \item #updated_k(offset) == true# [throw #QuanityNotSet#]
	  * \end{itemize} 
	  *
 	  * Postconditions:\begin{itemize}
	  * \item #updated_k(offset) == false#
	  * \end{itemize} 
	  */
	virtual void set_not_updated(int offset) = 0;

	/// Assert #has_storage_k(offset) == true# (throw #NoStorageAvailable#).
	void assert_has_storage_k(int offset) const {
		if(!has_storage_k(offset)) {
			std::ostringstream omsg;
			omsg	<< "IterQuantity::assert_has_storage_k(offset) : There is not storage available for "
					<< name() << " for the k"	<< std::showpos << offset << " iteration";
			throw NoStorageAvailable(omsg.str());
		}
	}

	/// Assert updated_k(offset) == true# (throw QuanityNotSet).
	void assert_updated_k(int offset) const {
		if(!updated_k(offset)) {
			std::ostringstream omsg;
			omsg	<< "IterQuantity::assert_updated_k(offset) : The interation quantity " << name() 
					<< " has not been been updated for the k"	<< std::showpos << offset
					<< " iteration yet.";
			throw QuanityNotSet(omsg.str());
		}
	}

	///
	/** Determine if the memory for the k #offset# quanity will be lost if
	  * #set_k(set_offset)# is called.
	  *
	  * This member function allows clients to know a little about the
	  * specific behavior of the subclass.  Clients can use this function
	  * to determine if it is safe to call #get_k(offset)# after #set_k(set_offset)#
	  * is called.  For example, imagine the case where you wanted to update a 
	  * vector in iteration k+1 given the elements in the k iteraiton.  For
	  * a subclass with only single storage (#info.will_loose_mem(0,+1) == true#)
	  * the following code would not work:\\
	  *
	  * #for(int i = 1; i <= n; ++i) info.set_k(+1)(i) = info.get_k(0)(i) * 2#\\
	  *
	  * For #i# == 1, #set_k(+1)# would cause a state transition and for #i# = 2
	  * #info.get_k(0)# would throw an exception.
	  *
	  * If the client knows that only single storage is needed then he should use
	  * the following code:\\
	  *
	  * #if(info.will_loose_mem(0,+1) {#\\
	  * #    info.set_k(+1) = info.get_k(0);#\\
	  * #    for(int i = 1; i <= n; ++i) info.set_k(+1)(i) = info.get_k(+1)(i) * 2#\\
	  * #}
	  * #else
	  * #    for(int i = 1; i <= n; ++i) info.set_k(+1)(i) = info.get_k(0)(i)#\\
	  *
	  * To implement the above code you would use temporary references but you get
	  * the basic idea.  The above code works for both single and multiple storage
	  * and will not result in any unnecessary copying since assingment to self
	  * should be detected.  In the above code #info.set_k(+1) = info.get_k(0);#
	  * is called to effect the state transistion.
	  *
	  * On the other hand if you need dual storage you will need a temporary
	  * copy in the event that #will_loose_mem(offset, set_offset)# returns true.
	  * For example you need dual storage for the code:\\
	  *
	  * #for(int i = 2; i <= n; ++i) info.set_k(+1)(i) = info.get_k(0)(i) * info.get_k(0)(i-1)#\\
	  *
	  * Even the above operation can be implemented without a temporary vector but you get the
	  * idea, the (i-1) quanity is modifed and is not the original for i > 2.
	  *
 	  * Preconditions:\begin{itemize}
	  * \item #updated_k(offset) == true# [throw #QuanityNotSet#]
	  * \end{itemize} 
	  */
	virtual bool will_loose_mem(int offset, int set_offset) const = 0;

	///
	/** Shift the reference point from the k to the k+1 iteration.
	  *
 	  * Postcondtions:\begin{itemize}
	  * \item #updated_k(offset)# before the call equals #updated_k(offset-1)# after return
	  * \item #&this->get_k(offset)# for #this->updated_k(offset) == true# before the call, equals
	  *			#&this->get_k(offset-1)# for #this->updated_k(offset-1) == true# after return
	  * \end{itemize} 
	  */
	virtual void next_iteration() = 0;

};	// end class IterQuantity 

}	// end namespace GeneralIterationPack

#endif	// ITER_QUANTITY_H