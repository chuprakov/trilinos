// //////////////////////////////////////////////////////////////////
// AlgorithmState.h

#ifndef ALGORITHM_STATE_H
#define ALGORITHM_STATE_H

#include <limits>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <iosfwd>

#include "IterQuantity.h"
#include "Misc/include/ref_count_ptr.h"

namespace GeneralIterationPack {

///
/** Abstacts a set of iteration quantities for an iterative algorithm.
  *
  * This object encapsulates a set of iteration quantity access objects.
  * The concrete types of the IterQuantity objects are expected to be subclasses
  * of IterQuantityAccess<...>.  It is therefore up the the clients to determine
  * the concrete types of these iteration quantity access objects and to use
  * dynamic_cast<...> (or static_cast<...> if you are sure) to access the
  * IterQuantityAccess<...> object and therefore the iteration quantities themselves.
  * Each iteration quantity access object (IQ) must have a unique name associated with it.
  * IQ objects are given to the state object by clients throught the set_iter_quant(...)
  * operation at with point the IQ object will be given a unique id that will never change
  * change until the IQ object is removed using erase_iter_quant(...).  Memory management
  * is performed using the ref_count_ptr<...> smart reference counting poiner class.
  * The id of any IQ object (iq_id) can be obtained from its name by calling
  * iq_id = get_iter_quant_id(iq_name).  If an IQ object with the name iq_name does not
  * exist then get_iter_quant_id(iq_name) == DOES_NOT_EXIST will be true.  The IQ objects
  * themselves can be accesed in O(log(num_iter_quant())) time using iter_quant(iq_name)
  * or in O(1) time using iter_quant(iq_id).  Therefore the access of IQ objects using iq_id
  * is an optimization for faster access and the client should never have to lookup iq_name
  * given iq_id.  The mapping only works from iq_name to iq_id, not the other way around.
  * It is garrentied that as long as erase_iter_quant(iq_id) is not called that each
  * &iter_quant(iq_id) == &iter_quant( get_iter_quant(iq_name) ) will be the same.
  * For iq_name if get_iter_quant_id(iq_name) == DOES_NOT_EXIST then iter_quant(iq_name)
  * will throw the exception DoesNotExist.
  *
  * The next_iteration(...) operation is called by the algorithm to call next_iteration()
  * on each of the IQ objects.
  *
  * The dump_iter_quant(out) operation prints out a list of all of the IQ objects of thier
  * iq_name, iq_name and concrete type.
  *
  * The default copy constructor, and assignment operator functions
  * are allowed since they have the proper sematics.
  */
class AlgorithmState {
public:

	/** @name Public types */
	//@{

	///
	typedef size_t													iq_id_type;
	///
	typedef ReferenceCountingPack::ref_count_ptr<IterQuantity>		IQ_ptr;
	///
	enum { DOES_NOT_EXIST = 1000 };	// should not ever be this many insertions.

	/// Thrown if name or id does not exist
	class DoesNotExist : public std::logic_error
	{public: DoesNotExist(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Thrown if name already exists
	class AlreadyExists : public std::logic_error
	{public: AlreadyExists(const std::string& what_arg) : std::logic_error(what_arg) {}};

	//@}

	///
	/** Construct with an initial guess for the number of iteration quantities.
	  *
	  * The iteration counter k is default constructed to zero.
	  */
	explicit AlgorithmState(size_t reserve_size = 0);

	/** @name iteration counter */
	//@{

	///
	void k(int k);
	///
	int k() const;
	///
	int incr_k();

	//@}

	/** @name iteration quantity encapsulation object setup */
	//@{
	
	/// Return the number of iteration quantities.
	virtual size_t num_iter_quant() const;

	///
	/** Inserts the iteration quantity through a ref_count_ptr<...> object.
	  *
	  * Time = O(log(num_iter_quant)), Space = O(1).
	  *
	  * If an iteration quantity already exists with the name #iq_name# then
	  * a #AlreadyExists# exception will be thrown.  Otherwise the function
	  * will return the iq_id assigned to the inserted interation quantity.
	  *
	  * Preconditions: \begin{itemize}
	  * \item #get_iter_quant_id(iq_name) == DOES_NOT_EXIST# (throw #AlreadyExists#)
	  * \end{itemize}
	  */
	virtual iq_id_type set_iter_quant(const std::string& iq_name, const IQ_ptr& iq);

	///
	/** Removes the iteration quantity with name iq_name.
	  *
	  * Time = O(log(num_iter_quant)), Space = O(1).
	  *
	  * If #get_iter_quant(iq_name).count() == 1# then the IterQuantity object
	  * pointed to will be deleted.  Subsequently, the iq_id returned from
	  * #set_iter_quant(...)# when #iq_name# was set is no longer valid.
	  *
	  * Preconditions: \begin{itemize}
	  * \item #get_iter_quant_id(iq_name) != DOES_NOT_EXIST# (throw #DoesNotExist#)
	  * \end{itemize}
	  */
	virtual void erase_iter_quant(const std::string& iq_name);

	///
	/** Return the iteration quantity id (iq_id) for the iteration quantity.
	  *
	  * If an iteration quantity with the name #iq_name# does not exist, then
	  * the value DOES_NOT_EXIST is returned.
	  *
	  * Time = O(log(num_iter_quant)), Space = O(1).
	  */
	virtual iq_id_type get_iter_quant_id(const std::string& iq_name) const;

	///
	/** Returns the ref_count_ptr<...> for the iteration quantiy with iq_id
	  *
	  * Time = O(1), Space = O(1).
	  */
	virtual IQ_ptr& get_iter_quant(iq_id_type iq_id);

	///
	virtual const IQ_ptr& get_iter_quant(iq_id_type iq_id) const;

	//@}

	///
	/** Iteration quantity encapsulation object access with via iq_name.
	  *
	  * Time = O(log(num_iter_quant())), Space = O(1).
	  *
	  * Preconditions: \begin{itemize}
	  * \item #get_iter_quant_id(iq_name) != DOES_NOT_EXIST# (throw #DoesNotExist#)
	  * \end{itemize}
	  */
	virtual IterQuantity& iter_quant(const std::string& iq_name );
	///
	virtual const IterQuantity& iter_quant(const std::string& iq_name ) const;
	///
	/** Iteration quantity encapsulation object access via iq_id.
	  *
	  * Time = O(1), Space = O(1).
	  *
	  * If the IQ object with iq_id does not exist then a #std::out_of_range#
	  * or #std::logic_error# will be thrown.
	  */
	virtual IterQuantity& iter_quant(iq_id_type iq_id);
	///
	virtual const IterQuantity& iter_quant(iq_id_type iq_id) const;

	///
	/** iteration quantity forwarding.
	  *
	  */
	virtual void next_iteration(bool incr_k = true);

	///
	/** iteration quantity information dumping.
	  *
	  * This function outputs a list with columns:
	  *
	  * iq_name		iq_id		concrete type
	  *
	  * for each iteration quantity object.
	  */
	virtual void dump_iter_quant(std::ostream& out) const;

private:
	// ///////////////////////////////////////////////////////////
	// Private types

	typedef std::vector<IQ_ptr>					iq_t;
	typedef std::map<std::string,iq_id_type>	iq_name_to_id_t;

	// ///////////////////////////////////////////////////////////
	// Private data members

	int k_;		// Iteration counter.

	iq_t					iq_;
	// Array of ref_count_ptr objects that point to set iteration quantities.
	// The index into this array is the iq_id for an IQ object.  This array
	// is filled sequantially from the beginning using push_back(...).
	// When erase_iter_quant(...) is called the iq_[iq_id] is set to null which
	// reduces the reference count of the IQ object (possible deleing it if
	// there are no other references).  Then if the user tries to access and
	// IQ object with this abandonded iq_id, the dereferencing operator for
	// ref_count_ptr<...> will throw an exception.

	iq_name_to_id_t			iq_name_to_id_;
	// Mapping of an IQ name to its id.

	// ///////////////////////////////////////////////////////////
	// Private member functions
	
	///
	iq_name_to_id_t::iterator find_and_assert(const std::string& iq_name);
	///
	iq_name_to_id_t::const_iterator find_and_assert(const std::string& iq_name) const;

};	// end class AlgorithmState

// /////////////////////////////////////////////////////////////////////////////////
// Inline member definitions for AlgorithmState

inline
AlgorithmState::AlgorithmState(size_t reserve_size)
	: k_(0)
{	iq_.reserve(reserve_size); }

// iteration counter

inline
void AlgorithmState::k(int k)
{	k_ = k; }

inline
int AlgorithmState::k() const
{	return k_; }

inline
int AlgorithmState::incr_k()
{	return ++k_; }

// 

inline
size_t AlgorithmState::num_iter_quant() const {
	return iq_name_to_id_.size();
}

inline
AlgorithmState::iq_id_type AlgorithmState::get_iter_quant_id(
	const std::string& iq_name) const
{
	const iq_name_to_id_t::const_iterator itr = iq_name_to_id_.find(iq_name);
	return itr == iq_name_to_id_.end() ? DOES_NOT_EXIST : itr->second;
}

inline
AlgorithmState::IQ_ptr& AlgorithmState::get_iter_quant(iq_id_type iq_id) {
	return iq_.at(iq_id);
}

inline
const AlgorithmState::IQ_ptr& AlgorithmState::get_iter_quant(
	iq_id_type iq_id) const
{
	return iq_.at(iq_id);
}

inline
IterQuantity& AlgorithmState::iter_quant(const std::string& iq_name ) {
	return *iq_[find_and_assert(iq_name)->second];
}

inline
const IterQuantity& AlgorithmState::iter_quant(const std::string& iq_name ) const {
	return *iq_[find_and_assert(iq_name)->second];
}

inline
IterQuantity& AlgorithmState::iter_quant(iq_id_type iq_id) {
	return *iq_.at(iq_id);
}

inline
const IterQuantity& AlgorithmState::iter_quant(iq_id_type iq_id) const {
	return *iq_.at(iq_id);
}

}	// end namespace GeneralIterationPack 

#endif // ALGORITHM_STATE_H