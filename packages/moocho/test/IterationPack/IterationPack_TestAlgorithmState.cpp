// /////////////////////////////////////////////////////////////
// TestAlgorithmState.cpp

// ToDo: RAB: 7/1/99: remove a dependancy of on LinAlgPack and
// SparseLinAlgPack from the testing software you dork.

#pragma warning(disable : 4786)	// too long class name for debugger warning

#include <ostream>
#include <iomanip>
#include <typeinfo>

#include "../test/TestGeneralIterationPack.h"
#include "../include/AlgorithmState.h"
#include "../include/IterQuantityAccessContinuous.h"
#include "../include/IterQuantityAccessDerivedToBase.h"
#include "Misc/include/update_success.h"

// explicit instantiation for testing compilation only
//template ReferenceCountingPack::ref_count_ptr<double>;
//class B {};
//class D : public B {};
//template GeneralIterationPack::IterQuantityAccessDerivedToBase<B,D>;

namespace {

// Create class hierarchy for derived to base conversion.

class B {
public:
	virtual B& operator=(const B&) = 0;
	virtual void alpha(double) = 0;
	virtual double alpha() const = 0; 
};	// end class B

class D : public B {
public:
	B& operator=(const B& b) {
		const D* d = dynamic_cast<const D*>(&b);
		if(!d)
			throw std::logic_error("D::operator=(...) : b is not of type D");
		alpha_ = d->alpha_;
		return *this;
	}
	void alpha(double alpha) {
		alpha_ = alpha;
	}
	double alpha() const {
		return alpha_;	
	}
private:
	double alpha_; 
};	// end class D

}	// end namespace

bool GeneralIterationPack::TestingPack::TestAlgorithmState(std::ostream* out) {

	using std::endl;
	using std::setw;
	
	using TestingHelperPack::update_success;

	try {

	bool success = true, result;
//	const int w = 15;
	if(out) *out << std::boolalpha;
		
	if(out)
		*out<< "\n\n******************************\n"
			<< "*** Testing AlgorithmState ***\n"
			<< "******************************\n"
			<< "\nWarning: this interface is weakly typed and mistakes"
			<< "can be very bad\n";

	typedef double													alpha_k_t;
	typedef std::vector<double>										x_k_t;
	typedef B														V_k_t;
	typedef IterQuantityAccessContinuous<alpha_k_t>					alpha_t;
	typedef IterQuantityAccessContinuous<x_k_t>						x_t;
	typedef IterQuantityAccessDerivedToBase<V_k_t,D>				V_t;

	if(out) *out << "\n*** Create state object ...\n";
	AlgorithmState state(4);

	// Set three types of iteration quantity access objects.

	if(out) *out << "\n*** Set three types of iteration quantity access objects.\n";

	typedef AlgorithmState::IQ_ptr		IQ_ptr;
	typedef V_t::IQA_ptr				IQA_ptr;

	if(out) *out << "set IterQuantityAccessContinuous<double>(2,\"alpha\")\n";
	state.set_iter_quant( "alpha", IQ_ptr(new alpha_t(2,"alpha")) );

	if(out) *out << "set IterQuantityAccessContinuous<std::vector<double> >(2,\"x\")\n";
	state.set_iter_quant( "x", IQ_ptr(new x_t(2,"x")) );

	if(out) *out << "set IterQuantityAccessDerivedToBase<B,D>(1,\"V\")\n";
	state.set_iter_quant( "V"
		, IQ_ptr(
		      new V_t( IQA_ptr( new IterQuantityAccessContinuous<D>(1,"V") ) )
		    )
	  );

	// Try to add the same name again.

	if(out) *out << "\nTry to add \"x\" (should throw execption)\n";
	try {
		state.set_iter_quant( "x",IQ_ptr(0) );
		success = false;
	}
	catch(const AlgorithmState::AlreadyExists& expt) {
		if(out) *out << "Caught a AlgorithmState::AlreadyExists execption : " << expt.what() << endl;
	}

	// dump the iteration quantity names, ids, and concrete types.

	if(out) {
		*out << "\n*** dump iteration quantitys\n";
		state.dump_iter_quant(*out);
	}

	update_success( result = state.k() == 0, &success );
	if(out) *out << "\nstate.k() == 0 : " << result << endl; 

	// Set the iteration quantities for the kth iteration and check them

	if(out) *out << "\n*** Set iteration quantities for the kth iteration\n";

	if(out) *out << "\nSet alpha_k(0) = 5.0 then get alpha_k(0) == 5.0 : ";
	alpha_t *alpha = dynamic_cast<alpha_t*>( &state.iter_quant("alpha") );
	alpha->set_k(0) = 5.0; 
	alpha = dynamic_cast<alpha_t*>( &state.iter_quant("alpha") );
	update_success( result = alpha->get_k(0) == 5.0, &success );
	if(out) *out << result << endl;

	if(out) *out << "\nSet x_k(0)[0] = 2.0 then get x_k(0)[0] == 2.0 : ";
	x_t *x = dynamic_cast<x_t*>( &state.iter_quant("x") );
	x->set_k(0).resize(1);
	x->set_k(0)[0] = 2.0;
	x = dynamic_cast<x_t*>( &state.iter_quant("x") );
	update_success( result = x->get_k(0)[0] == 2.0, &success );
	if(out) *out << result << endl;

	if(out) *out << "\nSet V_k(0).alpha() = 3.0 then get V_k(0).alpha() == 3.0 : ";
	V_t *V = dynamic_cast<V_t*>( &state.iter_quant("V") );
	V->set_k(0).alpha(3.0);
	V = dynamic_cast<V_t*>( &state.iter_quant("V") );
	update_success( result = V->get_k(0).alpha() == 3.0, &success );
	if(out) *out << result << endl;

	// Use an id to get at an iteration quantity.

	if(out) *out << "\n*** Use an id to get at an iteration quantity\n";

	if(out) *out << "\nid for \"x\" is : ";
	AlgorithmState::iq_id_type x_id = state.get_iter_quant_id("x");
	if(out) *out << x_id << endl;

	if(out) *out << "\nAccess \"x\" using id and x_k(0)[0] == 2.0 : ";
	x = dynamic_cast<x_t*>( &state.iter_quant(x_id) );
	update_success( result = x->get_k(0)[0] == 2.0, &success );
	if(out) *out << result << endl;

	// use a nonexistant name

	if(out) *out << "\n*** Use a nonexistant name or id to access iteration quantity\n";

	if(out) *out << "id for \"X\" is DOES_NOT_EXIST : ";
	update_success( result = state.get_iter_quant_id("X") == AlgorithmState::DOES_NOT_EXIST
		, & success );
	if(out) *out << result << endl;

	if(out) *out << "\nAccess nonexistant iteration quantity by name \"X\" throws a AlgorithmState::DoesNotExist exception : ";
	try {
		state.iter_quant("X");
		success = false;
		if(out) *out << false << endl;
	}
	catch(const AlgorithmState::DoesNotExist& expt) {
		if(out) *out << true << endl;
	}

	// use a nonexistant id

	if(out) *out << "\nUse of a nonexistant id = 100 throws a AlgorithmState::DoesNotExist exception : ";
	try {
		state.iter_quant(100);
		success = false;
		if(out) *out << false << endl;
	}
	catch(const AlgorithmState::DoesNotExist& expt) {
		if(out) *out << true << endl;
	}

	// update the iteration

	if(out) *out << "\n*** Update iteration quantities k+1 = k then check\n";

	if(out) *out << "alpha_k(+1) = alpha_k(0)...\n";
	alpha = dynamic_cast<alpha_t*>( &state.iter_quant("alpha") );
	{
		alpha_k_t &alpha_k = alpha->get_k(0);
		alpha->set_k(+1) = alpha_k;
	}

	if(out) *out << "x_k(+1) = x_k(0)...\n";
	x = dynamic_cast<x_t*>( &state.iter_quant("x") );
	{
		x_k_t &x_k = x->get_k(0);
		x->set_k(+1) = x_k;
	}

	if(out) *out << "V_k(+1) = V_k(0)...\n";
	V = dynamic_cast<V_t*>( &state.iter_quant("V") );
	{
		V_k_t &V_k = V->get_k(0);
		V->set_k(+1) = V_k;
	}

	if(out) *out << "shift reference from k to k+1...\n";
	state.next_iteration();

	if(out) *out << "\nalpha_k(-1) == 5.0 : ";
	alpha = dynamic_cast<alpha_t*>( &state.iter_quant("alpha") );
	update_success( result = alpha->get_k(-1) == 5.0, &success );
	if(out) *out << result << endl;
	if(out) *out << "alpha_k(0) == 5.0 : ";
	update_success( result = alpha->get_k(0) == 5.0, &success );
	if(out) *out << result << endl;

	if(out) *out << "\nx_k(-1)[0] == 2.0 : ";
	x = dynamic_cast<x_t*>( &state.iter_quant("x") );
	update_success( result = x->get_k(-1)[0] == 2.0, &success );
	if(out) *out << result << endl;
	if(out) *out << "x_k(0)[0] == 2.0 : ";
	update_success( result = x->get_k(0)[0] == 2.0, &success );
	if(out) *out << result << endl;

	if(out) *out << "\nV_k(0).alpha() == 3.0 : ";
	V = dynamic_cast<V_t*>( &state.iter_quant("V") );
	update_success( result = V->get_k(0).alpha() == 3.0, &success );
	if(out) *out << result << endl;

	// erase an iteration quantity then try to access it

	if(out) *out << "\n*** Erasing iteration quantity \"x\" then trying to access it throws"
					" a AlgorithmState::DoesNotExist exception : ";
	state.erase_iter_quant("x");
	try {
		state.iter_quant(x_id);
		if(out) *out << false << endl;
		update_success( false, &success );
	}
	catch(const AlgorithmState::DoesNotExist& expt) {
		if(out) *out << true << endl;
	}

	// final printout.

	if(out) {
		if(success) {
			*out << "\n*** Congradulations, all of the tests for AlgorihtmState"
					" returned the expected results\n";
		}
		else {
			*out << "\n*** Oops, at least one of the above tests for AlgorihtmState"
					" did not return the expected results ***\n";
		}
	}

	return success;

	} // end try
	catch(const std::exception& excpt) {
		if(out) *out << "\nCaught a std::exception: " << typeid(excpt).name() << " : " <<  excpt.what() << endl;
	}
	catch(...) {
		if(out) *out << "\nCaught an unknown exception\n";
	}

	return false;
}
