// ////////////////////////////////////////////////////////////////////////////////////
// rSQPState.h

#ifndef RSQP_STATE_H
#define RSQP_STATE_H

#include "ReducedSpaceSQPPackTypes.h"
#include "GeneralIterationPack/include/IterQuantityAccess.h"
#include "GeneralIterationPack/include/AlgorithmState.h"
#include "ConstrainedOptimizationPack/include/DecompositionSystem.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "LinAlgPack/include/IVector.h"
#include "LinAlgPack/include/Range1D.h"

namespace ReducedSpaceSQPPack {

// ///////////////////////////////////////////////////////////////
/** @name rSQPState iteration quantities names. */
//@{

// NLP Problem Info 
const std::string x_name				= "x";
const std::string f_name				= "f";
const std::string Gf_name				= "Gf";
const std::string Hf_name				= "Hf";
const std::string c_name				= "c";
const std::string norm_2_c_name			= "norm_2_c";
const std::string norm_inf_c_name		= "norm_inf_c";
const std::string Gc_name				= "Gc";
const std::string Hcj_name				= "Hcj";

// Constraint Gradient Null Space / Range Space Decomposition Info
const std::string Y_name				= "Y";
const std::string Z_name				= "Z";
const std::string U_name				= "U";
const std::string V_name				= "V";

// Search Direction Info
const std::string py_name				= "py";
const std::string Ypy_name				= "Ypy";
const std::string norm_2_Ypy_name		= "norm_2_Ypy";
const std::string norm_inf_Ypy_name		= "norm_inf_Ypy";
const std::string pz_name				= "pz";
const std::string Zpz_name				= "Zpz";
const std::string norm_2_Zpz_name		= "norm_2_ZpZ";
const std::string norm_inf_Zpz_name		= "norm_inf_ZpZ";
const std::string d_name				= "d";
const std::string norm_2_d_name			= "norm_2_d";
const std::string norm_inf_d_name		= "norm_inf_d";

// Reduced QP Subproblem Info
const std::string rGf_name				= "rGf";
const std::string rHL_name				= "rHL";
const std::string w_name				= "w";
const std::string zeta_name				= "zeta";
const std::string qp_grad_name			= "qp_grad";

// Global Convergence Info
const std::string alpha_name			= "alpha";
const std::string mu_name				= "mu";
const std::string Delta_name			= "Delta";
const std::string phi_name				= "phi";

// KKT Info
const std::string rGL_name				= "rGL";
const std::string norm_2_rGL_name		= "norm_2_rGL";
const std::string norm_inf_rGL_name		= "norm_inf_rGL";
const std::string lambda_name			= "lambda";
const std::string norm_2_lambda_name	= "norm_2_lambda";
const std::string norm_inf_lambda_name	= "norm_inf_lambda";
const std::string nu_name				= "nu";

//@}

///
/** Reduced space SQP state encapsulation interface.
  *
  * This in an interface to a set of data specific to a
  * reduced space SQP algorithm.  The iteration quantites
  * are abstracted within IterQuantityAccess objects.
  *
  * The interface functions use dynamic_cast<...> to cast
  * from IterQuantity to the specific types of IterQuantityAccess<...>.
  */
class rSQPState : public AlgorithmState {
public:

	// ////////////////////////////////////////////////////////////////////////////
	/** @name Public Types */
	//@{

	///
	typedef GeneralIterationPack::IterQuantityAccess<value_type>		IQA_value_type;
	///
	typedef GeneralIterationPack::IterQuantityAccess<Vector>			IQA_Vector;
	///
	typedef GeneralIterationPack::IterQuantityAccess<SpVector>			IQA_SpVector;
	///
	typedef GeneralIterationPack::IterQuantityAccess<MatrixWithOp>		IQA_MatrixWithOp;
	///
	typedef ReferenceCountingPack::ref_count_ptr<DecompositionSystem>	decomp_sys_ptr_t;

	/// Thrown if an iteration quantity is of an invalid type.
	class InvalidType : public std::logic_error
	{public: InvalidType(const std::string& what_arg) : std::logic_error(what_arg) {}};

	//@}

	/// Initializes num_basis() == 0;
	rSQPState();

	/// Virtual destructor
	virtual ~rSQPState() {}

	/// Set the interation information output level
	void iteration_info_output(EIterationInfoOutput iteration_info_output)
	{	iteration_info_output_ = iteration_info_output; }

	/// Get the interation information output level
	EIterationInfoOutput iteration_info_output() const
	{	return iteration_info_output_; }

	/// Set if the results are to be checked or not
	void check_results( bool check_results )
	{	check_results_ = check_results; }
	///
	bool check_results() const
	{	return check_results_; }

	// ////////////////////////////////////////
	// Overridden from AlgorithmState

	/// Causes fast access to iteration quantities to be lost.
	void erase_iter_quant(const std::string& iq_name);

	///
	void dump_iter_quant(std::ostream& out) const;

	/** @name «std comp» members for decomp_sys. */
	//@{

	///
	virtual void set_decomp_sys(const decomp_sys_ptr_t& decomp_sys);
	///
	virtual decomp_sys_ptr_t& get_decomp_sys();
	///
	virtual const decomp_sys_ptr_t& get_decomp_sys() const;
	///
	virtual DecompositionSystem& decomp_sys();
	///
	virtual const DecompositionSystem& decomp_sys() const;

	//@}

	/// Determine if the object is full setup and ready to use.
	virtual bool is_fully_setup() const;

	///
	/** Initializes fast access to iteration quantities after a configuration update.
	  *
	  * This function causes #*this# to access iteration quantities using id's instead
	  * by names.  This decreases lookup times from O(log(num_iter_quant()) to O(1).
	  */
	virtual void initialize_fast_access();
	
	/** @name Iteration Info */
	//@{

	/// num_basis:  Number of basis changes durring the course of the rSQP algorithm
	virtual void set_num_basis(int num_basis);
	///
	virtual int incr_num_basis();
	///
	virtual int num_basis() const;
	
	//@}

	/** @name NLP Problem Info 
	  *
	  * The problem is of the form:\\
	  *
	  * min     f(x)\\
	  * s.t.    c(x)\\
	  * .       xl < x < xu\\
	  * .       n = size(x), m = size(c)\\
	  *
	  * Gf is the gradient of f\\
	  * Gc = [Gc1, Gc2, ... , Gcm] = [Gc(indep)  Gc(dep)]\\
	  *    is the matrix of constriant gradients partitioned into independent and dependent\\
	  *    columns.
	  */
	//@{

	/// x:  The current NLP point
	virtual IQA_Vector& x();
	///
	virtual const IQA_Vector& x() const;

	/// f:  Objective function value
	virtual IQA_value_type& f();
	///
	virtual const IQA_value_type& f() const;

	/// Gf:  Gradient of the objective function sorted according to current basis selection ( n x 1 )
	virtual IQA_Vector& Gf();
	///
	virtual const IQA_Vector& Gf() const;

	/// Hf:  Hessian of the objective function ( n x n )
	virtual IQA_MatrixWithOp& Hf();
	///
	virtual const IQA_MatrixWithOp& Hf() const;

	/// c:  Vector of equality constraints ( m x 1 )
	virtual IQA_Vector& c();
	///
	virtual const IQA_Vector& c() const;

	/// norm_2_c:  ||c||2
	virtual IQA_value_type& norm_2_c();
	///
	virtual const IQA_value_type& norm_2_c() const;

	/// norm_inf_c:  ||c||infinity
	virtual IQA_value_type& norm_inf_c();
	///
	virtual const IQA_value_type& norm_inf_c() const;

	/// Gc:  Gradient of equality constraints ('c') matrix ( n x m )
	virtual IQA_MatrixWithOp& Gc();
	///
	virtual const IQA_MatrixWithOp& Gc() const;
	
	/// Hcj:  Hessain of the jth equality constraints ('c_j')( n x n )
	virtual IQA_MatrixWithOp& Hcj();
	///
	virtual const IQA_MatrixWithOp& Hcj() const;

	//@}

	/** @name Constraint Gradient Null Space / Range Space Decomposition Info. */
	//@{

	/// Y:  Range space matrix for Gc ([Y  Z] is non-singular) ( n x r )
	virtual IQA_MatrixWithOp& Y();
	///
	virtual const IQA_MatrixWithOp& Y() const;

	/// Z:  Null space matrix for Gc(con_indep)' (Gc(con_indep)' * Z) ( n x (n-r) )
	virtual IQA_MatrixWithOp& Z();
	///
	virtual const IQA_MatrixWithOp& Z() const;

	/// U:  Represents Gc(con_dep)' * Y ( (m-r) x r )
	virtual IQA_MatrixWithOp& U();
	///
	virtual const IQA_MatrixWithOp& U() const;

	/// V:  Represents Gc(con_dep)' * Z ( (m-r) x (m-r) )
	virtual IQA_MatrixWithOp& V();
	///
	virtual const IQA_MatrixWithOp& V() const;

	//@}

	/** @name Search Direction Info */
	//@{

	/// py:  Range space (dependent) QP solution component ( m x 1 )
	virtual IQA_Vector& py();
	///
	virtual const IQA_Vector& py() const;

	/// Ypy:  Range space (dependent) contribution to search direction (Ypy = Y * py) ( n x 1 )
	virtual IQA_Vector& Ypy();
	///
	virtual const IQA_Vector& Ypy() const;

	/// norm_2_Ypy:  ||Ypy||2
	virtual IQA_value_type& norm_2_Ypy();
	///
	virtual const IQA_value_type& norm_2_Ypy() const;

	/// norm_inf_Ypy:  ||Ypy||infinity
	virtual IQA_value_type& norm_inf_Ypy();
	///
	virtual const IQA_value_type& norm_inf_Ypy() const;

	/// pz:  Null space (independent) QP solution component ( (n-m) x 1 )
	virtual IQA_Vector& pz();
	///
	virtual const IQA_Vector& pz() const;

	/// Zpz:  Null space (independent) contribution to the search direction (Zpz = Z * pz) ( n x 1)
	virtual IQA_Vector& Zpz();
	///
	virtual const IQA_Vector& Zpz() const;

	/// norm_2_Zpz:  ||Zpz||2
	virtual IQA_value_type& norm_2_Zpz();
	///
	virtual const IQA_value_type& norm_2_Zpz() const;

	/// norm_inf_Zpz:  ||Zpz||infinity
	virtual IQA_value_type& norm_inf_Zpz();
	///
	virtual const IQA_value_type& norm_inf_Zpz() const;

	/// d:  Search direction (d = Zpz + Ypy) ( n x 1 )
	virtual IQA_Vector& d();
	///
	virtual const IQA_Vector& d() const;

	/// norm_2_d:  ||d||2
	virtual IQA_value_type& norm_2_d();
	///
	virtual const IQA_value_type& norm_2_d() const;

	/// norm_inf_d:  ||d||infinity
	virtual IQA_value_type& norm_inf_d();
	///
	virtual const IQA_value_type& norm_inf_d() const;

	//@}

	/** @name QP Subproblem Info */
	//@{

	/// rGf:  Reduced gradient of the objective function ( (n-r) x 1 )
	virtual IQA_Vector& rGf();
	///
	virtual const IQA_Vector& rGf() const;

	/// rHL:  Reduced Hessian of the Lagrangian function ( (n-r) x (n-r) )
	virtual IQA_MatrixWithOp& rHL();
	///
	virtual const IQA_MatrixWithOp& rHL() const;

	/// w:  QP gradient crossterm correction (Z' * HL * Y * py) ( (n-r) x 1 )
	virtual IQA_Vector& w();
	///
	virtual const IQA_Vector& w() const;

	/// zeta:  QP crossterm dampening parameter [0, 1]
	virtual IQA_value_type& zeta();
	///
	virtual const IQA_value_type& zeta() const;

	/// qp_grad:  QP gradient (qp_grad = rGf + zeta * ZtHLYpy) ( (n-m) x 1 )
	virtual IQA_Vector& qp_grad();
	///
	virtual const IQA_Vector& qp_grad() const;

	//@}

	/** @name Global Convergence Info */
	//@{

	/// alpha:  Line seach parameter
	virtual IQA_value_type& alpha();
	///
	virtual const IQA_value_type& alpha() const;

	/// mu:  Merit function penalty parameter
	virtual IQA_value_type& mu();
	///
	virtual const IQA_value_type& mu() const;

	/// Delta:  Trust region size ( n x 1 )
	virtual IQA_Vector& Delta();
	///
	virtual const IQA_Vector& Delta() const;

	/// phi:  Merrit function value
	virtual IQA_value_type& phi();
	///
	virtual const IQA_value_type& phi() const;


	//@}

	/** @name KKT Info.
	  *
	  * The Lagrangian is:\\
	  *
	  * L = f + lambda'*c + omega'*(xl - x) + gama'*(x - xu)\\
	  * L = f + lambda'*c + nu'*x + omega'*xl - gama'*xu\\
	  * where: nu = gama - omega\\
	  * 
	  * The gradient of the Lagrangian is:\\
	  *
	  * GL = Gf + Gc*lambda + nu\\
	  *
	  * The reduced gradient of the Lagrangian is:\\ 
	  *
	  * rGL = Z'*GL\\
	  * rGL = Z'*Gf + Z'*Gc*lambda + Z'*nu\\
	  * rGL = Z'*Gf + Z'*Gc(indep)*lambda(indep) + Z'*Gc(dep)*lambda(dep) + Z'*nu\\
	  *    with Gc(indep)' * Z, V = [Gc(dep)' * Z],  and rGf = Z'*Gf\\
	  * rGL = rGf + Z'*nu + V'*lambda(dep)\\
	  */
	//@{

	/// rGL:  Reduced gradient of the Lagrangian ( (n-m) x 1 )
	virtual IQA_Vector& rGL();
	///
	virtual const IQA_Vector& rGL() const;

	/// norm_2_rGL:  ||rGL||2
	virtual IQA_value_type& norm_2_rGL();
	///
	virtual const IQA_value_type& norm_2_rGL() const;

	/// norm_inf_rGL:  ||rGL||infinity
	virtual IQA_value_type& norm_inf_rGL();
	///
	virtual const IQA_value_type& norm_inf_rGL() const;

	/// lambda:  Lagrange multipliers for the equality constraints 'c' ( m x 1 )
	virtual IQA_Vector& lambda();
	///
	virtual const IQA_Vector& lambda() const;

	/// norm_2_lambda:  ||lambda||2
	virtual IQA_value_type& norm_2_lambda();
	///
	virtual const IQA_value_type& norm_2_lambda() const;

	/// norm_inf_lambda:  ||lambda||infinity
	virtual IQA_value_type& norm_inf_lambda();
	///
	virtual const IQA_value_type& norm_inf_lambda() const;

	/// nu:  Difference between Lagrange multipiers for the upper and lower bounds ( n x 1 )
	virtual IQA_SpVector& nu();
	///
	virtual const IQA_SpVector& nu() const;

	//@}

	/** @name Basis Pivot Info */
	//@{

	/// var_perm_new:  New variable permutation to form Gc' = [C  N] ( n x 1 )
	virtual IVector& var_perm_new();
	///
	virtual const IVector& var_perm_new() const;

	/// con_perm_new:  New constraint permutation to form Gc' = [C  N] ( m_full x 1 )
	virtual IVector& con_perm_new();
	///
	virtual const IVector& con_perm_new() const;

	/// var_perm_old:  Old variable permutation to form Gc' = [C  N] ( n x 1 )
	virtual IVector& var_perm_old();
	///
	virtual const IVector& var_perm_old() const;

	/// con_perm_old:  Old constraint permutation to form Gc' = [C  N] ( m_full x 1 )
	virtual IVector& con_perm_old();
	///
	virtual const IVector& con_perm_old() const;

	/// var_dep:  Range of dependent variables (columns of C, E)
	virtual Range1D& var_dep();
	///
	virtual const Range1D& var_dep() const;
	/// var_indep:  Range of independent variables (columns of N, F)
	virtual Range1D& var_indep();
	///
	virtual const Range1D& var_indep() const;
	
	/// con_indep:  Range of independent constraints (rows of C, N)
	virtual Range1D& con_indep();
	///
	virtual const Range1D& con_indep() const;
	///
	/** con_dep:  Range of dependent constraints. (rows of E, F)
	  *
	  * If there are no dependent constraints then con_dep = [m +1, m + 1]
	  */
	virtual Range1D& con_dep();
	///
	virtual const Range1D& con_dep() const;

	//@}

private:

	// ///////////////////////////////////////////////////////////////////////////
	// Private types.

	///
	enum EIQType { VALUE_TYPE, VECTOR, SP_VECTOR, MATRIX_WITH_OP }; 

	///
	enum { num_value_type_quantities = 17 };
	/// Enumeration for the value_type iteration quantities
	enum E_IterQuantities_value_type {
		Q_f,				Q_norm_2_c,			Q_norm_inf_c,			Q_norm_2_Ypy,
		Q_norm_inf_Ypy,		Q_norm_2_Zpz,		Q_norm_inf_Zpz,			Q_norm_2_d,
		Q_norm_inf_d,		Q_zeta,				Q_alpha,				Q_mu,
		Q_phi,				Q_norm_2_rGL,		Q_norm_inf_rGL,			Q_norm_2_lambda,
		Q_norm_inf_lambda
	};

	///
	enum { num_Vector_quantities = 14 };
	/// Enumeration for the Vector iteration quantities
	enum E_IterQuantities_Vector {
		Q_x,				Q_Gf,				Q_c,					Q_py,
		Q_Ypy,				Q_pz,				Q_Zpz,					Q_d,
		Q_rGf,				Q_w,				Q_qp_grad,				Q_Delta,
		Q_rGL,				Q_lambda
	};

	///
	enum { num_SpVector_quantities = 1 };
	/// Enumeration for the Vector iteration quantities
	enum E_IterQuantities_SpVector {
		Q_nu
	};

	///
	enum { num_MatrixWithOp_quantities = 8 };
	/// Enumeration for the Vector iteration quantities
	enum E_IterQuantities_MatrixWithOp {
		Q_Hf,				Q_Gc,				Q_Hcj,					Q_Y,
		Q_Z,				Q_U,				Q_V,					Q_rHL
	};

	///
	enum { num_quantites = num_value_type_quantities + num_Vector_quantities
		+ num_SpVector_quantities + num_MatrixWithOp_quantities };

	///
	enum { num_quantity_types = 4 };

	// ///////////////////////////////////////////////////////////////////////////
	// Private data members.

	// Non-iteration quantities.

	int num_basis_;
	IVector var_perm_new_, var_perm_old_, con_perm_new_, con_perm_old_;
	Range1D var_dep_, var_indep_, con_indep_, con_dep_;
	decomp_sys_ptr_t decomp_sys_;

	EIterationInfoOutput iteration_info_output_;
	// Output level.

	bool check_results_;
	// Check result flag

	bool iq_ids_initialized_;
	// Flag to keep track if the id's to the IQ objects are initialized

	static const std::string * const				iq_name_all_[num_quantites];
	// Array of pointers to the names of the IQ names.  Shared by all instances
	// of rSQPState.  These are partitioned into the different concrete types.

	AlgorithmState::iq_id_type						iq_id_all_[num_quantites];
	// Array of IQ id's that are updated when initialize_fast_access(...) is called.
	// These are partitioned into the different concrete types.

	static const std::string * const * const		iq_name_[num_quantity_types];
	// Indexed by EIQType to give array of names indexed by E_IterQuantities_«type»

	const AlgorithmState::iq_id_type *				iq_id_[num_quantity_types];
	// Indexed by EIQType to give array of IQ id's indexed by E_IterQuantities_«type»

	// ///////////////////////////////////////////////////////////////////////////
	// Private member functions.

	///
	IQA_value_type& iqa_value_type(E_IterQuantities_value_type);
	///
	const IQA_value_type& iqa_value_type(E_IterQuantities_value_type) const;

	///
	IQA_Vector& iqa_Vector(E_IterQuantities_Vector);
	///
	const IQA_Vector& iqa_Vector(E_IterQuantities_Vector) const;

	///
	IQA_SpVector& iqa_SpVector(E_IterQuantities_SpVector);
	///
	const IQA_SpVector& iqa_SpVector(E_IterQuantities_SpVector) const;

	///
	IQA_MatrixWithOp& iqa_MatrixWithOp(E_IterQuantities_MatrixWithOp);
	///
	const IQA_MatrixWithOp& iqa_MatrixWithOp(E_IterQuantities_MatrixWithOp) const;

	// not defined and not to be called
	rSQPState(const rSQPState&);
	rSQPState& operator=(const rSQPState&);

};	// end class rSQPState

}	// end namespace ReducedSpaceSQPPack

#endif	// RSQP_STATE_H