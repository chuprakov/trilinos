#ifdef EPETRA_ANSI_CPP
#include <typeinfo>
#include <sstream>
#endif

#include "TSFDefs.h"
#include "PetraVector.h"
#include "PetraVectorSpace.h"

#include "TSFUtils.h"
#include "TSFOut.h"
#include "TSFTimeMonitor.h"



using namespace TSF;


PetraVector::PetraVector(const TSFVectorSpace& space,
												 const TSFSmartPtr<Epetra_Vector>& localValues,
												 const TSFSmartPtr<Epetra_Vector>& allValues,
												 bool ghostValuesAreValid)
	: TSFMPIVector(space),
		allValues_(allValues),
		localValues_(localValues),
		ghostValuesAreValid_(ghostValuesAreValid)
{}

PetraVector::PetraVector(const TSFVectorSpace& space,
												 const TSFSmartPtr<Epetra_Vector>& localValues)
	: TSFMPIVector(space),
		allValues_(localValues),
		localValues_(localValues),
		ghostValuesAreValid_(false)
{}

#ifdef HAVE_RTOP

void PetraVector::apply_reduction(
	const RTOp_RTOp &op, int num_vecs, const TSFVectorBase* vecs_in[]
	,int num_targ_vecs, TSFVectorBase* targ_vecs_in[], RTOp_ReductTarget reduct_obj
	,const RTOp_index_type first_ele, const RTOp_index_type sub_dim, const RTOp_index_type global_offset
	) const
{
	// Convert from TSFVectorBase arrays to Epetra_MultiVector arrays
	typedef const Epetra_MultiVector*   const_vec_ptr_t;
	typedef Epetra_MultiVector*         vec_ptr_t;
	const_vec_ptr_t      *vecs      = NULL;
	vec_ptr_t            *targ_vecs = NULL;
	if(num_vecs) {
		vecs = new const_vec_ptr_t[num_vecs];
		for( int k = 0; k < num_vecs; ++k )
			vecs[k] = &(*dynamic_cast<const PetraVector*>(vecs_in[k])->localValues_);
	}
	if(num_targ_vecs) {
		targ_vecs = new vec_ptr_t[num_targ_vecs];
		for( int k = 0; k < num_targ_vecs; ++k ) {
			PetraVector *petra_v = dynamic_cast<PetraVector*>(targ_vecs_in[k]);
			targ_vecs[k] = &(*petra_v->localValues_);
			petra_v->invalidateGhostValues(); // These vector will be changes!
		}
	}
	// Call the implementation to invoke the operator
	localValues_->apply_reduction(
		op,num_vecs,vecs,num_targ_vecs,targ_vecs
		,&reduct_obj // First element in an array of one RTOp_ReductTarget objects
		,first_ele,sub_dim,global_offset
		);
	// Delete the arrays
	if(vecs)      delete [] vecs;
	if(targ_vecs) delete [] targ_vecs;
}

void PetraVector::apply_transformation(
	const RTOp_RTOp &op, int num_vecs, const TSFVectorBase* vecs_in[]
	,int num_targ_vecs, TSFVectorBase* targ_vecs_in[], RTOp_ReductTarget reduct_obj
	,const RTOp_index_type first_ele, const RTOp_index_type sub_dim, const RTOp_index_type global_offset
	)
{
	// Convert from TSFVectorBase arrays to Epetra_MultiVector arrays
	typedef const Epetra_MultiVector*   const_vec_ptr_t;
	typedef Epetra_MultiVector*         vec_ptr_t;
	const_vec_ptr_t      *vecs      = NULL;
	vec_ptr_t            *targ_vecs = NULL;
	if(num_vecs) {
		vecs = new const_vec_ptr_t[num_vecs];
		for( int k = 0; k < num_vecs; ++k )
			vecs[k] = &(*dynamic_cast<const PetraVector*>(vecs_in[k])->localValues_);
	}
	if(num_targ_vecs) {
		targ_vecs = new vec_ptr_t[num_targ_vecs];
		for( int k = 0; k < num_targ_vecs; ++k ) {
			PetraVector *petra_v = dynamic_cast<PetraVector*>(targ_vecs_in[k]);
			targ_vecs[k] = &(*petra_v->localValues_);
			petra_v->invalidateGhostValues(); // These vector will be changes!
		}
	}
	// Call the implementation to invoke the operator
	localValues_->apply_transforamtion(
		op,num_vecs,vecs,num_targ_vecs,targ_vecs
		,&reduct_obj // First element in an array of one RTOp_ReductTarget objects
		,first_ele,sub_dim,global_offset
		);
	this->invalidateGhostValues();
	// Delete the arrays
	if(vecs)      delete [] vecs;
	if(targ_vecs) delete [] targ_vecs;
}

#endif

TSFReal& PetraVector::setElement(int globalIndex) 
{
	return allValues_->operator[](getLocalIndex(globalIndex));
}

const TSFReal& PetraVector::getElement(int globalIndex) const 
{
	return allValues_->operator[](getLocalIndex(globalIndex));
}

void PetraVector::setElements(int n, const int* globalIndices,
															const TSFReal* values) 
{
	for (int i=0; i<n; i++)
		{
			setElement(globalIndices[i]) = values[i];
		}
}

void PetraVector::getElements(int n, const int* globalIndices,
															TSFReal* values) const 
{
	for (int i=0; i<n; i++)
		{
			values[i] = getElement(globalIndices[i]);
		}
}

void PetraVector::addToElements(int n, const int* globalIndices,
																const TSFReal* values) 
{
	for (int i=0; i<n; i++)
		{
			setElement(globalIndices[i]) += values[i];
		}
}

const Epetra_Vector& PetraVector::values() const
{
	if (ghostValuesAreValid_) return *allValues_;
	else return *localValues_;
}

Epetra_Vector& PetraVector::values() 
{
	if (ghostValuesAreValid_) return *allValues_;
	else return *localValues_;
}

const Epetra_Vector& PetraVector::values(bool otherValuesAreValid) const
{
	if (ghostValuesAreValid_ && otherValuesAreValid) return *allValues_;
	else return *localValues_;
}

Epetra_Vector& PetraVector::values(bool otherValuesAreValid) 
{
	if (ghostValuesAreValid_ && otherValuesAreValid) return *allValues_;
	else return *localValues_;
}

const Epetra_Import& PetraVector::importer() const 
{
	return *(PetraVectorSpace::getImporter(space()));
}

const Epetra_Map& PetraVector::getLocalMap() const 
{
	return *(PetraVectorSpace::getLocalMap(space()));
}

const Epetra_Map& PetraVector::getMap() const 
{
	return *(PetraVectorSpace::getMap(space()));
}

TSFVectorBase* PetraVector::deepCopy() const 
{
	if (&(*allValues_) == &(*localValues_))
		{
			Epetra_Vector* v = new Epetra_Vector(*localValues_);
			return new PetraVector(space(), v);
		}
	else
		{
			Epetra_Vector* allValues = new Epetra_Vector(*allValues_);
			double* values;
			allValues->ExtractView(&values);
			Epetra_Vector* localValues 
				= new Epetra_Vector(View, getLocalMap(), values);
			return new PetraVector(space(), localValues, allValues, 
														 ghostValuesAreValid_);
		}
}


ostream& PetraVector::print(ostream& os) const 
{
	os << endl << *allValues_;
	return os;
}

int PetraVector::getLocalIndex(int globalIndex) const 
{
	return allValues_->Map().LID(globalIndex);
}

const Epetra_Vector& PetraVector::getLocalValues(const TSFVector& v)
{
	return getConcrete(v).localValues();
}


Epetra_Vector& PetraVector::getLocalValues(TSFVector& v)
{
	return getConcrete(v).localValues();
}

const PetraVector& PetraVector::getConcrete(const TSFVector& x)
{
	const PetraVector* v = dynamic_cast<const PetraVector*>(x.ptr());
	if (v==0) {
#ifdef EPETRA_ANSI_CPP
	TSFOStringStream omsg;
	omsg << "PetraVector::getConcrete(x): Error, bad cast!"
		 << "  The concrete object of type \'" << typeid(*x.ptr()).name()
		 << "\' does not support the type \'PetraVector\'."
		;
	TSFError::raise(omsg.str().c_str());
#else
	TSFError::raise("PetraVector::getConcrete bad cast");
#endif
	}
	return *v;
}

PetraVector& PetraVector::getConcrete(TSFVector& x)
{
	PetraVector* v = dynamic_cast<PetraVector*>(x.ptr());
	if (v==0) TSFError::raise("PetraVector::getConcrete bad cast");
	return *v;
}

void PetraVector::synchronizeGhostValues() const 
{
	if (ghostValuesAreValid_ || PetraVectorSpace::getImporter(space()).get()==0) return;
	Epetra_Vector& v = const_cast<Epetra_Vector&>(*allValues_);
	const Epetra_Vector& lv = localValues();
	const Epetra_Import& i = importer();
	{
		TSFTimeMonitor t(commTimer());
		v.Import(lv, i, Insert);
	}
	ghostValuesAreValid_ = true;
}

























