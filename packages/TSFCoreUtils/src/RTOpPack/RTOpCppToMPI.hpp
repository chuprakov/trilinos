// @HEADER
// ***********************************************************************
// 
//      TSFCoreUtils: Trilinos Solver Framework Utilities Package 
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef RTOP_CPP_TO_MPI_H
#define RTOP_CPP_TO_MPI_H

#include "RTOpCpp.hpp"

namespace RTOpPack {

/** \defgroup RTOpCppToMPI_grp Helper fuctions for applying reduction/transformation operations with MPI.
  */
//@{

///
/** Apply a reduction/transformation operation over a set of local
 *	sub-vectors using MPI.
 *
 *	@param	comm
 *				[in] MPI communicator.
 *	@param	op	[in] Reduction/transformation operator object.
 *	@param	root_rank
 *				[in] See below.
 *	@param	num_vecs
 *				[in] See <tt>RTOpPack::RTOp::apply_op()</tt>
 *	@param	sub_vecs
 *				[in] Array (size <tt>num_vecs</tt>) of nonmutable subvectors.
 *              Can be <tt>NULL</tt> if <tt>num_vecs==0</tt> or if there are no local
 *              vector elements that participate in the RTOp operation.
 *	@param	num_targ_vecs
 *				[in] See <tt>%RTOpPack::RTOp::apply_op()</tt>
 *	@param	sub_targ_vecs
 *				[in] Array (size <tt>num_targ_vecs*num_cols</tt>) of mutable subvectors.
 *              Can be <tt>NULL</tt> if <tt>num_targ_vecs==0</tt> or if there are no local
 *              vector elements that participate in the RTOp operation.
 *	@param	reduct_obj
 *				[in/out] Reduction object.  If there is no reduction then 
 *              <tt>reduct_obj</tt> must be <tt>RTOp_REDUCT_OBJ_NULL</tt>.
 *
 * This function encapsulates a lot of the details of using MPI to
 * perform reduction operations.  For this function to work properly,
 * each processor must contain a single sub-vector for each vector
 * argument and every process (in the communicator) must call this
 * function.  Each process's sub-vector for each vector argument may
 * be a different size.  It is also allowed for a process to have
 * empty sub-vectors, but each process still needs to participate in
 * the collective operation (and may want the value of the reduction
 * object returned anyway).  The main responsibly of the client is to
 * setup the <tt>sub_vecs[]</tt> and <tt>sub_targ_vecs[]</tt>
 * arguments and to decide what type of reduction operation to perform
 * (<tt>MPI_Allreduce()</tt> or <tt>MPI_Reduce()</tt>).
 *
 * Let <tt>rank</tt> be the value returned from
 * <tt>MPI_Comm_rank(comm,&rank)</tt> in this process.  The expected
 * arguments to this function depend on the argument
 * <tt>root_rank</tt> and <tt>rank</tt> as follows:
 *
 * <ul>
 *	<li> <tt>root_rank < 0</tt> :  In this case we are performing a
 *		<tt>MPI_Allreduce(...)</tt> reduction operation over all of
 *		the processes with the results collected in all of
 *		the processes.  In this case <tt>reduct_obj!=RTOp_REDUCT_OBJ_NULL</tt>
 *		must be true in all of the processes (if a reduction operation
 *      is being performed, otherwise <tt>reduct_obj</tt> must be <tt>RTOp_REDUCT_OBJ_NULL</tt>
 *      of course).
 *		<ul>
 *		<li> <tt>root_rank < 0</tt> in all processes
 *		<li> <tt>reduct_obj != RTOp_REDUCT_OBJ_NULL</tt> in all processes
 *           if a reduction is being performed.
 *		</ul>
 *	<li> <tt>root_rank >= 0</tt> :  In this case we are performing a
 *		<tt>MPI_Reduce(...)</tt> reduction operation over all of the processes
 *		with the result only collected in the process with rank
 *		<tt>root_rank</tt>.  Here all of the processes must pass the same
 *		value for <tt>root_rank</tt>.  The reduction target object is only
 *		passed in for the processes with <tt>rank == root_rank</tt>.
 *		<ul>
 *		<li> <tt>root_rank</tt> same in all of the processes
 *		<li> If <tt>rank == root_rank</tt> then <tt>reduct_obj != RTOp_REDUCT_OBJ_NULL</tt>
 *		<li> If <tt>rank != root_rank</tt> then <tt>reduct_obj == RTOp_REDUCT_OBJ_NULL</tt>
 *		</ul>
 *	</ul>
 *
 * Note that if a reduction operation is not performed then no synchronization of the the processes
 * are performed.  This is the fastest behavior and should be fine for a single thread
 * but may cause problems with multiple threads.
 */
void  MPI_apply_op(
	MPI_Comm comm, const RTOp& op, int root_rank
	,const int num_vecs, const RTOpPack::SubVector sub_vecs[]
	,const int num_targ_vecs, const RTOpPack::MutableSubVector targ_sub_vecs[]
	,RTOp_ReductTarget reduct_obj
	);

///
/** Apply a reduction/transformation operation over a set of local
 *	sub-multi-vectors using MPI.
 *
 * ToDo: Finish documentation!
 */
void MPI_apply_op(
	MPI_Comm comm, const RTOp& op, int root_rank
	,const int num_cols
	,const int num_multi_vecs, const RTOpPack::SubMultiVector sub_multi_vecs[]
	,const int num_targ_multi_vecs, const RTOpPack::MutableSubMultiVector targ_sub_multi_vecs[]
	,RTOp_ReductTarget reduct_objs[]
	);

//@}

} // end namespace RTOpPack

#endif // RTOP_CPP_TO_MPI_H
