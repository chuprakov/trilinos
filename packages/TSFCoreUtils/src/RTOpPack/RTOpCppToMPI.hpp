// //////////////////////////////////////////////
// RTOpCppToMPI.hpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#ifndef RTOP_CPP_TO_MPI_H
#define RTOP_CPP_TO_MPI_H

#include "RTOpCpp.hpp"

namespace RTOpPack {

/** \defgroup RTOpCppToMPI_grp Helper fuctions for applying reduction/transformation operations with MPI.
  */
//@{

///
/** Apply a reduction operation over a set of local sub-vectors using MPI.
 *
 *
 *	@param	comm
 *				[in] MPI communicator
 *	@param	op	[in] reduction/transformation operator object
 *	@param	root_rank
 *				[in] See below
 *	@param	num_cols
 *				[in] The number of columns for each vector
 *	@param	num_vecs
 *				[in] See <tt>RTOpPack::RTOp::apply_op()</tt>
 *	@param	sub_vecs
 *				[in] Array (size <tt>num_vecs*num_cols</tt>)
 *				of nonmutable subvectors.  The vectors for each column <tt>kc</tt>
 *              begin at <tt>sub_vecs+kc*num_cols</tt> where <tt>kc=0...num_cols-1</tt>.
 *              Can be <tt>NULL</tt> if <tt>num_vecs==0</tt> or if there are no local
 *              vector elements that participate in the RTOp operation.
 *	@param	num_targ_vecs
 *				[in] See <tt>%RTOpPack::RTOp::apply_op()</tt>
 *	@param	sub_targ_vecs
 *				[in] Array (size <tt>num_targ_vecs*num_cols</tt>)
 *				of mutable subvectors.  The vectors for each column <tt>kc</tt>
 *              begin at <tt>sub_targ_vecs+kc*num_cols</tt> where <tt>kc=0...num_cols-1</tt>
 *              Can be <tt>NULL</tt> if <tt>num_targ_vecs==0</tt> or if there are no local
 *              vector elements that participate in the RTOp operation.
 *	@param	reduct_objs
 *				[in/out] Array (size <tt>num_cols</tt>) See below.
 *				If <tt>reduct_objs != NULL</tt>
 *				then on output, <tt>reduct_objs[i]</tt> will contain the reduction target
 *				over all the processes along with the reduction on input if
 *              <tt>reduct_objs[i]</tt> has already been through one or more reductions
 *              already and not reinitialized.
 *
 * This function encapsulates a lot of the details of using MPI to perform
 * reduction operations.  For this function to work properly, each
 * processor must contain a single sub-vector for each vector
 * argument and every process (in the communicator) must call this function.
 * Each process's sub-vector for each vector argument
 * may be a different size and even dense on some processes
 * and sparse in others.  It is also allowed for a process to have
 * empty sub-vectors, but still needs to participate in the collective
 * operation and may want the value of the reduction object returned.
 * The main responsibly of the client is to setup the <tt>sub_vecs[]</tt>
 * and <tt>sub_targ_vecs[]</tt> arguments and to deside what type of
 * reduction operation to perform (\c MPI_Allreduce() or MPI_Reduce()).
 *
 * Let <tt>rank</tt> be the value returned from <tt>MPI_Comm_rank(comm,&rank)</tt> in
 * this process.  The expected arguments to this function depend
 * on the argument <tt>root_rank</tt> and <tt>rank</tt> as follows:
 *
 * <ul>
 *	<li> <tt>root_rank < 0</tt> :  In this case we are performing a
 *		<tt>MPI_Allreduce(...)</tt> reduction operation over all of
 *		the processes with the results collected in all of
 *		the processes.  In this case <tt>reduct_obj!=RTOp_REDUCT_OBJ_NULL</tt>
 *		must be true in all of the processes.
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
	,const int num_cols
	,const int num_vecs, const RTOpPack::SubVector sub_vecs[]
	,const int num_targ_vecs, const RTOpPack::MutableSubVector targ_sub_vecs[]
	,RTOp_ReductTarget reduct_objs[]
	);

//@}

} // end namespace RTOpPack

#endif // RTOP_CPP_TO_MPI_H
