/* //////////////////////////////////////////////
// RTOp_mpi.c
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
//
// Selected hollow MPI function definitions for a sinlge process
// implementation.
*/

#include <assert.h>

#ifndef RTOp_USE_MPI

#include "RTOp_mpi.h"

int MPI_Init(int *argc, char ***argv)
{
  return 0;
}

int MPI_Finalize(void)
{
  return 0;
}

int MPI_Comm_size(MPI_Comm comm, int *size)
{
  *size = 1;
  return 0;
}

int MPI_Comm_rank(MPI_Comm comm, int *rank)
{
  *rank = 0;
  return 0;
}

int MPI_Type_struct(int count , int *array_of_blocklengths, MPI_Aint *array_of_displacements
  , MPI_Datatype *array_of_types, MPI_Datatype *data_type)
{
  /* Make the mpi datatype just the extent (needed latter!) */
  int len = 0, extent = 0, k = 0;
  for( k = 0; k < count; ++k ) {
    switch( array_of_types[k] ) {
      case MPI_CHAR:
        len = sizeof(char);
        break;
      case MPI_INT:
        len = sizeof(int);
        break;
      case MPI_DOUBLE:
        len = sizeof(double);
        break;
      default:
        assert(0);
    }
    len = array_of_displacements[k] + array_of_blocklengths[k] * len;
    if( len > extent )
      extent = len;
  }
  *data_type = extent;
  return 0;
}

int MPI_Type_commit(MPI_Datatype *datatype)
{
  return 0;
}

int MPI_Type_free(MPI_Datatype *op)
{
  *op = MPI_DATATYPE_NULL;
  return 0;
}

int MPI_Op_create(MPI_User_function *func, int communitive, MPI_Op *op)
{
  *op = (MPI_Op)*func;
  return 0;
}

int MPI_Op_free( MPI_Op *op)
{
  *op = MPI_OP_NULL;
  return 0;
}

int MPI_Allreduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype
  , MPI_Op op, MPI_Comm comm)
{
  char
    *_sendbuf = sendbuf,
    *_recvbuf = recvbuf;
  int k;
  for( k = 0; k < count * datatype; ++k )
    _recvbuf[k] =_sendbuf[k]; /* just copy bit for bit */
  return 0;
}

int MPI_Reduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op
  , int root, MPI_Comm comm)
{
  char
    *_sendbuf = sendbuf,
    *_recvbuf = recvbuf;
  int k;
  for( k = 0; k < count * datatype; ++k )
    _recvbuf[k] =_sendbuf[k]; /* just copy bit for bit */
  return 0;
}

int MPI_Barrier(MPI_Comm comm)
{
  return 0;
}

int MPI_Bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm )
{
  return 0;
}

int MPI_Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype
         , void* recvbuf, int recvcount, MPI_Datatype recvtype, int root , MPI_Comm comm )
{
  char
    *_sendbuf = sendbuf,
    *_recvbuf = recvbuf;
  int k;
  assert(sendtype == recvtype);
  assert(sendcount == recvcount);
  for( k = 0; k < sendcount * sendtype; ++k )
    _recvbuf[k] =_sendbuf[k]; /* just copy bit for bit */
  return 0;
}

#endif /* RTOp_USE_MPI */
