/* ///////////////////////////////////////////// */
/* RTOp_mpi.h */
/* */
/* Copyright (C) 2001 Roscoe Ainsworth Bartlett */
/* */
/* This is free software; you can redistribute it and/or modify it */
/* under the terms of the "Artistic License" (see the web site */
/*   http://www.opensource.org/licenses/artistic-license.html). */
/* This license is spelled out in the file COPYING. */
/* */
/* This software is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* above mentioned "Artistic License" for more details. */
/* */
/* MPI declarations used by RTOp example program. */
/* These where taken from mpich for Windows NT. */
/* */

#ifndef MPI_H
#define MPI_H

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////// */
/* MPI declarations */

#define MPI_Aint int
typedef int MPI_Datatype;
#define MPI_CHAR           ((MPI_Datatype)1)
#define MPI_INT            ((MPI_Datatype)6)
#define MPI_DOUBLE         ((MPI_Datatype)11)
typedef int MPI_Comm;
#define MPI_COMM_WORLD 91
#define MPI_COMM_NULL      ((MPI_Comm)0)
typedef int MPI_Op;
#define MPI_OP_NULL        ((MPI_Op)0)
#define MPI_DATATYPE_NULL  ((MPI_Datatype)0)
typedef void (MPI_User_function) ( void *, void *, int *, MPI_Datatype * ); 

/* ////////////////////////// */
/* MPI functions */

#define EXPORT_MPI_API
EXPORT_MPI_API int MPI_Init(int *, char ***);
EXPORT_MPI_API int MPI_Finalize(void);
EXPORT_MPI_API int MPI_Comm_size(MPI_Comm, int *);
EXPORT_MPI_API int MPI_Comm_rank(MPI_Comm, int *);
EXPORT_MPI_API int MPI_Type_struct(int, int *, MPI_Aint *, MPI_Datatype *, MPI_Datatype *);
EXPORT_MPI_API int MPI_Type_commit(MPI_Datatype *);
EXPORT_MPI_API int MPI_Type_free(MPI_Datatype *);
EXPORT_MPI_API int MPI_Reduce(void* , void*, int, MPI_Datatype, MPI_Op, int, MPI_Comm);
EXPORT_MPI_API int MPI_Op_create(MPI_User_function *, int, MPI_Op *);
EXPORT_MPI_API int MPI_Op_free( MPI_Op *);
EXPORT_MPI_API int MPI_Allreduce(void* , void*, int, MPI_Datatype, MPI_Op, MPI_Comm);
EXPORT_MPI_API int MPI_Barrier(MPI_Comm);
EXPORT_MPI_API int MPI_Bcast(void*, int, MPI_Datatype, int, MPI_Comm );
EXPORT_MPI_API int MPI_Gather(void* , int, MPI_Datatype, void*, int, MPI_Datatype, int, MPI_Comm); 

#ifdef __cplusplus
}
#endif

#endif /* MPI_H */
