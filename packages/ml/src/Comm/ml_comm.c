/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Functions to communicate between processors                          */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : Sept, 1998                                           */
/* ******************************************************************** */

#include "ml_comm.h"

/* ******************************************************************** */
/* Create a communicator for ML                                         */
/* -------------------------------------------------------------------- */
ML_Comm *global_comm;

int ML_Comm_Create( ML_Comm ** com ) 
{
   ML_Comm *com_ptr;

   ML_memory_alloc( (void **) com, sizeof (ML_Comm), "CM1" );
   com_ptr              = (*com);
   com_ptr->ML_id       = ML_ID_COMM;
   com_ptr->ML_mypid    = 0;
   com_ptr->ML_nprocs   = 1;
   com_ptr->USR_comm    = 0;
   com_ptr->USR_sendbytes  = ML_Comm_Send;
   com_ptr->USR_irecvbytes = ML_Comm_Irecv;
   com_ptr->USR_waitbytes  = ML_Comm_Wait;
#ifdef ML_MPI
   MPI_Comm_size(MPI_COMM_WORLD, &(com_ptr->ML_nprocs));
   MPI_Comm_rank(MPI_COMM_WORLD, &(com_ptr->ML_mypid));
   com_ptr->USR_sendbytes  = ML_Comm_Send;
   com_ptr->USR_irecvbytes = ML_Comm_Irecv;
   com_ptr->USR_waitbytes  = ML_Comm_Wait;
#endif


   return 0;
}

/* ******************************************************************** */
/* destroy a data structure for grid communication functions            */
/* -------------------------------------------------------------------- */

int ML_Comm_Destroy( ML_Comm ** com )
{
   if ( (*com)->ML_id != ML_ID_COMM )
   {
      printf("Wrong Comm object to destroy. \n");
      return -1;
   }
   (*com)->ML_id = -1;
   ML_memory_free( (void **) com );
   return 0;
}

/* ******************************************************************** */
/* Check that all communicator functions and variables have been set    */
/* (return -1 if not)                                                   */
/* -------------------------------------------------------------------- */

int ML_Comm_Check( ML_Comm *com_ptr )
{
   int ready_flag = 1;

   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("Wrong Comm object to check. \n");
      return -1;
   }
   if ( com_ptr->USR_irecvbytes == NULL ) ready_flag = 0;
   if ( com_ptr->USR_sendbytes  == NULL ) ready_flag = 0;
   if ( com_ptr->USR_waitbytes  == NULL ) ready_flag = 0;
   if ( com_ptr->ML_mypid  < 0 )          ready_flag = 0;
   if ( com_ptr->ML_nprocs < 0 )          ready_flag = 0;
   if ( com_ptr->USR_comm == 0 )          ready_flag = 0;
   if ( ready_flag == 1 ) return 0;
   else                   return -1;
}

/* ******************************************************************** */
/* Functions to set communicator parameters                             */
/* -------------------------------------------------------------------- */

int ML_Comm_Set_UsrComm( ML_Comm *com_ptr, USR_COMM com )
{
   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("Wrong Comm object to set processor ID (%d).\n",com_ptr->ML_id);
      exit(1);
   }
   com_ptr->USR_comm = com;
   return 0;
}

int ML_Comm_Set_Mypid( ML_Comm *com_ptr, int mypid )
{
   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("Wrong Comm object to set processor ID. \n");
      exit(1);
   }
   com_ptr->ML_mypid = mypid;
   return 0;
}

int ML_Comm_Set_Nprocs( ML_Comm *com_ptr, int nprocs )
{
   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("Wrong Comm object to set number of processors. \n");
      exit(1);
   }
   com_ptr->ML_nprocs = nprocs;
   return 0;
}

int ML_Comm_Set_SendFcn( ML_Comm *com_ptr, int (*func)())
{
   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("Wrong Comm object to set send function. \n");
      exit(1);
   }
   com_ptr->USR_sendbytes = func;
   return 0;
}

int ML_Comm_Set_RecvFcn( ML_Comm *com_ptr, int (*func)())
{
   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("Wrong Comm object to set receive function. \n");
      exit(1);
   }
   com_ptr->USR_irecvbytes = func;
   return 0;
}

int ML_Comm_Set_WaitFcn( ML_Comm *com_ptr, int (*func)())
{
   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("Wrong Comm object to set receive function. \n");
      exit(1);
   }
   com_ptr->USR_waitbytes = func;
   return 0;
}

/************************************************************************/
/* find the maximum integer among the ones in each processor            */
/* -------------------------------------------------------------------- */

int ML_Comm_GmaxInt(ML_Comm *com_ptr, int idata)

{
   int     mask, partner, hbit, msgtype, msgbase=245, nprocs, mypid;
   int     i, k, indata, outdata;
   USR_REQ Request;

   /* ----------------------------------------------------------------- */
   /* check validity of the communication                               */
   /* ----------------------------------------------------------------- */
   
   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_GmaxInt : wrong Comm object. \n");
      exit(1);
   }

   /* ----------------------------------------------------------------- */
   /* get processor information                                         */
   /* ----------------------------------------------------------------- */

   mypid  = com_ptr->ML_mypid;
   nprocs = com_ptr->ML_nprocs;

   /* ----------------------------------------------------------------- */
   /* Find next higher power of 2.                                      */
   /* ----------------------------------------------------------------- */

   for (hbit = 0; (nprocs >> hbit) != 0; hbit++);
   if (nprocs > (1 << hbit)) hbit++;

   /* ----------------------------------------------------------------- */
   /* do a binary collapae (in processor number ascending order)        */
   /* ----------------------------------------------------------------- */

   mask = 0;
   outdata = idata;
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            k = sizeof(int);
            com_ptr->USR_irecvbytes((void*) &indata, k, &partner, 
                                    &msgtype, com_ptr->USR_comm, 
                                    (void *) &Request );
            com_ptr->USR_waitbytes((void*) &indata, k, &partner, 
                                   &msgtype, com_ptr->USR_comm, 
                                   (void *) &Request );
            if ( indata > outdata ) outdata = indata;
         } 
         else if (partner < nprocs)
         {
            k = sizeof(int);
            com_ptr->USR_sendbytes((void*) &outdata, k, partner,
                                   msgtype, com_ptr->USR_comm );
         }
      }
      mask    = mask | (1 << i);
   }

   /* ----------------------------------------------------------------- */
   /* Finally, broadcast this information to every processor in a tree  */
   /* manner.                                                           */
   /* ----------------------------------------------------------------- */

   msgbase = 538;
   mask    = 32767;
   k       = sizeof(int);
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      mask    = mask << 1;
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            com_ptr->USR_sendbytes((void*) &outdata, k, partner, msgtype,
                                   com_ptr->USR_comm );
         } 
         else if (partner < nprocs)
         {
            com_ptr->USR_irecvbytes((void*) &outdata, k, &partner, &msgtype,
                                     com_ptr->USR_comm, (void *) &Request );
            com_ptr->USR_waitbytes((void*) &outdata, k, &partner, &msgtype,
                                   com_ptr->USR_comm, (void *) &Request );
         }
      }
   }
   return outdata;
}

/************************************************************************/
/* find the sum of integers residing in each processor                  */
/* -------------------------------------------------------------------- */

int ML_Comm_GsumInt(ML_Comm *com_ptr, int idata)
{
   int     mask, partner, hbit, msgtype, msgbase=246;
   int     i, k, indata, outdata, nprocs, mypid;
   USR_REQ Request;

   /* ----------------------------------------------------------------- */
   /* check validity of the communication                               */
   /* ----------------------------------------------------------------- */
   
   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_GmaxInt : wrong Comm object. \n");
      exit(1);
   }

   /* ----------------------------------------------------------------- */
   /* get processor information                                         */
   /* ----------------------------------------------------------------- */

   mypid  = com_ptr->ML_mypid;
   nprocs = com_ptr->ML_nprocs;

   /* ----------------------------------------------------------------- */
   /* Find next higher power of 2.                                      */
   /* ----------------------------------------------------------------- */

   for (hbit = 0; (nprocs >> hbit) != 0; hbit++);
   if (nprocs > (1 << hbit)) hbit++;

   /* ----------------------------------------------------------------- */
   /* do a binary collapae (in processor number ascending order)        */
   /* ----------------------------------------------------------------- */

   mask = 0;
   outdata = idata;
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            k = sizeof(int);
            com_ptr->USR_irecvbytes((void*) &indata, k, &partner, &msgtype, 
                                    com_ptr->USR_comm, (void *) &Request );
            com_ptr->USR_waitbytes((void*) &indata, k, &partner, &msgtype, 
                                   com_ptr->USR_comm, (void *) &Request );
            outdata = outdata + indata;
         } 
         else if (partner < nprocs)
         {
            k = sizeof(int);
            com_ptr->USR_sendbytes((void*) &outdata, k, partner, msgtype,
                                   com_ptr->USR_comm );
         }
      }
      mask    = mask | (1 << i);
   }

   /* ----------------------------------------------------------------- */
   /* Finally, broadcast this information to every processor in a tree  */
   /* manner.                                                           */
   /* ----------------------------------------------------------------- */

   msgbase = 539;
   mask    = 32767;
   k       = sizeof(int);
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      mask    = mask << 1;
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            com_ptr->USR_sendbytes((void*) &outdata, k, partner, msgtype,
                                   com_ptr->USR_comm );
         } 
         else if (partner < nprocs)
         {
            com_ptr->USR_irecvbytes((void*) &outdata, k, &partner, &msgtype,
                                    com_ptr->USR_comm, (void *) &Request );
            com_ptr->USR_waitbytes((void*) &outdata, k, &partner, &msgtype,
                                   com_ptr->USR_comm, (void *) &Request );
         }
      }
   }
   return outdata;
}

/************************************************************************/
/* find the sum of doubles  residing in each processor                  */
/* -------------------------------------------------------------------- */

double ML_Comm_GsumDouble(ML_Comm *com_ptr, double ddata)
{
   int     mask, partner, hbit, msgtype, msgbase=247;
   int     i, k, nprocs, mypid;
   double  indata, outdata;
   USR_REQ Request;

   /* ----------------------------------------------------------------- */
   /* check validity of the communication                               */
   /* ----------------------------------------------------------------- */
   
   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_GmaxDouble : wrong Comm object. \n");
      exit(1);
   }

   /* ----------------------------------------------------------------- */
   /* get processor information                                         */
   /* ----------------------------------------------------------------- */

   mypid  = com_ptr->ML_mypid;
   nprocs = com_ptr->ML_nprocs;

   /* ----------------------------------------------------------------- */
   /* Find next higher power of 2.                                      */
   /* ----------------------------------------------------------------- */

   for (hbit = 0; (nprocs >> hbit) != 0; hbit++);
   if (nprocs > (1 << hbit)) hbit++;

   /* ----------------------------------------------------------------- */
   /* do a binary collapae (in processor number ascending order)        */
   /* ----------------------------------------------------------------- */

   mask = 0;
   outdata = ddata;
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            k = sizeof(double);
            com_ptr->USR_irecvbytes((void*) &indata, k, &partner, &msgtype, 
                                    com_ptr->USR_comm, (void *) &Request );
            com_ptr->USR_waitbytes((void*) &indata, k, &partner, &msgtype, 
                                   com_ptr->USR_comm, (void *) &Request );
            outdata = outdata + indata;
         } 
         else if (partner < nprocs)
         {
            k = sizeof(double);
            com_ptr->USR_sendbytes((void*) &outdata, k, partner, msgtype,
                                   com_ptr->USR_comm );
         }
      }
      mask    = mask | (1 << i);
   }

   /* ----------------------------------------------------------------- */
   /* Finally, broadcast this information to every processor in a tree  */
   /* manner.                                                           */
   /* ----------------------------------------------------------------- */

   msgbase = 539;
   mask    = 32767;
   k       = sizeof(double);
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      mask    = mask << 1;
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            com_ptr->USR_sendbytes((void*) &outdata, k, partner, msgtype,
                                   com_ptr->USR_comm );
         } 
         else if (partner < nprocs)
         {
            com_ptr->USR_irecvbytes((void*) &outdata, k, &partner, &msgtype,
                                    com_ptr->USR_comm, (void *) &Request );
            com_ptr->USR_waitbytes((void*) &outdata, k, &partner, &msgtype,
                                   com_ptr->USR_comm, (void *) &Request );
         }
      }
   }
   return outdata;
}

/************************************************************************/
/* sum an integer vector across all processors                          */
/*----------------------------------------------------------------------*/

int ML_Comm_GsumVecInt(ML_Comm *com_ptr, int *idata, int *itmp, int leng)
{
   int     mask, partner, hbit, msgtype, msgbase=247;
   int     i, j, k, nprocs, mypid;
   USR_REQ Request;

   /* ----------------------------------------------------------------- */
   /* check validity of the communication                               */
   /* ----------------------------------------------------------------- */
   
   if ( com_ptr->ML_id != ML_ID_COMM ) {
      printf("ML_Comm_GmaxInt : wrong Comm object. \n");
      exit(1);
   }

   /* ----------------------------------------------------------------- */
   /* get processor information                                         */
   /* ----------------------------------------------------------------- */

   mypid  = com_ptr->ML_mypid;
   nprocs = com_ptr->ML_nprocs;

   /* ----------------------------------------------------------------- */
   /* Find next higher power of 2.                                      */
   /* ----------------------------------------------------------------- */

   for (hbit = 0; (nprocs >> hbit) != 0; hbit++);
   if (nprocs > (1 << hbit)) hbit++;

   /* ----------------------------------------------------------------- */
   /* do a binary collapae (in processor number ascending order)        */
   /* ----------------------------------------------------------------- */

   mask = 0;

   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            k = leng * sizeof(int);
            com_ptr->USR_irecvbytes((void*) itmp, k, &partner, &msgtype,
                                    com_ptr->USR_comm, (void *) &Request);
            com_ptr->USR_waitbytes((void*) itmp, k, &partner, &msgtype,
                                   com_ptr->USR_comm, (void *) &Request);
            for ( j = 0; j < leng; j++ ) idata[j] += itmp[j];
         } 
         else if (partner < nprocs)
         {
            k = leng * sizeof(int);
            com_ptr->USR_sendbytes((void*) idata, k, partner, msgtype,
                                   com_ptr->USR_comm );
         }
      }
      mask    = mask | (1 << i);
   }

   /* ----------------------------------------------------------------- */
   /* Finally, broadcast this information to every processor in a tree  */
   /* manner.                                                           */
   /* ----------------------------------------------------------------- */

   msgbase = 540;
   mask    = 32767;
   k       = leng * sizeof(int);
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      mask    = mask << 1;
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            com_ptr->USR_sendbytes((void*) idata, k, partner, msgtype,
                                   com_ptr->USR_comm );
         } 
         else if (partner < nprocs)
         {
            com_ptr->USR_irecvbytes((void*) idata, k, &partner, &msgtype,
                                    com_ptr->USR_comm, (void *) &Request);
            com_ptr->USR_waitbytes((void*) idata, k, &partner, &msgtype,
                                   com_ptr->USR_comm, (void *) &Request);
         }
      }
   }
   return 0;
}

/************************************************************************/
/* This is a modification from Tuminaro's AZ_gappend_int subroutine.    */
/* The modification is done so that the data are arranged according to  */
/* processor number.                                                    */
/*----------------------------------------------------------------------*/

int ML_Comm_GappendInt(ML_Comm *com_ptr, int *vals, int *cur_length, 
                    int total_length)
{
   int     mask, partner, hbit, msgtype, msgbase=145;
   int     i, k, nbytes, mypid, nprocs;
   USR_REQ Request;

   /* ----------------------------------------------------------------- */
   /* check validity of the communication                               */
   /* ----------------------------------------------------------------- */
   
   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_GmaxInt : wrong Comm object. \n");
      exit(1);
   }

   /* ----------------------------------------------------------------- */
   /* get processor information                                         */
   /* ----------------------------------------------------------------- */

   mypid  = com_ptr->ML_mypid;
   nprocs = com_ptr->ML_nprocs;

   /* ----------------------------------------------------------------- */
   /* Find next higher power of 2.                                      */
   /* ----------------------------------------------------------------- */

   for (hbit = 0; (nprocs >> hbit) != 0; hbit++);
   if (nprocs > (1 << hbit)) hbit++;

   /* ----------------------------------------------------------------- */
   /* do a binary collapae (in processor number ascending order)        */
   /* ----------------------------------------------------------------- */

   mask = 0;
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            k = (total_length - (*cur_length)) * sizeof(int);
            com_ptr->USR_irecvbytes((void*)&(vals[*cur_length]), k, 
                                   &partner, &msgtype, 
                                   com_ptr->USR_comm, (void *) &Request);
            nbytes = com_ptr->USR_waitbytes((void*)&(vals[*cur_length]), k, 
                                   &partner, &msgtype, 
                                   com_ptr->USR_comm, (void *) &Request);
            (*cur_length) += (nbytes / sizeof(int));
         } 
         else if (partner < nprocs)
         {
            k = (*cur_length) * sizeof(int);
            com_ptr->USR_sendbytes((void*) vals, k, partner, msgtype, 
                                   com_ptr->USR_comm );
         }
      }
      mask    = mask | (1 << i);
   }

   /* ----------------------------------------------------------------- */
   /* Finally, broadcast this information to every processor in a tree  */
   /* manner.                                                           */
   /* ----------------------------------------------------------------- */

   msgbase = 438;
   mask    = 32767;
   k       = total_length * sizeof(int);
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      mask    = mask << 1;
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            com_ptr->USR_sendbytes((void*) vals, k, partner, msgtype, 
                                   com_ptr->USR_comm );
         }
         else if (partner < nprocs)
         {
            com_ptr->USR_irecvbytes((void*) vals, k, &partner, &msgtype, 
                                   com_ptr->USR_comm, (void *) &Request);
            com_ptr->USR_waitbytes((void*) vals, k, &partner, &msgtype, 
                                   com_ptr->USR_comm, (void *) &Request);
         }
      }
   }
   (*cur_length) = total_length;
   return 0;
}

/************************************************************************/
/* This is a modification from Tuminaro's AZ_gappend_double subroutine. */
/* The modification is done so that the data are arranged according to  */
/* processor number.                                                    */
/*----------------------------------------------------------------------*/

int ML_Comm_GappendDouble(ML_Comm *com_ptr, double *vals, int *cur_length, 
                            int total_length)
{
   int     mask, partner, hbit, msgtype, msgbase=245;
   int     i, k, nbytes, nprocs, mypid;
   USR_REQ Request;

   /* ----------------------------------------------------------------- */
   /* check validity of the communication                               */
   /* ----------------------------------------------------------------- */
   
   if ( com_ptr->ML_id != ML_ID_COMM )
   {
      printf("ML_Comm_GmaxInt : wrong Comm object. \n");
      exit(1);
   }

   /* ----------------------------------------------------------------- */
   /* get processor information                                         */
   /* ----------------------------------------------------------------- */

   mypid  = com_ptr->ML_mypid;
   nprocs = com_ptr->ML_nprocs;

   /* ----------------------------------------------------------------- */
   /* Find next higher power of 2.                                      */
   /* ----------------------------------------------------------------- */

   for (hbit = 0; (nprocs >> hbit) != 0; hbit++);
   if (nprocs > (1 << hbit)) hbit++;

   /* ----------------------------------------------------------------- */
   /* do a binary collapae (in processor number ascending order)        */
   /* ----------------------------------------------------------------- */

   mask = 0;
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            k = (total_length - (*cur_length)) * sizeof(double);
            com_ptr->USR_irecvbytes((void*)&(vals[*cur_length]), k, 
                                    &partner, &msgtype, 
                                    com_ptr->USR_comm, (void *) &Request);
            nbytes = com_ptr->USR_waitbytes((void*)&(vals[*cur_length]), k, 
                                    &partner, &msgtype, 
                                    com_ptr->USR_comm, (void *) &Request);
            (*cur_length) += (nbytes / sizeof(double));
         }
         else if (partner < nprocs)
         {
            k = (*cur_length) * sizeof(double);
            com_ptr->USR_sendbytes((void*) vals, k, partner, msgtype, 
                                   com_ptr->USR_comm );
         }
      }
      mask    = mask | (1 << i);
   }

   /* ----------------------------------------------------------------- */
   /* Finally, broadcast this information to every processor in a tree  */
   /* manner.                                                           */
   /* ----------------------------------------------------------------- */

   msgbase = 538;
   mask    = 32767;
   k       = total_length * sizeof(double);
   for ( i = 0; i < hbit; i++ )
   {
      msgtype = msgbase + i;
      partner = mypid ^ (1 << i);
      mask    = mask << 1;
      if ((mypid & mask) == 0)
      {
         if (((mypid & (1 << i)) == 0) && (partner < nprocs))
         {
            com_ptr->USR_sendbytes((void*) vals, k, partner, msgtype, 
                                   com_ptr->USR_comm );
         }
         else if (partner < nprocs)
         {
            com_ptr->USR_irecvbytes((void*) vals, k, &partner, &msgtype, 
                                    com_ptr->USR_comm, (void *) &Request );
            com_ptr->USR_waitbytes((void*) vals, k, &partner, &msgtype, 
                                    com_ptr->USR_comm, (void *) &Request );
         }
      }
   }
   (*cur_length) = total_length;
   return 0;
}

/**************************************************************************/
/* communication subroutines                                              */
/*------------------------------------------------------------------------*/
int ML_Comm_Irecv(void* buf, unsigned int count, int *src,
                  int *mid, USR_COMM comm, USR_REQ *request )
{
   int err = 0;
#ifdef ML_MPI
   if (*mid == -1) *mid = MPI_ANY_TAG;
   if (*src == -1) *src = MPI_ANY_SOURCE;
   err = MPI_Irecv(buf,count,MPI_BYTE,*src,*mid,MPI_COMM_WORLD,request);
#endif
   return err;
}

int ML_Comm_Wait (void* buf, unsigned int count, int *src,
                  int *mid, USR_COMM comm, USR_REQ *request )
{
   int        return_cnt = 0;
#ifdef ML_MPI
   MPI_Status status;
   MPI_Wait(request, &status);
   MPI_Get_count(&status, MPI_BYTE, &return_cnt);
   *src = status.MPI_SOURCE;
   *mid = status.MPI_TAG;
#endif
   return return_cnt;
}

int ML_Comm_Send(void* buf, unsigned int count, int dest, int mid,
                 USR_COMM comm )
{
   int err = 0;
#ifdef ML_MPI
   err = MPI_Send(buf, count, MPI_BYTE, dest, mid, MPI_COMM_WORLD);
#endif
   return err;
}

