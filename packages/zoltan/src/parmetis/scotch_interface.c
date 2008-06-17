/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2008 Sandia National Laboratories.                          *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <ctype.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "all_allo_const.h"
#include "params_const.h"
#include "order_const.h"
#include "third_library.h"
#include "scotch_interface.h"

  /**********  parameters structure for Scotch methods **********/
static PARAM_VARS Scotch_params[] = {
  { "SCOTCH_METHOD", NULL, "STRING", 0 },
  { "SCOTCH_STRAT", NULL, "STRING", 0 },
  { "SCOTCH_STRAT_FILE", NULL, "STRING", 0 },
  { NULL, NULL, NULL, 0 } };

static int Zotlan_Scotch_Bind_Param(ZZ * zz, char *alg, char **strat);


/***************************************************************************
 *  The Scotch ordering routine piggy-backs on the Scotch
 *  partitioning routines.
 **************************************************************************/

/*
 * TODO: at this time, only distributed computations are allowed.
 */

int Zoltan_Scotch_Order(
  ZZ *zz,               /* Zoltan structure */
  int num_obj,		/* Number of (local) objects to order. */
  ZOLTAN_ID_PTR gids,   /* List of global ids (local to this proc) */
  /* The application must allocate enough space */
  ZOLTAN_ID_PTR lids,   /* List of local ids (local to this proc) */
/* The application must allocate enough space */
  int *rank,		/* rank[i] is the rank of gids[i] */
  ZOOS *order_opt 	/* Ordering options, parsed by Zoltan_Order */
)
{
  static char *yo = "Zoltan_Scotch_Order";
  int n, ierr;
  ZOLTAN_Output_Order ord;
  ZOLTAN_Third_Graph gr;
  ZOLTAN_Third_Vsize vsp;
  SCOTCH_Strat        stradat;
  SCOTCH_Dgraph       grafdat;
  SCOTCH_Dordering    ordedat;
  int edgelocnbr = 0;


  MPI_Comm comm = zz->Communicator;/* want to risk letting external packages */
  int timer_p = 0;
  int get_times = 0;
  int use_timers = 0;
  double times[5];

  char alg[MAX_PARAM_STRING_LEN+1];
  char *strat = NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);

  memset(&gr, 0, sizeof(ZOLTAN_Third_Graph));
  memset(&vsp, 0, sizeof(ZOLTAN_Third_Vsize));
  memset(&ord, 0, sizeof(ZOLTAN_Output_Order));

  strcpy (alg, "NODEND");
  ierr = Zotlan_Scotch_Bind_Param(zz, alg, &strat);
  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ierr);
  }

  ord.order_info = &(zz->Order);
  ord.order_opt = order_opt;

  if (!order_opt){
    /* If for some reason order_opt is NULL, allocate a new ZOOS here. */
    /* This should really never happen. */
    order_opt = (ZOOS *) ZOLTAN_MALLOC(sizeof(ZOOS));
    strcpy(order_opt->method,"SCOTCH");
  }

  /* Scotch only computes the rank vector */
  order_opt->return_args = RETURN_RANK;

  /* Check that num_obj equals the number of objects on this proc. */
  /* This constraint may be removed in the future. */
  n = zz->Get_Num_Obj(zz->Get_Num_Obj_Data, &ierr);
  if ((ierr!= ZOLTAN_OK) && (ierr!= ZOLTAN_WARN)){
    /* Return error code */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Get_Num_Obj returned error.");
    return(ZOLTAN_FATAL);
  }
  if (n != num_obj){
    /* Currently this is a fatal error. */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Input num_obj does not equal the number of objects.");
    return(ZOLTAN_FATAL);
  }

  /* Do not use weights for ordering */
  gr.obj_wgt_dim = -1;
  gr.edge_wgt_dim = -1;
  gr.num_obj = num_obj;

  /* Check what ordering type is requested */
  if (order_opt){
    if (strcmp(order_opt->order_type, "LOCAL") == 0)
      gr.graph_type = - LOCAL_GRAPH;
    else if (strcmp(order_opt->order_type, "GLOBAL") == 0)
      gr.graph_type =  - GLOBAL_GRAPH;
    else
      gr.graph_type = - NO_GRAPH;
  }
  gr.get_data = 1;

  /* If reorder is true, we already have the id lists. Ignore weights. */
  if ((order_opt && order_opt->reorder))
    gr.id_known = 1;                        /* We already have global_ids and local_ids */

  timer_p = Zoltan_Preprocess_Timer(zz, &use_timers);

    /* Start timer */
  get_times = (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME);
  if (get_times){
    MPI_Barrier(zz->Communicator);
    times[0] = Zoltan_Time(zz->Timer);
  }

  ierr = Zoltan_Preprocess_Graph(zz, &gids, &lids,  &gr, NULL, NULL, &vsp);

  if (SCOTCH_dgraphInit (&grafdat, comm) != 0) {
    Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
    ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL, "Cannot initialize Scotch graph.");
  }

  edgelocnbr =  gr.xadj[gr.num_obj];
  if (SCOTCH_dgraphBuild (&grafdat, 0, gr.num_obj, gr.num_obj, gr.xadj, gr.xadj + 1,
			  NULL, NULL,edgelocnbr, edgelocnbr, gr.adjncy, NULL, NULL) != 0) {
    Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
    ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL, "Cannot construct Scotch graph.");
  }

  SCOTCH_stratInit (&stradat);
  if (strat != NULL) {
    if ((SCOTCH_stratDgraphOrder (&stradat, strat)) != 0) {
      Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
      ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL, "Invalid Scotch strat.");
    }
  }

  if (SCOTCH_dgraphOrderInit (&grafdat, &ordedat) != 0) {
    Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
    ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL, "Cannot construct Scotch graph.");
  }

  /* Get a time here */
  if (get_times) times[1] = Zoltan_Time(zz->Timer);

/*   if (gr.graph_type==GLOBAL_GRAPH){ */
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the PT-Scotch library");
    ierr = SCOTCH_dgraphOrderCompute (&grafdat, &ordedat, &stradat);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the PT-Scotch library");
/*   } */
/*   else { */
/*     ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the Scotch library"); */
/*     METIS_NodeND (&gr.num_obj, gr.xadj, gr.adjncy, */
/* 		  &numflag, options, ord.rank, ord.iperm); */
/*     ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the Scotch library"); */
/*   } */

  if (ierr != 0) {
    Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
    ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL, "Cannot compute Scotch ordering.");
  }

  /* Get a time here */
  if (get_times) times[2] = Zoltan_Time(zz->Timer);

  /* Allocate space for rank array */
  ord.rank = (indextype *) ZOLTAN_MALLOC(gr.num_obj*sizeof(indextype));
  if (!ord.rank){
    /* Not enough memory */
    Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
    ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Out of memory.");
  }

  if (SCOTCH_dgraphOrderPerm (&grafdat, &ordedat, ord.rank) != 0) {
    Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
    ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL, "Cannot compute Scotch rank array");
  }

  ierr = Zoltan_Postprocess_Graph (zz, gids, lids, &gr, NULL, NULL, &vsp, &ord, NULL);

  /* Get a time here */
  if (get_times) times[3] = Zoltan_Time(zz->Timer);

  SCOTCH_dgraphOrderExit (&grafdat, &ordedat);
  SCOTCH_stratExit (&stradat);
  SCOTCH_dgraphExit (&grafdat);

  if (get_times) Zoltan_Third_DisplayTime(zz, times);

  if (use_timers)
    ZOLTAN_TIMER_STOP(zz->ZTime, timer_p, zz->Communicator);

  if (order_opt->return_args&RETURN_RANK)
    memcpy(rank, ord.rank, gr.num_obj*sizeof(indextype));
  ZOLTAN_FREE(&ord.rank);
  ZOLTAN_FREE(strat);
  Zoltan_Third_Exit(&gr, NULL, NULL, &vsp, NULL, NULL);

  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ZOLTAN_OK);
}


/************************
 * Auxiliary function used to parse scotch specific parameters
 ***************************************/

static int Zotlan_Scotch_Bind_Param(ZZ* zz, char *alg, char **strat)
{
  char stratsmall[MAX_PARAM_STRING_LEN+1];
  char stratfilename[MAX_PARAM_STRING_LEN+1];

  *strat = NULL;
  stratsmall[0] = stratfilename[0] = '\0';
  Zoltan_Bind_Param(Scotch_params, "SCOTCH_METHOD",
		    (void *) alg);
  Zoltan_Bind_Param(Scotch_params, "SCOTCH_STRAT",
		    (void *) stratsmall);
  Zoltan_Bind_Param(Scotch_params, "SCOTCH_STRAT_FILE",
		    (void *) stratfilename);
  Zoltan_Assign_Param_Vals(zz->Params, Scotch_params, zz->Debug_Level,
			   zz->Proc, zz->Debug_Proc);

  if ((strlen(stratsmall) > 0) && (strlen(stratfilename) > 0)) {
    ZOLTAN_THIRD_ERROR(ZOLTAN_WARN,
		       "SCOTCH_STRAT and SCOTCH_STRAT_FILE both defined: ignoring\n");
  }

  if (strlen(stratsmall) > 0) {
    *strat = (char *) ZOLTAN_MALLOC((strlen(stratsmall)+1)*sizeof(char));
    strcpy (*strat, stratsmall);
    return (ZOLTAN_OK);
  }

  if (strlen(stratfilename) > 0) {
    long size;
    FILE *stratfile;

    stratfile = fopen(stratfilename, "r");
    if (stratfile == NULL) {
      ZOLTAN_THIRD_ERROR(ZOLTAN_WARN,
			 "Cannot open Scotch strategy file\n");
    }

    fseek(stratfile, (long)0, SEEK_END);
    size = ftell(stratfile);
    *strat = (char *) ZOLTAN_MALLOC((size+1)*sizeof(char));
    if (*strat == NULL) {
      ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }
    fseek(stratfile, (long)0, SEEK_SET);
    fgets (*strat, size, stratfile);
    strat[size] = '\0';
    fclose(stratfile);

    return (ZOLTAN_OK);
  }

  return (ZOLTAN_OK);
}


/*********************************************************************/
/* Scotch parameter routine                                          */
/*********************************************************************/

int Zoltan_Scotch_Set_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
  int status, i;
  PARAM_UTYPE result;         /* value returned from Check_Param */
  int index;                  /* index returned from Check_Param */
  char *valid_methods[] = {
    "NODEND", /* for nested dissection ordering */
    NULL };

  status = Zoltan_Check_Param(name, val, Scotch_params, &result, &index);
  if (status == 0){
    /* OK so far, do sanity check of parameter values */

    if (strcmp(name, "SCOTCH_METHOD") == 0){
      status = 2;
      for (i=0; valid_methods[i] != NULL; i++){
	if (strcmp(val, valid_methods[i]) == 0){
	  status = 0;
	  break;
	}
      }
    }
  }
  return(status);
}

#ifdef __cplusplus
}
#endif

