/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * This software is distributed under the GNU Lesser General Public License. *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <values.h>
#include "iohb.h"


main (int argc, char *argv[])
    {
    int nRow, nCol, nz, nEdge, nPin;
    int *rowindex = NULL, *colstart = NULL;
    int maxEdgeLen = 0, minEdgeLen = MAXINT;
    int single = 0;  /* Num of cols with only one non-zero. */
    double *val = NULL;
    int i;
    int column;
    char filename[100];
    FILE *hg;

    readHB_newmat_double(argv[1], &nRow, &nCol, &nz, &colstart, &rowindex, &val);

    /* write hypergraph */
    sprintf (filename, "./%s.hg", argv[1]);
    hg = fopen (filename, "w");
    fprintf (hg, "%60s \n", " ");  /* temporary header line */

    nEdge = 0;
    nPin = 0;
    for (column = 0; column < nCol; column++)
       {
       if ((colstart[column+1] - colstart[column]) > maxEdgeLen) 
          maxEdgeLen = colstart[column+1] - colstart[column];
       if ((colstart[column+1] - colstart[column]) < minEdgeLen) 
          minEdgeLen = colstart[column+1] - colstart[column];

       if  ((colstart[column+1] - colstart[column]) < 2)
          {
          single++;
#ifdef REMOVE_ONE_NODE_HYPEREDGES
          continue;  /* surpress self edges */
#endif /* REMOVE_ONE_NODE_HYPEREDGES */
          }
       nEdge++;

       for (i = colstart[column]; i < colstart[column+1]; i++)
          {
          fprintf (hg, "%d ", rowindex[i]+1);  /* write vertex */
          nPin++;
          }
       fprintf (hg, "\n");
       }
    rewind (hg);
    fprintf (hg, "%d %d %d 00", nRow, nEdge, nPin); /* header line for hg */
    fclose (hg);
    free(colstart);
    free(rowindex);
    if (val) free(val);
    printf("maxEdgeLen = %d  minEdgeLen = %d  single = %d\n", 
      maxEdgeLen, minEdgeLen, single);
    }
