

#ifndef CTRILINOS_ENUMS_H
#define CTRILINOS_ENUMS_H


#ifdef __cplusplus
extern "C" {
#endif


typedef enum { FALSE=0, TRUE=1 } boolean;

typedef enum {
    CT_Invalid_ID, 
    CT_Epetra_Comm_ID, 
    CT_Epetra_Map_ID, 
    CT_Epetra_MultiVector_ID, 
    CT_Epetra_MpiComm_ID, 
    CT_Epetra_SerialComm_ID, 
    CT_Epetra_BlockMap_ID, 
    CT_Epetra_Vector_ID
} CTrilinos_Type_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
} CTrilinos_Object_ID_t;

/* To keep the function prototypes readable... */
typedef CTrilinos_Object_ID_t CT_Epetra_Comm_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_Map_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_MultiVector_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_MpiComm_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_SerialComm_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_BlockMap_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_Vector_ID_t ;


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif
