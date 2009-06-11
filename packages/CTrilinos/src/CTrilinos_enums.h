

#ifndef CTRILINOS_ENUMS_H
#define CTRILINOS_ENUMS_H


#ifdef __cplusplus
extern "C" {
#endif


#ifdef __cplusplus
typedef bool boolean;
#define FALSE false
#define TRUE true
#else
typedef enum { FALSE=0, TRUE=1 } boolean;
#endif

typedef enum {
    CT_Invalid_ID, 
    CT_Epetra_Distributor_ID, 
    CT_Epetra_SerialCommData_ID, 
    CT_Epetra_SerialComm_ID, 
    CT_Epetra_BLAS_ID, 
    CT_Epetra_Comm_ID, 
    CT_Epetra_Operator_ID, 
    CT_Epetra_MultiVector_ID, 
    CT_Epetra_OffsetIndex_ID, 
    CT_Epetra_Object_ID, 
    CT_Epetra_Data_ID, 
    CT_Epetra_RowMatrix_ID, 
    CT_Epetra_CompObject_ID, 
    CT_Epetra_Directory_ID, 
    CT_Epetra_BlockMapData_ID, 
    CT_Epetra_Flops_ID, 
    CT_Epetra_SrcDistObject_ID, 
    CT_Epetra_MpiComm_ID, 
    CT_Epetra_MpiCommData_ID, 
    CT_Epetra_CrsMatrix_ID, 
    CT_Epetra_CrsGraph_ID, 
    CT_Epetra_DistObject_ID, 
    CT_Epetra_Export_ID, 
    CT_Epetra_Vector_ID, 
    CT_Epetra_Map_ID, 
    CT_Epetra_CrsGraphData_ID, 
    CT_Epetra_BlockMap_ID, 
    CT_Epetra_Import_ID
} CTrilinos_Type_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
} CTrilinos_Object_ID_t;

/* To keep the function prototypes readable... */
typedef CTrilinos_Object_ID_t CT_Epetra_Distributor_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_SerialCommData_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_SerialComm_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_BLAS_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_Comm_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_Operator_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_MultiVector_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_OffsetIndex_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_Object_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_Data_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_RowMatrix_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_CompObject_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_Directory_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_BlockMapData_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_Flops_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_SrcDistObject_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_MpiComm_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_MpiCommData_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_CrsMatrix_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_CrsGraph_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_DistObject_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_Export_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_Vector_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_Map_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_CrsGraphData_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_BlockMap_ID_t ;
typedef CTrilinos_Object_ID_t CT_Epetra_Import_ID_t ;


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif
