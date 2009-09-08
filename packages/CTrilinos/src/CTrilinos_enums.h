

#ifndef CTRILINOS_ENUMS_H
#define CTRILINOS_ENUMS_H


#include "CTrilinos_config.h"


#ifdef __cplusplus
extern "C" {
#endif


/* C does not support the C++ bool type, so the enum below will act as a custom
 * boolean type for users using a C compiler.  NOTE, however, that CTrilinos
 * currently expects users to #include directly the Epetra_DataAccess.h and
 * Epetra_CombineMode.h when needed, which may not actually compile with C
 * compilers due to the Doxygen comments! */

#ifdef __cplusplus
typedef bool boolean;

#ifndef FALSE
#define FALSE false
#endif

#ifndef TRUE
#define TRUE true
#endif

#else 
#ifndef TRUE
typedef enum { FALSE=0, TRUE=1 } boolean;
#else
typedef int boolean;
#endif
#endif /* __cplusplus */

/* The enum below lists all the classes that CTrilinos supports.  Classes
 * that are derived from a base class listed below can also be used with some
 * restrictions: (1) the only methods that can be called on such classes are
 * the ones defined by the base class and (2) objects of the derived type can
 * only be passed to methods expecting an argument of the base-class type. */

typedef enum {
    CT_Invalid_ID, 
    CT_Epetra_Distributor_ID, 
    CT_Epetra_SerialComm_ID, 
    CT_Epetra_BLAS_ID, 
    CT_Epetra_Comm_ID, 
    CT_Epetra_Operator_ID, 
    CT_Epetra_MultiVector_ID, 
    CT_Epetra_OffsetIndex_ID, 
    CT_Epetra_Object_ID, 
    CT_Epetra_RowMatrix_ID, 
    CT_Epetra_CompObject_ID, 
    CT_Epetra_Directory_ID, 
    CT_Epetra_Flops_ID, 
    CT_Epetra_SrcDistObject_ID, 
    CT_Epetra_MpiComm_ID, 
    CT_Epetra_CrsMatrix_ID, 
    CT_Epetra_CrsGraph_ID, 
    CT_Epetra_DistObject_ID, 
    CT_Epetra_Vector_ID, 
    CT_Epetra_Export_ID, 
    CT_Epetra_Map_ID, 
    CT_Epetra_BlockMap_ID, 
    CT_Epetra_Import_ID, 
    CT_Epetra_Time_ID, 
    CT_Epetra_JadMatrix_ID, 
    CT_Epetra_LinearProblem_ID, 
    CT_Epetra_LAPACK_ID, 
    CT_Teuchos_CommandLineProcessor_ID, 
    CT_Teuchos_ParameterList_ID, 
    CT_Teuchos_ParameterEntry_ID, 
    CT_Teuchos_any_ID, 
    CT_Amesos_BaseSolver_ID, 
    CT_Amesos_ID
} CTrilinos_Type_ID_t;

/* The type in the struct below is actually used to identify the table in
 * which the object is stored and not the type of the underlying object.  The
 * end user should NEVER change the content of this struct as it will result
 * in segfaults and/or memory leaks.  When the user asks that the object
 * be removed from the table, the type will be changed to CT_Invalid_ID by the
 * destructor function. */

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CTrilinos_Object_ID_t;

/* All structs for specific types are identical but are declared separately
 * rather than as a typedef of the generic one so that they are not inter-
 * changeable.  This forces the compiler to perform type-checking on calls
 * to the wrapper functions and will make debugging easier for the user. */

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_Distributor_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_SerialComm_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_BLAS_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_Comm_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_Operator_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_MultiVector_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_OffsetIndex_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_Object_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_RowMatrix_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_CompObject_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_Directory_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_Flops_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_SrcDistObject_ID_t;

#ifdef HAVE_MPI
typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_MpiComm_ID_t;
#endif /* HAVE_MPI */

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_CrsMatrix_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_CrsGraph_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_DistObject_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_Vector_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_Export_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_Map_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_BlockMap_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_Import_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_Time_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_JadMatrix_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_LinearProblem_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Epetra_LAPACK_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Teuchos_CommandLineProcessor_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Teuchos_ParameterList_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Teuchos_ParameterEntry_ID_t;

typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Teuchos_any_ID_t;

#ifdef HAVE_CTRILINOS_AMESOS
typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Amesos_BaseSolver_ID_t;
#endif /* HAVE_CTRILINOS_AMESOS */

#ifdef HAVE_CTRILINOS_AMESOS
typedef struct {
    CTrilinos_Type_ID_t type;	/* Data type of the object */
    int index;			/* Array index of the object */
    boolean is_const;		/* Whether or not object was declared const */
} CT_Amesos_ID_t;
#endif /* HAVE_CTRILINOS_AMESOS */



typedef enum {
    CT_Epetra_DataAccess_E_Copy,
    CT_Epetra_DataAccess_E_View
} CT_Epetra_DataAccess_E_t;

typedef enum {
    CT_Epetra_CombineMode_E_Add,
    CT_Epetra_CombineMode_E_Zero,
    CT_Epetra_CombineMode_E_Insert,
    CT_Epetra_CombineMode_E_InsertAdd,
    CT_Epetra_CombineMode_E_Average,
    CT_Epetra_CombineMode_E_AbsMax
} CT_Epetra_CombineMode_E_t;

typedef enum {
    CT_ProblemDifficultyLevel_E_easy,
    CT_ProblemDifficultyLevel_E_moderate,
    CT_ProblemDifficultyLevel_E_hard,
    CT_ProblemDifficultyLevel_E_unsure
} CT_ProblemDifficultyLevel_E_t;

typedef enum {
    CT_EParseCommandLineReturn_E_PARSE_SUCCESSFUL = 0,
    CT_EParseCommandLineReturn_E_PARSE_HELP_PRINTED = 1,
    CT_EParseCommandLineReturn_E_PARSE_UNRECOGNIZED_OPTION = 2
} CT_EParseCommandLineReturn_E_t;

typedef enum {
    CT_EOptType_E_OPT_NONE,
    CT_EOptType_E_OPT_BOOL_TRUE,
    CT_EOptType_E_OPT_BOOL_FALSE,
    CT_EOptType_E_OPT_INT,
    CT_EOptType_E_OPT_DOUBLE,
    CT_EOptType_E_OPT_STRING,
    CT_EOptType_E_OPT_ENUM_INT
} CT_EOptType_E_t;

typedef enum {
    CT_EValidateUsed_E_VALIDATE_USED_ENABLED,
    CT_EValidateUsed_E_VALIDATE_USED_DISABLED
} CT_EValidateUsed_E_t;

typedef enum {
    CT_EValidateDefaults_E_VALIDATE_DEFAULTS_ENABLED,
    CT_EValidateDefaults_E_VALIDATE_DEFAULTS_DISABLED
} CT_EValidateDefaults_E_t;


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif
