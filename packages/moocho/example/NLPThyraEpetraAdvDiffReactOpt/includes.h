class Epetra_Comm;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_FECrsMatrix;
class Epetra_FEVector;
class Epetra_SerialDenseMatrix;
class Epetra_SerialDenseVector;
class Epetra_IntSerialDenseVector;
class Epetra_IntSerialDenseMatrix;
namespace Teuchos {
 template <class T> class RefCountPtr;
}
class Usr_Par;

const double GLp_pi = 3.14159265358979323846;

bool CrsMatrix2MATLAB(const Epetra_CrsMatrix &, ostream &);

bool Vector2MATLAB( const Epetra_Vector &, ostream &);

bool FEVector2MATLAB( const Epetra_FEVector &, ostream &);

int quadrature(const int, const int, Epetra_SerialDenseMatrix &,
               Epetra_SerialDenseVector &);

int meshreader(const Epetra_Comm &,
               Epetra_IntSerialDenseVector &,
               Epetra_SerialDenseMatrix &,
               Epetra_IntSerialDenseVector &,
               Epetra_SerialDenseMatrix &,
               Epetra_IntSerialDenseMatrix &,
               Epetra_IntSerialDenseMatrix &,
               const char geomFileBase[],
               const bool trace = true,
               const bool dumpAll = false
               );

int lassembly(const Epetra_SerialDenseMatrix &,
              const Epetra_SerialDenseVector &,
              const Epetra_SerialDenseMatrix &,
              const Epetra_SerialDenseVector &,
              const Epetra_SerialDenseVector &,
              const Usr_Par &,
              Epetra_SerialDenseMatrix &,
              Epetra_SerialDenseVector &);

int assemblyFECrs(const Epetra_Comm &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseMatrix &,
                  const Epetra_IntSerialDenseMatrix &,
                  Teuchos::RefCountPtr<Epetra_FECrsMatrix> &,
                  Teuchos::RefCountPtr<Epetra_FEVector> &);

int assemblyFECrs(const Epetra_Comm &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseMatrix &,
                  const Epetra_IntSerialDenseMatrix &,
                  Teuchos::RefCountPtr<Epetra_FECrsMatrix> &,
                  Teuchos::RefCountPtr<Epetra_FEVector> &,
                  bool);

int assemble(const Epetra_Comm &,
             const Epetra_IntSerialDenseVector &,
             const Epetra_SerialDenseMatrix &,
             const Epetra_IntSerialDenseVector &,
             const Epetra_SerialDenseMatrix &,
             const Epetra_IntSerialDenseMatrix &,
             const Epetra_IntSerialDenseMatrix &,
             Teuchos::RefCountPtr<Epetra_FECrsMatrix> &,
             Teuchos::RefCountPtr<Epetra_FECrsMatrix> &,
             Teuchos::RefCountPtr<Epetra_FEVector> &);

int assemble_bdry(const Epetra_Comm &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseMatrix &,
                  const Epetra_IntSerialDenseMatrix &,
                  Teuchos::RefCountPtr<Epetra_FECrsMatrix> &,
                  Teuchos::RefCountPtr<Epetra_FECrsMatrix> &);

int nonlinvec(const Epetra_Comm &,
              const Epetra_IntSerialDenseVector &,
              const Epetra_SerialDenseMatrix &,
              const Epetra_IntSerialDenseVector &,
              const Epetra_SerialDenseMatrix &,
              const Epetra_IntSerialDenseMatrix &,
              const Teuchos::RefCountPtr<const Epetra_MultiVector> &,
              Teuchos::RefCountPtr<Epetra_FEVector> &);

int nonlinjac(const Epetra_Comm &,
              const Epetra_IntSerialDenseVector &,
              const Epetra_SerialDenseMatrix &,
              const Epetra_IntSerialDenseVector &,
              const Epetra_SerialDenseMatrix &,
              const Epetra_IntSerialDenseMatrix &,
              const Teuchos::RefCountPtr<const Epetra_MultiVector> &,
              Teuchos::RefCountPtr<Epetra_FECrsMatrix> &);

int nonlinhessvec(const Epetra_Comm &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseMatrix &,
                  const Teuchos::RefCountPtr<const Epetra_MultiVector> &,
                  const Teuchos::RefCountPtr<const Epetra_MultiVector> &,
                  const Teuchos::RefCountPtr<const Epetra_MultiVector> &,
                  Teuchos::RefCountPtr<Epetra_FEVector> &);


ostream& operator <<(ostream &, const Usr_Par &);
