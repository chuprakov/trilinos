#ifndef TSF_H
#define TSF_H

#include "TSFDefs.h"

#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFLinearSolver.h"
#include "TSFCommandLine.h"
#include "PetraMatrix.h"
#include "LAPACKGeneralMatrix.h"
#include "TSFLinearProblem.h"
#include "TSFMPI.h"
#include "TSFOut.h"
#include "ILUKPreconditionerFactory.h"
#include "TSFPreconditionerFactory.h"
#include "TSFVectorType.h"
#include "PetraVectorType.h"
#include "DenseSerialVectorType.h"
#include "AZTECSolver.h"
#include "TSFDeferredLinearCombination.h"
#include "DLARAN.h"
#include "SystemRand.h"
#include "CGSolver.h"
#include "PetraMatrix.h"
#include "PetraVectorType.h"
#include "DenseSerialVectorType.h"
#include "DirectSolver.h"
#include "TSFPreconditionerFactory.h"
#include "ILUKPreconditionerFactory.h"
#include "TSFMatrixReader.h"
#include "MatrixMarketReader.h"
#include "TSFUtils.h"

using namespace TSF;

namespace TSF
{
  void init(int* argc, void*** argv);

  void finalize();

  void handleError(exception& e, const string& filename);
}


#endif
