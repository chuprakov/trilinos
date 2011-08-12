#ifndef CTHULHU_LOOKUPSTATUS_HPP
#define CTHULHU_LOOKUPSTATUS_HPP

#include "Cthulhu_ConfigDefs.hpp"
#ifdef HAVE_CTHULHU_TPETRA  
#include "Tpetra_ConfigDefs.hpp"
#endif

namespace Cthulhu {

#ifdef HAVE_CTHULHU_TPETRA  
  const Tpetra::LookupStatus toTpetra(Cthulhu::LookupStatus); //unused
  const Cthulhu::LookupStatus toCthulhu(Tpetra::LookupStatus);

  const Tpetra::ProfileType toTpetra(Cthulhu::ProfileType);
  const Tpetra::OptimizeOption toTpetra(Cthulhu::OptimizeOption);
#endif // HAVE_CTHULHU_TPETRA

#ifdef HAVE_CTHULHU_EPETRA  
  // const int toEpetra(Cthulhu::LookupStatus);
  const Cthulhu::LookupStatus toCthulhu(int);
#endif // HAVE_CTHULHU_TPETRA

}

#endif // CTHULHU_LOOKUPSTATUS_HPP
