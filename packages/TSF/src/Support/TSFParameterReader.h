#ifndef TSFPARAMETERREADER_H
#define TSFPARAMETERREADER_H

#include "TSFDefs.h"
#include <string>
#include "TSFArray.h"

namespace TSF
{

  using namespace std;

  class TSFParameterReader
    {
    public:
      TSFParameterReader(const TSFSmartPtr<TSFParameterReaderBase>& ptr);

      TSFParameterList read() const ;

    private:
      TSFSmartPtr<TSFParameterReader> ptr_;
    };

}

#endif
