#ifndef TSFTIMEMONITOR_H
#define TSFTIMEMONITOR_H

#include "TSFDefs.h"
#include "TSFTimer.h"

namespace TSF
{

  using std::string;

  /** \ingroup IO
   *
   */
  class TSFTimeMonitor
    {
    public:
      /** */
      TSFTimeMonitor(TSFTimer& timer)
        : timer_(timer), isRoot_(!timer.isRunning())
        {
          if (isRoot_) timer_.start();
        }

      /** */
      inline ~TSFTimeMonitor()
        {
          if (isRoot_) timer_.stop();
        }
    private:
      TSFTimer& timer_;
      bool isRoot_;
    };

}
#endif
