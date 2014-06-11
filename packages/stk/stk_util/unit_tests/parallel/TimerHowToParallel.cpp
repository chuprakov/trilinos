#include <gtest/gtest.h>
#include <stk_util/diag/PrintTimer.hpp>
#include <stk_util/diag/Timer.hpp>
#include <comparison/stringAndNumberComparisons.h>

namespace
{

const double tolerance = 5e-2;

void doWork()
{
    ::usleep(1e5);
}

TEST(StkDiagTimerHowTo, useTimersInParallel)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);
    if(numProcs == 2)
    {
        enum {CHILDMASK1 = 1};
        stk::diag::TimerSet enabledTimerSet(CHILDMASK1);
        stk::diag::Timer rootTimer = createRootTimer("totalTestRuntime", enabledTimerSet);
        rootTimer.start();

        stk::diag::Timer childTimer1("childTimer1", CHILDMASK1, rootTimer);

        {
            stk::diag::TimeBlockSynchronized timerStartSynchronizedAcrossProcessors(childTimer1, communicator);
            doWork();
        }

        std::ostringstream outputStream;
        bool printTimingsOnlySinceLastPrint = false;
        stk::diag::printTimersTable(outputStream, rootTimer, stk::diag::METRICS_ALL, printTimingsOnlySinceLastPrint, communicator);

        int procId = -1;
        MPI_Comm_rank(communicator, &procId);
        if(procId == 0)
        {
            std::string expectedOutput = "                                                  \
                                CPU Time            CPU Time            CPU Time            \
                               Wall Time           Wall Time           Wall Time            \
     Timer      Count      Sum (% of System)   Min (% of System)   Max (% of System)        \
                           Sum (% of System)   Min (% of System)   Max (% of System)        \
SKIP ----- --------------------- --------------------- SKIP SKIP SKIP --------------------- \
totalTestRuntime  2            0.001 SKIP          0.000 SKIP          0.001 SKIP           \
                               0.200 SKIP          0.100 SKIP          0.100 SKIP           \
  childTimer1     2            0.001 SKIP          0.000 SKIP          0.001 SKIP           \
                               0.200 SKIP          0.100 SKIP          0.100 SKIP           \
                    ";
            EXPECT_TRUE(unitTestUtils::areStringsEqualWithToleranceForNumbers(expectedOutput, outputStream.str(), tolerance));
        }

        stk::diag::deleteRootTimer(rootTimer);
    }
}

}
