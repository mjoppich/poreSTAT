// Define C functions for the C++ class - as ctypes can only talk to C...

#include "./align_stats/AlignmentStatistic.h"



extern "C"
{
    AlignmentStatistics* AlignmentStatistics_new() {return new AlignmentStatistics();}
    void AlignmentStatistics_process(AlignmentStatistics* astat, int ncount, char** pfiles=NULL) {return astat->processFiles(ncount, pfiles);};
}