// Define C functions for the C++ class - as ctypes can only talk to C...

#include "./align_stats/AlignmentStatistic.h"
#include <vector>
#include <iostream>

extern "C"
{
    AlignmentStatistics* AlignmentStatistics_new()
    {
        return new AlignmentStatistics();
    }
    void AlignmentStatistics_process(AlignmentStatistics* astat, int ncount, char** pfiles=NULL)
    {
        return astat->processFiles(ncount, pfiles);
    };

    void AlignmentStatistics_load_fasta(AlignmentStatistics* astat, char* pFastaFile)
    {
        std::cout << pFastaFile << std::endl;
        return astat->loadFASTA(pFastaFile);
    };



    void* AlignmentStatistics_ReadStats(AlignmentStatistics* astat)
    {
        return (void*)astat->readStats->data();
    };

    uint32_t AlignmentStatistics_ReadStats_size(AlignmentStatistics* astat)
    {
        std::cout << astat->readStats->size() << std::endl;
        return astat->readStats->size();
    };
}