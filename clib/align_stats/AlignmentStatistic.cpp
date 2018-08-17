#include "./AlignmentStatistic.h"

#include "AlignmentStatisticProcessor.h"

AlignmentStatistics::AlignmentStatistics()
{
}

void AlignmentStatistics::processFiles(int n_samples, char** pFilenames)
{
    std::cout << "Processing files " << n_samples << " " << pFilenames << std::endl;

    for (int i = 0; i < n_samples; ++i)
    {

        std::cout << i << "\t" << pFilenames[i] << std::endl;

        std::string sFileName(pFilenames[i]);
        AlignmentStatisticProcessor* pReader = new AlignmentStatisticProcessor(&sFileName, NULL);

        std::vector<std::string>* pSeqNames = pReader->getSeqNames();

        for (size_t si = 0; si < pSeqNames->size(); ++si)
        {
            std::cout << pSeqNames->at(si) << std::endl;
        }

        pReader->startAllParallel(10, NULL);

        std::cout << pReader->seenReads << std::endl;
    }

    
}