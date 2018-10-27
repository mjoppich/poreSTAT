#include "./AlignmentStatistic.h"

#include "AlignmentStatisticProcessor.h"

AlignmentStatistics::AlignmentStatistics()
{
    this->readStats = new std::vector<PythonReadStats>();
}

void AlignmentStatistics::loadFASTA(char* pFastaFile)
{
    if (this->pReader != NULL)
    {
        delete this->pReader;
    }

    if (this->m_pFASTAFile != NULL)
    {
        delete this->m_pFASTAFile;
    }

    this->m_pFASTAFile = new std::string(pFastaFile);

    std::cout << "Loading fasta file:" << *(this->m_pFASTAFile) << std::endl;

    this->pReader = new FASTAreader(this->m_pFASTAFile, NULL);
    
}

void AlignmentStatistics::processFiles(int n_samples, char** pFilenames)
{
    std::cout << "Processing files " << n_samples << " " << pFilenames << std::endl;

    for (int i = 0; i < n_samples; ++i)
    {

        std::cout << i << "\t" << pFilenames[i] << std::endl;

        std::string sFileName(pFilenames[i]);
        AlignmentStatisticProcessor* pSAMreader = new AlignmentStatisticProcessor(&sFileName, NULL, this->pReader);
        std::vector<std::string>* pSeqNames = pSAMreader->getSeqNames();

        for (size_t si = 0; si < pSeqNames->size(); ++si)
        {
            std::cout << pSeqNames->at(si) << std::endl;
        }

        this->readStats->clear();
        pSAMreader->startAllParallel(100, this->readStats);


        std::cout << "proc reads" << this->readStats->size() << " " << pSAMreader->seenReads << std::endl;

        std::vector<size_t> vsids = pSAMreader->getSeqIDs();
        for (size_t vi = 0; vi < vsids.size(); ++vi)
        {
            std::cout << "SEQ ID " << vi << " " << pSAMreader->getSeqName(vi);
        }

        delete pSAMreader;

    }

    std::cout << "Processing finished" << std::endl;
    std::cout << std::string(this->readStats->at(0).COVERAGES[0].seq) << std::endl;

}