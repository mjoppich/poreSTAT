#include <iostream>
#include <vector>
#include "../utils/FASTAreader.h"
#include "ReadStats.h"


// A simple class with a constuctor and some methods...
class AlignmentStatistics
{
    public:
        AlignmentStatistics();

        void processFiles(int n_samples, char** pFilenames);

        void loadFASTA(char* pFastaFile);

        std::vector<PythonReadStats>* readStats;

    protected:

        FASTAreader* pReader = NULL;
        std::string* m_pFASTAFile = NULL;
};





