#include "../utils/XAMReadProcessor.h"

class AlignmentStatisticProcessor : public XAMReadProcessor {
public:


    AlignmentStatisticProcessor(std::string* pBAMFile, std::string* pBAMidxFile)
    : XAMReadProcessor(pBAMFile, pBAMidxFile)
    {


    }


    int seenReads = 0;
    
    
protected:

    virtual void process(bam1_t* pRead, uint32_t iSeqID, void* pData)
    {

        if (this->readIsMapped(pRead))
        {
            this->seenReads += 1;
        }
    }



private:


};