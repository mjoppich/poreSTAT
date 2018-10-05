//
// Created by joppich on 11/27/15.
//

#ifndef PROJECT_READREGION_H
#define PROJECT_READREGION_H

#include "GenomicRegion.h"

class ReadRegion : public GenomicRegion {

public:

    ReadRegion(int iStart, int iEnd, char cCIGAR, uint32_t iCIGARLength)
            : GenomicRegion(iStart, iEnd)
    {
        this->m_cCIGAR = cCIGAR;
        this->m_iCIGARLength = iCIGARLength;
    }

    ~ReadRegion()
    {
        if (this->pReadRegion != NULL)
        {
            delete this->pReadRegion;
        }
    }

    char getCIGAR()
    {
        return m_cCIGAR;
    }

    void setCIGAR(char cCIGAR)
    {
        this->m_cCIGAR = cCIGAR;
    }

    uint32_t getCIGARLength()
    {
        return m_iCIGARLength;
    }

    void setCIGARLenght(uint32_t iCIGARLength)
    {
        this->m_iCIGARLength = iCIGARLength;
    }


    void setSequence(std::string sSequence)
    {
        this->m_sSequence = sSequence;
    }

    std::string getSequence()
    {
        return m_sSequence;
    }


    void setQueryName(std::string sQueryName)
    {
        this->m_sQueryName = sQueryName;
    }

    std::string getQueryName()
    {
        return m_sQueryName;
    }


    void setDebug(std::string sDebug)
    {
        this->m_sDebug = sDebug;
    }

    std::string getDebug()
    {
        return m_sDebug;
    }

    GenomicRegion* pReadRegion = NULL;

protected:

    char m_cCIGAR = 0;
    uint32_t m_iCIGARLength = 0;
    std::string m_sSequence = "";
    std::string m_sQueryName = "";
    std::string m_sDebug = "";

};


#endif //PROJECT_READREGION_H
