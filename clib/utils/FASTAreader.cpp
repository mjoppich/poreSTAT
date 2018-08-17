//
// Created by joppich on 09.11.15.
//
#include "Utils.h"

#include "FASTAreader.h"

size_t FASTAreader::getSequenceLength(const std::string *pIndex) {

    std::map<std::string, size_t>::iterator oIt = m_pIndexToSeqLength->find( *pIndex );

    if (oIt == m_pIndexToSeqLength->end())
        return -1;

    return oIt->second;

}

size_t FASTAreader::getPosition(const std::string* pIndex) {

    std::map<std::string, size_t>::iterator oIt = m_pIndexToPosition->find( *pIndex );

    if (oIt == m_pIndexToPosition->end())
        return -1;

    return oIt->second;

}

size_t FASTAreader::getPosition(uint32_t iIndex) {

    std::map<std::string, size_t>::iterator oIt = m_pIndexToPosition->find( std::to_string(iIndex) );

    if (oIt == m_pIndexToPosition->end())
        return -1;

    return oIt->second;

}

size_t FASTAreader::removeNewLines(char* pBuffer, size_t iBufferSize) {

    size_t iPosition = 0;
    size_t i = 0;
    for (; (i < iBufferSize) && (iPosition < iBufferSize); ++i)
    {

        while ((pBuffer[iPosition] < 32) && (iPosition < iBufferSize))
        {
            ++iPosition;
        }

        if (!(iPosition < iBufferSize))
            break;

        pBuffer[i] = pBuffer[iPosition];
        ++iPosition;

    }

    size_t iAdded = i;
    for (; i < iBufferSize; ++i)
        pBuffer[i] = 0;

    return iAdded;

}

void FASTAreader::parseIndexFile() {
    std::vector<std::string>* pLines = Utils::readByLine( m_pIndexFile );

    uint32_t iLineLength = 0;

    for (uint32_t i = 0; i < pLines->size(); ++i)
    {

        std::vector<std::string> vLineContent = StringUtils::split(pLines->at(i), '\t');

        if (vLineContent.size() != 5)
            continue;

        std::string sIndex = vLineContent.at(0);
        size_t iSeqStart = std::stoul(vLineContent[2]);
        size_t iSeqLength = std::stoul(vLineContent[1]);
        iLineLength = std::stoul(vLineContent[3]);

        std::pair<std::string, size_t> oSeqStart( sIndex, iSeqStart );
        m_pIndexToPosition->insert(oSeqStart);

        std::pair<std::string, size_t> oSeqLength(sIndex, iSeqLength);
        m_pIndexToSeqLength->insert(oSeqLength);

    }

    this->setLineLength(iLineLength);

}

std::string *FASTAreader::readFromFile(size_t iPositionInFile, size_t iChars) {

    size_t iExtraChars = std::ceil((float)iChars / (float) m_iLineLength);
//    size_t iOriginalChars = iChars;
    iChars = iChars + iExtraChars;

    const size_t iBufferSize = 256;
    char* aBuffer = (char*) malloc(sizeof(char) * iBufferSize);

    FILE * pFile;
    pFile = fopen ( m_pFileName->c_str() , "r" );
    std::string* pReturn = new std::string();
    fseek ( pFile , iPositionInFile , SEEK_SET );


    for (size_t i = 0; i < iChars; i += iBufferSize)
    {

        fread(aBuffer, sizeof(char), iBufferSize, pFile);

        size_t iToCopy = this->removeNewLines(aBuffer, iBufferSize);

        /* check for buffer overruns
        for (uint32_t j = 0; j < iToCopy; ++j)
        {
            if (!((aBuffer[j] == 65) || (aBuffer[j] == 84) || (aBuffer[j] == 67) || (aBuffer[j] == 71)))
            {
                std::cerr << "ERROR (invalid seq) Fasta reader|" << (uint8_t) aBuffer[j] << "|" << std::endl;
            }
        }
        */

        pReturn->append(aBuffer,iToCopy);
    }

    free(aBuffer);
    fclose ( pFile );
    return pReturn;

}

std::string FASTAreader::retrieveSequence(GenomicRegion* pRegion, std::string& sChromIdentifier) {

    size_t iFileStart = this->getPosition( &sChromIdentifier );

    if (iFileStart == (size_t) -1)
        return "";

    return this->retrieveSequence(iFileStart, pRegion->getStart(), pRegion->getLength());
}

/**
 * \brief translates position to character-position in fasta file
 *
 * \param iPosition genomic! sequence start
 * \return sequence start in fasta file
 *
 */
size_t FASTAreader::getSequenceStart(uint32_t iPosition) {

    if (iPosition == 0)
    {
        std::cerr << "input is not one-based genomic position!" << std::endl;
        return 0;
    }

    iPosition = iPosition -1;

    size_t iLines = iPosition / m_iLineLength;

    size_t iSeqStart = iLines * (m_iLineLength+1) + iPosition % m_iLineLength;

    return iSeqStart;

}
