/*
 * Splitter Application Package - some toolkit to work with gtf/gff files
 * Copyright (C) 2015  Markus Joppich
 *
 * The Splitter Application Package is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The Splitter Application Package is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef PROJECT_FASTAREADER_H
#define PROJECT_FASTAREADER_H

#include <stdio.h>
#include <string>
#include <map>
#include <math.h>
#include <vector>
#include <cmath>
#include "GenomicRegion.h"
#include "GenomicRegionUtils.h"

#include "Sequence.h"
#include "SequenceUtils.h"

#include <htslib/faidx.h>

/**
 * \brief Can open and extract sequences from a FASTA file
 *
 * All position information are 1-based!
 *
 */
class FASTAreader {
public:

    FASTAreader(std::string* pFileName, std::string* pIndexFile)
    {

        m_pFileName = pFileName;
        m_pIndexFile = pIndexFile;
        m_iLineLength = 60;

        m_pIndexToPosition = new std::map<std::string, size_t>();
        m_pIndexToSeqLength = new std::map<std::string, size_t>();
        m_pSequences = new std::map<std::string, std::string>();


        if (pIndexFile == NULL)
        {
            int iFAISuccess = fai_build(pFileName->c_str());

            if (iFAISuccess != 0)
            {
                std::cerr << "Error building FAI index for: " << *pFileName << std::endl;
            }

            m_pIndexFile = new std::string((*pFileName) + ".fai");
        }

        this->parseIndexFile();

    }

    FASTAreader(std::string sFileName, std::string sIndexFile)
        : FASTAreader(new std::string(sFileName), new std::string(sIndexFile))
    {
    }

    ~FASTAreader()
    {

        SAFEDEL(m_pIndexToSeqLength)
        SAFEDEL(m_pIndexToPosition)
        SAFEDEL(m_pSequences)

    }

    void loadAllSequences()
    {
        std::vector<std::string>* pLines = Utils::readByLine( m_pIndexFile );

        std::cout << "Fasta Index File" << pLines->size() << std::endl;

        for (uint32_t i = 0; i < pLines->size(); ++i)
        {

            std::vector<std::string> vLineContent = StringUtils::split(pLines->at(i), '\t');

            if (vLineContent.size() != 5)
                continue;

            std::string sIndex = vLineContent.at(0);
            size_t iSeqStart = std::stoul(vLineContent[2]);
            size_t iSeqLength = std::stoul(vLineContent[1]);
            
            std::string sSeq = this->retrieveSequence(iSeqStart,0, iSeqLength);

            std::pair<std::string, std::string> oSeq( sIndex, sSeq );
            m_pSequences->insert(oSeq);



        }

    }

    size_t getSequenceLength(const std::string* pIndex);
    std::vector<std::string>* getSequenceIndeces()
    {

        std::vector<std::string>* pReturn = new std::vector<std::string>();

        std::map<std::string, size_t>::iterator oIt = m_pIndexToSeqLength->begin();

        for (; oIt != m_pIndexToSeqLength->end(); ++oIt )
        {

            pReturn->push_back( oIt->first );

        }

        return pReturn;

    }

    std::string getFileName()
    {
        return std::string(m_pFileName->c_str());
    }

    std::string getIndexFileName()
    {
        return std::string(m_pIndexFile->c_str());
    }

    size_t getPosition(const std::string* pIndex);
    size_t getPosition(uint32_t iIndex);

    void setLineLength(uint32_t iLength)
    {
        m_iLineLength = iLength;
    }

    size_t removeNewLines(char* pBuffer, size_t iBufferSize);

    std::string * readFromFile(size_t iPositionInFile, size_t iChars);

    void parseIndexFile();

    size_t getSequenceStart(uint32_t iPosition);


    /*
    std::string * retrieveSequence(GenomicRegion* pRegion, XAMReader* pReader = NULL)
    {

        size_t iFileStart = -1;

        if ( pReader != NULL) {

            std::string sSeqName = pReader->getSeqName( pRegion->getChromosomeID() );
            iFileStart = this->getPosition(&sSeqName);

        } else {

            iFileStart = this->getPosition( pRegion->getChromosomeID() );
        }

        if (iFileStart == -1)
            return NULL;


        return this->retrieveSequence(iFileStart, pRegion->getStart(), pRegion->getLength());

    }
    */

    std::string retrieveSequence(GenomicRegion* pRegion, std::string& sChromIdentifier);

    std::string retrieveSequences(std::vector<GenomicRegion*>* pRegions, std::string& sChromIdentifier)
    {
        std::sort( pRegions->begin(), pRegions->end(), GenomicRegionUtils::sortStartAsc );

        std::string sSequence = "";

        for (size_t i = 0; i < pRegions->size(); ++i)
        {

            std::string sRegionSeq = this->retrieveSequence(pRegions->at(i), sChromIdentifier);

            sSequence.append(sRegionSeq);
        }

        return sSequence;

    }

    std::string retrieveSequence(size_t iFileStart, uint32_t iSequenceStart, uint32_t iSequenceLength )
    {


        if (iFileStart == (size_t) -1)
            return "";

        size_t iSeqStart = this->getSequenceStart( iSequenceStart );

        std::string* pReturn = this->readFromFile(iFileStart + iSeqStart, iSequenceLength);

        std::string oResult( pReturn->substr(0, iSequenceLength) );

        delete pReturn;

        return oResult;

    }


    /**
     *
     * \brief retrieves a Sequence object from Start with length |iLength|. If iLength < 0, [start-ilength, start].
     * \param iStart start base
     * \param iLength length to retrieve
     * \param pIntrons any intronic sequences that should be skipped. Introns may not overlap!
     */
    Sequence retrieveSequence(std::string& sChromIdentifier, uint32_t iStart, int32_t iLength, std::vector<GenomicRegion*>* pIntrons = NULL)
    {

        if (pIntrons != NULL)
        {
            if (iLength < 0)
            {

                std::sort(pIntrons->begin(), pIntrons->end(), GenomicRegionUtils::sortStartAsc);

            } else {

                std::sort(pIntrons->begin(), pIntrons->end(), GenomicRegionUtils::sortEndDesc);

            }
        }


        size_t iCurrentIntron = 0;
//        uint32_t iSeqPosition = iStart;

        uint32_t iMaxLength = std::abs(iLength) + GenomicRegionUtils::getLengthVec(pIntrons, 0);



        size_t iFileStart = this->getPosition( &sChromIdentifier );

        Sequence oSequence("");
        size_t iIdx = 0;

        if (iLength > 0)
        {
            std::string sChromSequence = this->retrieveSequence(iFileStart, iStart, iMaxLength);

            while ((iIdx < sChromSequence.size()) && (oSequence.size() < (size_t) std::abs(iLength)))
            {

                if ((pIntrons != NULL) && (pIntrons->size() > iCurrentIntron) && (iIdx + iStart == pIntrons->at(iCurrentIntron)->getStart() ))
                {

                    iIdx = iIdx + pIntrons->at(iCurrentIntron)->getLength();
                    ++iCurrentIntron;

                }


                oSequence.push_back( sChromSequence[iIdx] );

                ++iIdx;
            }

        } else {

            std::string sChromSequence = this->retrieveSequence(iFileStart, iStart - iMaxLength+1, iMaxLength);
            Sequence oSeq(sChromSequence);
            oSeq.reverse();

            while ((iIdx < sChromSequence.size()) && (oSequence.size() < (size_t) std::abs(iLength)))
            {

                if ((pIntrons != NULL) && (pIntrons->size() > iCurrentIntron) && (iStart - iIdx == pIntrons->at(iCurrentIntron)->getEnd() ))
                {

                    iIdx = iIdx - pIntrons->at(iCurrentIntron)->getLength();
                    ++iCurrentIntron;

                }


                oSequence.push_back( sChromSequence[iIdx] );

                ++iIdx;
            }

        }

        return oSequence;

    }


protected:

    std::string* m_pFileName;
    std::string* m_pIndexFile;
    uint32_t m_iLineLength;

    std::map<std::string, size_t>* m_pIndexToPosition;
    std::map<std::string, size_t>* m_pIndexToSeqLength;

    std::map<std::string, std::string>* m_pSequences;


};


#endif //PROJECT_FASTAREADER_H
