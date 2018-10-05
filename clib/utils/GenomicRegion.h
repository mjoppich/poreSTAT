#ifndef GENOMICREGION_H_
#define GENOMICREGION_H_

/**
    Splitter Application Package - some toolkit to work with gtf/gff files
    Copyright (C) 2015  Markus Joppich

    The Splitter Application Package is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The Splitter Application Package is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <inttypes.h>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <algorithm>
#include <functional>
#include "Utils.h"
#include "StringUtils.h"

/**
 *
 * \brief A class to represent a genomic region
 *
 * This class implements a GenomicRegion (1-based). It is useful for any sort of range/region on a genome/transcriptome
 *
 * A region from [10 15] contains 6 bases (namely 10,11,12,13,14,15)
 *
 * \author Markus Joppich
 * \date 2015.11.18
 *
 */

class GenomicRegion
{
public:


    /**
     *
     * @param iStart
     * @param iEnd
     * @param iChromID
     * @param iStrandInfo Remember: Strand Info. 1: 5'->3', 0: 3'->5' -1/255: no info
     * @return GenomicRegion starting at iStart, ending at Iend on Chromsosome iChromID and strand iStrandInfo
     */
    GenomicRegion(uint32_t iStart, uint32_t iEnd, uint32_t iChromID, uint8_t iStrandInfo)
            : GenomicRegion(iStart, iEnd, iChromID)
    {
        this->setStrandInfo(iStrandInfo);
    }

    /**
     *
     * @param iStart
     * @param iEnd
     * @param iChromID
     * @return Genomic Region starting at iStart, ending at iEnd and on Chromosome iChromID
     */
    GenomicRegion(uint32_t iStart, uint32_t iEnd, uint32_t iChromID)
        : GenomicRegion(iStart, iEnd)
    {
        m_iChromID = iChromID;
    }

    /**
     *
     * @param iStart
     * @param iEnd
     * @return GenomicRegion starting at iStart and ending at iEnd - no chromosome info nor strand info is known
     */
    GenomicRegion(uint32_t iStart, uint32_t iEnd)
    {
        m_iStart = iStart;
        m_iEnd = iEnd;
    }

    /**
     *
     * @param pOther
     * @return copy constructor
     */
    GenomicRegion(GenomicRegion* pOther)
    {
        const GenomicRegion& oOther = *pOther;

        this->update(oOther);
    }

    virtual ~GenomicRegion() {


    }

    /**
     *
     * @param pOther deferred copy constructor from other GenomicRegion
     */
    void update(const GenomicRegion& oOther)
    {

        m_iStart = oOther.getStart();
        m_iEnd = oOther.getEnd();
        m_iChromID = oOther.getChromosomeID();
        m_iStrandInfo = oOther.getStrandInfo();

    }


    uint32_t getStart() const
    {
        return this->m_iStart;
    }

    uint32_t getEnd() const
    {
        return this->m_iEnd;
    }

    uint32_t getLength() const
    {

        return this->m_iEnd - this->m_iStart + 1;

    }

    void setStart(int iStart)
    {
        this->m_iStart = iStart;
    }
    void setEnd(int iEnd)
    {
        this->m_iEnd = iEnd;
    }

    uint8_t charToStrand(char cStrandInfo) {
        uint8_t iStrandInfo;
        if (cStrandInfo == '.')
            iStrandInfo = -1;
        else
            iStrandInfo = (cStrandInfo == '+') ? 1 : 0;

        return iStrandInfo;
    }

    char getStrand() const
    {
        if (m_iStrandInfo == -1)
            return '.';

        return (m_iStrandInfo == 1) ? '+' : '-';
    }

    uint32_t getChromosomeID() const
    {
        return m_iChromID;
    }

    void setStrand(char cStrandInfo)
    {
        m_iStrandInfo = this->charToStrand(cStrandInfo);
    }

    void setStrandInfo(uint8_t iStrandInfo)
    {
        m_iStrandInfo = iStrandInfo;
    }

    /**
     *
     * @return Strand Info. 1: 5'->3', 0: 3'->5' -1/255: no info
     */
    uint8_t getStrandInfo() const
    {
        return m_iStrandInfo;
    }

    uint8_t getReversedStrandInfo()
    {
        uint8_t iStrandInfo = -1;
        if (m_iStrandInfo == 1)
        {
            iStrandInfo = 0;
        }
        else {
            if (m_iStrandInfo == 0)
                iStrandInfo = 1;
        }

        return iStrandInfo;

    }

    /**
     * \brief checks whether a position is within this region
     * \param iPosition position to check
     * \return true if [start <= iPosition <= end]
     */
    bool contains(uint32_t iPosition)
    {
        return ((m_iStart <= iPosition) && (iPosition <= m_iEnd));
    }

    /**
     *
     * \brief checks whether a region is inside this region
     * \param oOther region to check for
     * \return true if [start <= other.start <= other.end <= end]
     */
    bool contains(GenomicRegion& oOther)
    {

        return this->contains(&oOther);
    }

    /**
     *
     * \brief checks whether a region is inside this region (or if bPartially==true overlaps at borders)
     * \param oOther region to check for
     * \param bPartially also allows [start <= other.start <= end <= other.end] or [other.start <= start <= other.end <= end]
     * \return true if [start <= other.start <= other.end <= end]
     */
    bool contains(GenomicRegion* pOther, bool bPartially = false)
    {

        bool bStartContained = this->contains(pOther->m_iStart);
        bool bEndContained = this->contains(pOther->m_iEnd);

        if (!bPartially)
            return bStartContained && bEndContained;

        return bStartContained || bEndContained;
    }

    /**
     *
     * \brief checks whether a region overlaps with another region
     * \param oOther region to check overlap for
     * \return true if this includes other region [start, o.start, o.end, end] or [o.start,start,end,o.end] or [start,o.start,end,o.end] or [o.start,start,o.end,end]
     */
    bool overlaps(GenomicRegion* pOther)
    {

        // this includes other
        if ((this->getStart() <= pOther->getStart()) && (pOther->getEnd() <= this->getEnd()))
            return true;

        // this is included by other
        if ((pOther->getStart() <= this->getStart()) && (this->getEnd() <= pOther->getEnd()))
            return true;

        if ( (this->getStart() <= pOther->getStart()) && (this->getEnd() <= pOther->getEnd()) && (this->getEnd() >= pOther->getStart()) )
            return true;

        if ( (pOther->getStart() <= this->getStart()) && (pOther->getEnd() <= this->getEnd()) && (this->getStart() <= pOther->getEnd()) )
            return true;

        return false;
    }


    enum REL_POSITION {

        LEFT=0, RIGHT, IN5P, IN3P,
        OVERLAP5P, OVERLAP3P, OVERLAPLEFT, OVERLAPRIGHT,
        EQUAL,
        CONTAINS, // [this [other]]
        CONTAINED, // [other [this]]
        ISOLATED, //[this] .... [other]
        ERROR
    };



    /**
     * \param oOther other element
     * \return this is .... compared to pOther
     */
    REL_POSITION relateTo(GenomicRegion* pOther)
    {

        if ((this->getStart() == pOther->getStart()) && (this->getEnd() == pOther->getEnd() ))
            return REL_POSITION::EQUAL;

        if (this->contains(pOther))
            return REL_POSITION::CONTAINS;

        if (pOther->contains(this))
            return REL_POSITION::CONTAINED;

        bool bStranded = (this->getStrandInfo() != (uint8_t)-1);

        if (pOther->getStart() < this->getStart())
        {

            if (pOther->getEnd() < this->getStart())
            {

                if (bStranded)
                {

                    if (this->getStrandInfo() == 1)
                    {
                        return REL_POSITION::IN3P;
                    } else {
                        return REL_POSITION::IN5P;
                    }

                }

                return REL_POSITION::RIGHT;
            } else {

                if (bStranded)
                {

                    if (this->getStrandInfo() == 1)
                    {
                        return REL_POSITION::OVERLAP3P;
                    } else {
                        return REL_POSITION::OVERLAP5P;
                    }

                }

                return REL_POSITION::OVERLAPRIGHT;
            }

        } else {

            if (this->getEnd() < pOther->getStart()) {


                if (bStranded)
                {

                    if (this->getStrandInfo() == 1)
                    {
                        return REL_POSITION::IN5P;
                    } else {
                        return REL_POSITION::IN3P;
                    }

                }

                return REL_POSITION::LEFT;
            } else {

                if (bStranded)
                {

                    if (this->getStrandInfo() == 1)
                    {
                        return REL_POSITION::OVERLAP5P;
                    } else {
                        return REL_POSITION::OVERLAP3P;
                    }

                }

                return REL_POSITION::OVERLAPLEFT;
            }

        }

        return REL_POSITION::ERROR;

    }


    /**
     * \brief returns string description: [start end] cDel length
     *
     * */
    virtual std::string toString(char cDel = ' ')
    {

        std::stringstream oStringStream;
        oStringStream << this->getStart() << cDel << this->getEnd() << cDel << this->getLength();
        oStringStream.flush();

        return oStringStream.str();

    }


/**
 * \brief calculates all non covered regions in pNegElemnets within [start, end] of this region
 *
 * for region [100, 500] and pNegElements={[200,300], [400,450]} this would return {[100,199],[301,399],[451,500]}
 *
 * \param pNegElements Covered Elements within this region
 * \return vector of genomicregion* filling any gaps in pNegElements
 */
    std::vector<GenomicRegion*>* fillEmptyRegions( std::vector<GenomicRegion*>* pSortedNegElements ) {
        if ((pSortedNegElements == NULL) || (pSortedNegElements->size() < 1)) {
            return NULL;
        }

        std::vector<GenomicRegion *> *pReturn = new std::vector<GenomicRegion *>();

        GenomicRegion *pOldElement = pSortedNegElements->at(0);

        if (this->getStart() < pOldElement->getStart()) {
            // add new intron
            GenomicRegion *pGap = new GenomicRegion(this->getStart(),pOldElement->getStart() - 1);
            pGap->setChromosomeID(this->getChromosomeID());
            pGap->setStrandInfo(this->getStrandInfo());

            if (pGap->getStart() <= pGap->getEnd())
                pReturn->push_back(pGap);
        }

        uint32_t iNext = 1;
        GenomicRegion *pElement = NULL;
        while (iNext < pSortedNegElements->size()) {
            pElement = pSortedNegElements->at(iNext);
            ++iNext;

            uint32_t iDist = pElement->getStart() - pOldElement->getEnd();

            if (iDist > 1) {
                GenomicRegion *pGap = new GenomicRegion(pOldElement->getEnd() + 1, pElement->getStart() - 1);
                pGap->setChromosomeID(this->getChromosomeID());
                pGap->setStrandInfo(this->getStrandInfo());

                if (pGap->getStart() <= pGap->getEnd())
                    pReturn->push_back(pGap);
            }

            pOldElement = pElement;

        }

        if (pOldElement->getEnd() < this->getEnd()) {
            GenomicRegion *pGap = new GenomicRegion(pOldElement->getEnd() + 1, this->getEnd());
            pGap->setChromosomeID(this->getChromosomeID());
            pGap->setStrandInfo(this->getStrandInfo());

            if (pGap->getStart() <= pGap->getEnd())
                pReturn->push_back(pGap);
        }

        return pReturn;
    }





    void setChromosomeID(uint32_t iID)
    {
        m_iChromID = iID;
    }

    /**
     * \brief tests for euqality with other GenomicRegion
     * \param pOther other GenomicRegion
     * \return true if start == o.start, end == o.end and if ChromID is set, chromid == o.chromid
     *
     */
    bool equals(GenomicRegion* pOther)
    {

        bool bEquals = true;
        if (this->m_iChromID != (uint32_t) -1)
            bEquals = (this->m_iChromID == pOther->m_iChromID);

        if (this->m_iStrandInfo != (uint32_t) -1)
            bEquals = bEquals && (this->m_iStrandInfo == pOther->m_iStrandInfo);

        return bEquals && (this->getStart() == pOther->getStart()) && (this->getLength() == pOther->getLength());
    }

    /**
     * \see equals(GenomicRegion* pOther)
     */
    bool operator==(const GenomicRegion& oOther)
    {

        return this->equals((GenomicRegion*)&oOther);

       if ( this->m_iStart == oOther.m_iStart)
           if (this->m_iEnd == oOther.m_iStart)
               if (this->m_iChromID == oOther.m_iChromID)
                   return true;

        return false;

    }

    /**
     * \see operator==(const GenomicRegion& oOther)
     * */
    bool operator==(const GenomicRegion* pOther)
    {
        return this->equals((GenomicRegion*) pOther);

    }


    void copyAttributes(GenomicRegion* pOther)
    {

        this->m_mAttributes.insert( pOther->m_mAttributes.begin(), pOther->m_mAttributes.end() );

    }

    // merge of attributes with gffentry
    bool hasAttribute(std::string sAttribute)
    {
        std::map<std::string, std::string>::iterator oIt = m_mAttributes.find( sAttribute );

        if (oIt == m_mAttributes.end())
            return false;

        return true;
    }

    bool hasAttribute(std::string sAttribute, std::string sValue)
    {
        std::map<std::string, std::string>::iterator oIt = m_mAttributes.find( sAttribute );

        if (oIt == m_mAttributes.end())
            return false;

        return (oIt->second.compare(sValue) == 0);
    }

    std::string getAttributeByValue(std::string sAttribute)
    {

        std::string* pAttrib = this->getAttribute(sAttribute);

        if (pAttrib == NULL)
            return std::string("");

        std::string sReturn = *pAttrib;
        delete pAttrib;
        return sReturn;

    }


    bool hasAttribute(std::string sAttribute, char cValue)
    {
        std::string sOneElem = "A";
        sOneElem[0] = cValue;

        return this->hasAttribute(sAttribute, sOneElem);
    }

    void setAttribute(std::string sKey, char cValue)
    {

        std::string sOneElem = "A";
        sOneElem[0] = cValue;

        this->setAttribute(sKey, sOneElem);

    }

    std::string* getAttribute(std::string sKey) {

        std::transform(sKey.begin(), sKey.end(), sKey.begin(), ::toupper);

        std::map<std::string, std::string>::iterator oIt = m_mAttributes.find(sKey);

        if (oIt != m_mAttributes.end())
            return new std::string(oIt->second);

        return NULL;
    }

    void printAttribute(std::string* pAttrib = NULL) {

        if (pAttrib == NULL)
        {
            std::map<std::string, std::string>::iterator oIt = m_mAttributes.begin();

            for ( ; oIt != m_mAttributes.end(); ++oIt)
            {
                std::cout << oIt->first << '\t' << oIt->second << std::endl;
            }


            return;
        }


        std::transform(pAttrib->begin(), pAttrib->end(), pAttrib->begin(), ::toupper);

        std::map<std::string, std::string>::iterator oIt = m_mAttributes.find(*pAttrib);

        if (oIt != m_mAttributes.end())
        {
            std::cout << oIt->first << '\t' << oIt->second << std::endl;
        }


        return;
    }

    void setAttribute(std::string* pKey, std::string* pValue) {
        std::transform(pKey->begin(), pKey->end(), pKey->begin(), ::toupper);

        std::pair < std::map<std::string, std::string>::iterator, bool> oReturn = m_mAttributes.insert(std::pair<std::string, std::string>(*pKey, *pValue));

        if (oReturn.second == false)
        {
            // replace existing element
            oReturn.first->second = *pValue;
        }

    }

    void setAttribute(std::string sKey, std::string sValue)
    {
        this->setAttribute(&sKey, &sValue);
    }


    std::vector<std::string>* getAllAttributes()
    {

        std::map<std::string, std::string>::iterator oIt = m_mAttributes.begin();

        std::vector<std::string>* pAllAttributes = new std::vector<std::string>();

        for ( ; oIt != m_mAttributes.end(); ++oIt)
            pAllAttributes->push_back( oIt->first );

        return pAllAttributes;


    }


    static GenomicRegion *parseFromString(std::string sLine, char cDel = '\t')
    {

        std::vector<std::string> vElems = StringUtils::split(&sLine, cDel);

        uint32_t iChromID = -1;
        uint32_t iStart = std::stoi(vElems[1]);
        uint32_t iEnd = std::stoi(vElems[2]);

        try {

            iChromID = std::stoi(vElems[0]);
        }
        catch (const std::invalid_argument& ia) {
            iChromID = -1;
        }

        GenomicRegion* pReturn = new GenomicRegion(iStart, iEnd, iChromID);

        if (iChromID == (uint32_t) -1)
            pReturn->setAttribute("seqname", vElems[0]);

        pReturn->setStrand( vElems[3].at(0) );

        return pReturn;

    }


protected:

    /**
     * \brief Start of Region
     */
    uint32_t m_iStart = -1;
    /**
     * \brief End of Region
     */
    uint32_t m_iEnd = -1;
    /**
     * \brief Chromosome ID
     */
    uint32_t m_iChromID = -1;

    /**
     * \brief Strand Info. 1: 5'->3', 0: 3'->5' -1/255: no info
     */
    uint8_t m_iStrandInfo = -1;

    /** \brief Attributes Map
     *
     * Attributes Map holding attributes as upper-cased Element -> Value
     * */
    std::map<std::string, std::string> m_mAttributes;


};


#endif /*GENOMICREGION_H_*/
