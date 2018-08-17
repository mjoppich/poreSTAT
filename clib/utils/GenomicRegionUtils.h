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

#ifndef PROJECT_GENOMICREGIONUTILS_H
#define PROJECT_GENOMICREGIONUTILS_H

#include "macros.h"
#include "GenomicRegion.h"


/**
 *
 * \brief Utility class for GenomicRegion
 *
 * This class implements some utility functions mainly for vectors of GenomicRegion*
 *
 *
 * \author Markus Joppich
 * \date 2015.11.18
 *
 */
class GenomicRegionUtils {
public:

/**
     * @return true if pRegion1 < pRegion2
     * */
    static bool sortStartAsc(GenomicRegion* pRegion1, GenomicRegion* pRegion2)
    {

        if (pRegion1->getStart() < pRegion2->getStart())
            return true;

        if (pRegion1->getStart() > pRegion2->getStart())
            return false;

        return (pRegion1->getLength() < pRegion2->getLength());
    }

    static bool sortStartDesc(GenomicRegion* pRegion1, GenomicRegion* pRegion2)
    {
        if (pRegion1->getStart() > pRegion2->getStart())
            return true;

        if (pRegion1->getStart() < pRegion2->getStart())
            return false;

        return (pRegion1->getLength() > pRegion2->getLength());

    }

    static bool sortEndDesc(GenomicRegion* pRegion1, GenomicRegion* pRegion2)
    {
        return pRegion1->getEnd() > pRegion2->getEnd();
    }

    static bool sortEndAsc(GenomicRegion* pRegion1, GenomicRegion* pRegion2)
    {
        return pRegion1->getEnd() < pRegion2->getEnd();
    }



    /**
     * \brief returns a new vector of genomic regions containing all regions of pRegions with StrandInfo == iStrandInfo
     *
     * \param pRegions input vector of region
     * \param iStrandInfo strand information all output regions must have. if iStrandInfo > 1, all regions of pRegions are returned
     *
     * \return new vector containing references to regions in pRegions
     *
     */
    static std::vector<GenomicRegion*>* getStrandedRegions(std::vector<GenomicRegion*>* pRegions, uint8_t iStrandInfo)
    {
        std::vector<GenomicRegion*>* pReturn = new std::vector<GenomicRegion*>();
        if (iStrandInfo > 1)
        {
            pReturn->reserve(pRegions->size());
            pReturn->insert(pReturn->end(), pRegions->begin(), pRegions->end());
            return pReturn;
        }

        // approximation
        pReturn->reserve(pRegions->size() / 2);
        for (size_t i = 0; i < pRegions->size(); ++i)
        {
            GenomicRegion* pRegion = pRegions->at(i);
            if (pRegion->getStrandInfo() == iStrandInfo)
                pReturn->push_back( pRegion );

        }

        return pReturn;

    }

    static void deleteVec(std::vector<GenomicRegion*>* pVec)
    {
        for (size_t i = 0; i < pVec->size(); ++i)
        {
            GenomicRegion* pElem = pVec->at(i);
            delete pElem;
        }

        delete pVec;
    }


    static std::vector<GenomicRegion*>* copyVec(std::vector<GenomicRegion*>* pVec)
    {
        if (pVec == NULL)
            return NULL;


        std::vector<GenomicRegion*>* pReturn = new std::vector<GenomicRegion*>();

        pVec->reserve( pVec->size() );

        // TODO make this more efficient
        for (uint32_t i = 0; i < pVec->size(); ++i)
        {
            pReturn->push_back( new GenomicRegion(pVec->at(i)) );
        }

        return pReturn;
    }


    static uint32_t getLengthVec(std::vector<GenomicRegion*>* pVec, uint32_t iNullValue = -1)
    {
        if (pVec == NULL)
            return iNullValue;

        uint32_t iLength = 0;
        const uint32_t iVecLength = pVec->size();
        for (uint32_t i = 0; i < iVecLength; ++i) {

            iLength += pVec->at(i)->getLength();

        }

        return iLength;
    }

    static void deleteVecContent(std::vector<GenomicRegion *> *pElements) {

        if (pElements == NULL)
            return;

        for (uint32_t i = 0; i < pElements->size(); ++i) {

            GenomicRegion *pElem = pElements->at(i);

            delete pElem;

        }

    }

    static void testJoinRegions()
    {

        std::vector<GenomicRegion*> vRegions;
        std::cout << "Before Sort: " << std::endl;
        vRegions.push_back( new GenomicRegion(100,500) );
        vRegions.push_back( new GenomicRegion(150,250) );
        vRegions.push_back( new GenomicRegion(200,300) );
        vRegions.push_back( new GenomicRegion(250,600) );

        vRegions.push_back( new GenomicRegion(800,900) );
        vRegions.push_back( new GenomicRegion(850,950) );

        GenomicRegionUtils::print(&vRegions);

        GenomicRegionUtils::joinRegions(&vRegions);

        std::cout << "After Sort: " << std::endl;

        GenomicRegionUtils::print(&vRegions);

    }


    static bool allContained(GenomicRegion* pRegion, std::vector<GenomicRegion*>* pRegions) {

        bool bContained = true;

        for (uint32_t i = 0; i < pRegions->size(); ++i)
            bContained &= pRegion->contains(pRegions->at(i), true);

        return bContained;
    }

    static std::vector<GenomicRegion*>* joinRegions(std::vector<GenomicRegion*>* pRegions, std::function<void (GenomicRegion*, GenomicRegion*)>* pFunc = NULL)
    {

        if (pRegions == NULL)
            return NULL;

        std::sort(pRegions->begin(), pRegions->end(), GenomicRegionUtils::sortStartAsc);

        std::vector<GenomicRegion*>::iterator oIt, oJt;
        for (oIt = pRegions->begin(); oIt != pRegions->end(); ++oIt)
        {
            oJt = oIt + 1;

            // oIt is last element, no overlap possible
            if (oJt >= pRegions->end())
                break;

            GenomicRegion* pRegion0 = *oIt;
            GenomicRegion* pRegion1 = *oJt;

            while (pRegion0->overlaps(pRegion1)) // rightoverlaps?
            {

                pRegion0->setStart( std::min(pRegion0->getStart(), pRegion1->getStart()) );
                pRegion0->setEnd( std::max(pRegion0->getEnd(), pRegion1->getEnd()));

                delete pRegion1;
                oJt = pRegions->erase(oJt);

                if (oJt >= pRegions->end())
                    break;

                pRegion0 = *oIt;
                pRegion1 = *oJt;

            }

        }
/*
        // REMOVE INNER REGIONS e.g. [100, 200], [125, 175] -> [100, 200]
        std::vector<GenomicRegion*>::iterator oIt, oJt;
        for (oIt = pRegions->begin(); oIt != pRegions->end(); ++oIt)
        {
            oJt = oIt + 1;

            if (oJt == pRegions->end())
                break;

            GenomicRegion* pRegion0 = *oIt;
            GenomicRegion* pRegion1 = *oJt;

            while ( pRegion0->contains(pRegion1) )
            {
                delete pRegion1;
                oJt = pRegions->erase(oJt);

                if (oJt >= pRegions->end())
                    break;

                pRegion0 = *oIt;
                pRegion1 = *oJt;

            }
        }

        // MERGE REGIONS: [100, 200], [150,250] -> [150, 250]
        for (oIt = pRegions->begin(); oIt != pRegions->end(); ++oIt)
        {

            oJt = oIt + 1;

            if (oJt == pRegions->end())
                break;

            GenomicRegion* pRegion0 = *oIt;
            GenomicRegion* pRegion1 = *oJt;

            while ( pRegion0->getEnd() > pRegion1->getStart() )
            {

                if (pFunc != NULL)
                    (*pFunc)(pRegion0, pRegion1);

                pRegion0->setEnd( pRegion1->getEnd() );

                delete pRegion1;
                pRegions->erase(oJt);

                oJt = oIt + 1;

                if (oJt == pRegions->end())
                    break;

                pRegion0 = *oIt;
                pRegion1 = *oJt;

            }

        }
*/
        return pRegions;
    }

    /**
     *
     * @param pReference reference regions (sorted! start asc, not overlapping)
     * @param pRegions regions to check (sorted!, start asc, not overlapping)
     * @return
     */
    static uint32_t getContainedBases(std::vector<GenomicRegion*>* pReference, std::vector<GenomicRegion*>* pRegions)
    {
        // TODO does this work as intended?

        uint32_t iRefMax = GenomicRegionUtils::maxPosition(pReference);

        std::vector<GenomicRegion*>::iterator oRefIt = pReference->begin();
        std::vector<GenomicRegion*>::iterator oIt = pRegions->begin();

        // reg .... ref
        if ( (*oIt)->getEnd() < (*oRefIt)->getStart())
            return 0;

        if ( iRefMax < (*oIt)->getStart() )
            return 0;

        // spool forward
        while ( (*oRefIt)->getEnd() < (*oIt)->getStart() )
        {
            ++oRefIt;
        }

        uint32_t iAccumOverlap = 0;

        while(true)
        {

            if (oRefIt == pReference->end())
                break;

            if (oIt == pRegions->end())
                break;

            uint32_t iOverlap = GenomicRegionUtils::getOverlap( *oRefIt, *oIt );

            iAccumOverlap+= iOverlap;

            /*
             *
             * three cases?
             *
             * 1: ref is longer -> ++regions
             * 2: region is longer -> ++ref
             *
             * 3:
             *
             */
            if ( (*oRefIt)->getEnd() < (*oIt)->getEnd() )
            {
                ++oRefIt;
            } else if ( (*oIt)->getEnd() < (*oRefIt)->getEnd() ) {
                ++oIt;
            } else if ( (*oIt)->getEnd() == (*oRefIt)->getEnd() ) {
                ++oRefIt;
                ++oIt;
            }


        }


        return iAccumOverlap;

    }

    static uint32_t getContained(GenomicRegion* pContained, std::vector<GenomicRegion*>* pRegions)
    {

        std::sort(pRegions->begin(), pRegions->end(), GenomicRegionUtils::sortStartAsc);

        uint32_t iContained = 0;

        std::vector<GenomicRegion*>::iterator oIt, oJt;
        for (oIt = pRegions->begin(); oIt != pRegions->end(); ++oIt)
        {

            if ( (*oIt)->contains(pContained) == true)
            {
                ++iContained;
            }

        }

        return iContained;

    }

    /**
     *
     * @param pRegions vector of GenomicRegion to iterate over to find
     * @return minimum
     */
    static uint32_t minPosition(std::vector<GenomicRegion*>* pRegions)
    {
        if (pRegions == NULL)
            return -1;

        uint32_t iMin = pRegions->at(0)->getStart();

        for (size_t i = 0; i < pRegions->size(); ++i)
            if (iMin > pRegions->at(i)->getStart())
                iMin = pRegions->at(i)->getStart();

        return iMin;

    }
    /**
     *
     * @param pRegions vector of GenomicRegion to iterate over to find
     * @return maximum
     */
    static uint32_t maxPosition(std::vector<GenomicRegion*>* pRegions)
    {
        if (pRegions == NULL)
            return -1;

        uint32_t iMax = pRegions->at(0)->getEnd();

        for (size_t i = 0; i < pRegions->size(); ++i)
            if (iMax < pRegions->at(i)->getEnd())
                iMax = pRegions->at(i)->getEnd();

        return iMax;

    }

    /**
     *
     * \brief checks whether a region overlaps with another region
     * \param oOther region to check overlap for
     * \return true if this includes other region [start, o.start, o.end, end] or [o.start,start,end,o.end] or [start,o.start,end,o.end] or [o.start,start,o.end,end]
     */
    static uint32_t getOverlap(GenomicRegion* pRegion1, GenomicRegion* pRegion2)
    {

        // this includes other
        if ((pRegion1->getStart() <= pRegion2->getStart()) && (pRegion2->getEnd() <= pRegion1->getEnd()))
            return pRegion2->getLength();

        // this is included by other
        if ((pRegion2->getStart() <= pRegion1->getStart()) && (pRegion1->getEnd() <= pRegion2->getEnd()))
            return pRegion1->getLength();

        if ( (pRegion1->getStart() <= pRegion2->getStart()) && (pRegion1->getEnd() <= pRegion2->getEnd()) && (pRegion1->getEnd() >= pRegion2->getStart()) )
            return (pRegion1->getEnd()-pRegion2->getStart()+1);

        if ( (pRegion2->getStart() <= pRegion1->getStart()) && (pRegion2->getEnd() <= pRegion1->getEnd()) && (pRegion1->getStart() <= pRegion2->getEnd()) )
            return (pRegion2->getEnd()-pRegion1->getStart()+1);

        return 0;
    }

    static void print(std::vector<GenomicRegion*>* pRegions)
    {

        std::sort(pRegions->begin(), pRegions->end(), GenomicRegionUtils::sortStartAsc);

        std::vector<GenomicRegion*>::iterator oIt, oJt;
        for (oIt = pRegions->begin(); oIt != pRegions->end(); ++oIt)
        {

            std::cout << (*oIt)->toString() << std::endl;

        }

    }
    static std::vector<GenomicRegion*>* makeUnique(std::vector<GenomicRegion*>* pRegions)
    {

        std::sort(pRegions->begin(), pRegions->end(), GenomicRegionUtils::sortStartAsc);

        std::vector<GenomicRegion*>::iterator oIt, oJt;
        for (oIt = pRegions->begin(); oIt != pRegions->end(); ++oIt)
        {

            oJt = oIt + 1;

            if (oJt == pRegions->end())
                break;

            GenomicRegion* pRegion0 = *oIt;
            GenomicRegion* pRegion1 = *oJt;

            while ( pRegion0->equals(pRegion1) == true)
            {
                delete pRegion1;

                pRegions->erase(oJt);

                oJt = oIt + 1;

                if (oJt == pRegions->end())
                    break;

                pRegion0 = *oIt;
                pRegion1 = *oJt;

            }

        }

        return pRegions;
    }

    static std::vector<GenomicRegion*>* select(std::vector<GenomicRegion*>* pRegions, GenomicRegion* pRegion, uint8_t iStrandSelection = -1)
    {
        return GenomicRegionUtils::select(pRegions, pRegion->getStart(), pRegion->getEnd(), iStrandSelection);
    }

    static std::vector<GenomicRegion*>* select(std::vector<GenomicRegion*>* pRegions, uint32_t iStart, uint32_t iEnd, uint8_t iStrandSelection = -1)
    {

        if (pRegions == NULL)
            return NULL;

        std::vector<GenomicRegion*>* pSelected = new std::vector<GenomicRegion*>();

//        uint32_t iLength = iEnd - iStart+1; // unused

        const uint32_t iElemsLength = pRegions->size();

        for (uint32_t i = 0; i < iElemsLength; ++i) {

            GenomicRegion* pElem = (*pRegions)[i];

            if (iStrandSelection < (uint8_t)-1)
            {
                if (pElem->getStrandInfo() != iStrandSelection)
                    continue;
            }

            // under no circumstances can this region be within start/end!
            if (pElem->getEnd() < iStart)
                continue;

            if (pElem->getStart() > iEnd)
                continue;

            pSelected->push_back(pElem);
        }


        return pSelected;
    }

    static std::vector<uint32_t> getCoverage(std::vector<GenomicRegion*>* pElements, GenomicRegion* pRegion) {
        return GenomicRegionUtils::getCoverage(pElements, pRegion->getStart(), pRegion->getEnd());
    }

    /**
     * \brief calculates coverage from GenomicRegion elements between iStart and iEnd
     *
     * \return vector of int of size iEnd-istart+1 which vec[i] being amount of reads at genomic position iStart+i
     */
    static std::vector<uint32_t> getCoverage(std::vector<GenomicRegion*>* pElements, uint32_t iStart, uint32_t iEnd) {

        std::vector<uint32_t> oCoverage;
        if (pElements == NULL)
            return oCoverage;

        uint32_t iLength = iEnd - iStart+1;
        oCoverage.reserve( iEnd - iStart );
        oCoverage.resize(iLength, 0);

        const uint32_t iElemsLength = pElements->size();

        for (uint32_t i = 0; i < iElemsLength; ++i) {

            GenomicRegion* pElem = (*pElements)[i];


            // under no circumstances can this region be within start/end!
            if (pElem->getEnd() < iStart)
                continue;

            if (pElem->getStart() > iEnd)
                continue;

            const uint32_t iRegionLength = pElem->getLength();
            const uint32_t iRegionStart = pElem->getStart();
            for (uint32_t j = 0; j < iRegionLength; ++j) {
                uint32_t iPosition = iRegionStart + j;

                if (iPosition >= iStart) {
                    uint32_t iIndex = iPosition - iStart;
                    if (iIndex < iLength)
                        oCoverage[iIndex] += 1;
                }
            }
        }


        return oCoverage;
    }

    static GenomicRegion* getSpanningRegion(std::vector<GenomicRegion*>* pElements)
    {
        return GenomicRegionUtils::getBoundaries(pElements);
    }

    static GenomicRegion* getBoundaries(std::vector<GenomicRegion*>* pElements) {

        uint32_t iStart = (uint32_t) -1;
        uint32_t iEnd = 0;


        for (uint32_t i = 0; i < pElements->size(); ++i) {

            GenomicRegion* pElem = pElements->at(i);

            if (pElem->getStart() < iStart)
                iStart = pElem->getStart();

            if (pElem->getEnd() > iEnd)
                iEnd = pElem->getEnd();

        }

        GenomicRegion* pBounds = new GenomicRegion(iStart, iEnd);
        return pBounds;
    }

    static GenomicRegion* find(std::vector<GenomicRegion*>* pElements, uint32_t iPosition)
    {

        for (uint32_t i = 0; i < pElements->size(); ++i) {

            GenomicRegion* pElem = pElements->at(i);

            if (pElem->contains(iPosition))
                return pElem;

        }

        return NULL;

    }

    static uint32_t getDistance(GenomicRegion* pLeft, GenomicRegion* pRight)
    {
        if (pLeft == NULL)
            return -1;
        if (pRight == NULL)
            return -1;

        if (pLeft->overlaps(pRight))
            return 0;

        // pLeft and pRight are seperate
        // [pRight] yyy [pLeft]
        if (pLeft->getStart() > pRight->getStart())
        {
            GenomicRegion* pTmp = pLeft;
            pLeft = pRight;
            pRight = pTmp;
        }

        // [pLeft] xxxxx [pRight]
        uint32_t iDist = pRight->getStart() - pLeft->getEnd();
        return iDist;

    }

    static std::string enum_REL_POSITION(GenomicRegion::REL_POSITION ePos)
    {

        switch (ePos)
        {

            case GenomicRegion::REL_POSITION::LEFT: return "LEFT";
            case GenomicRegion::REL_POSITION::RIGHT: return "RIGHT";
            case GenomicRegion::REL_POSITION::IN5P: return "IN5P";
            case GenomicRegion::REL_POSITION::IN3P: return "IN3P";
            case GenomicRegion::REL_POSITION::OVERLAP5P: return "OVERLAP5P";
            case GenomicRegion::REL_POSITION::OVERLAP3P: return "OVERLAP3P";
            case GenomicRegion::REL_POSITION::OVERLAPLEFT: return "OVERLAPLEFT";
            case GenomicRegion::REL_POSITION::OVERLAPRIGHT: return "OVERLAPRIGHT";
            case GenomicRegion::REL_POSITION::EQUAL: return "EQUAL";
            case GenomicRegion::REL_POSITION::CONTAINS: return "CONTAINS";
            case GenomicRegion::REL_POSITION::CONTAINED: return "CONTAINED";
            case GenomicRegion::REL_POSITION::ISOLATED: return "ISOLATED";

                // ERROR case
            default: return "ERROR";

        }

        return "ERROR";

    }


    static void printCoverage(std::vector<uint32_t> *pVector);
};


#endif // PROJECT_GENOMICREGIONUTILS_H
