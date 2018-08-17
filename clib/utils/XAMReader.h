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

#ifndef PROJECT_XAMREADER_H
#define PROJECT_XAMREADER_H

#include <string>
#include <htslib/sam.h>
#include <vector>
#include <iostream>
#include <malloc.h>
#include <cstring>

#include "Maths.h"
#include "XAMReadProcessor.h"
#include "GenomicRegionUtils.h"

/**
 *
 * \brief a class which parses sam/bam files and return all reads within a specific region (if asked)
 *
 */
class XAMReader : public XAMReadProcessor {
public:



    XAMReader(std::string& sBAMFile, std::string& sBAMidxFile)
    : XAMReader(&sBAMFile, &sBAMidxFile)
    {

    }

    XAMReader(std::string* pBAMFile, std::string* pBAMidxFile)
    : XAMReadProcessor(pBAMFile, pBAMidxFile)
    {


    }

    std::vector<bam1_t*>* getReads(GenomicRegion* pRegion)
    {

        uint32_t iStart = pRegion->getStart();
        uint32_t iEnd = pRegion->getEnd();
        uint32_t iSeqID = pRegion->getChromosomeID();

        return this->getReads(iSeqID, iStart, iEnd);

    }

    std::vector<bam1_t*>* getReads(GenomicRegion* pRegion, std::string* pSeqName)
    {

        uint32_t iStart = pRegion->getStart();
        uint32_t iEnd = pRegion->getEnd();

        uint32_t iSeqID = this->getSeqID( pSeqName );

        return this->getReads(iSeqID, iStart, iEnd);

    }

    void process(bam1_t* pRead, uint32_t iSeqID, void* pData)
    {

        bam1_t* pCopy = bam_dup1(pRead);

        ((std::vector<bam1_t*>*)pData)->push_back(pCopy);

    }

    std::vector<bam1_t*>* getReads(uint32_t iSeqID, uint32_t iStart, uint32_t iEnd)
    {
        std::vector<bam1_t*>* pReturn = new std::vector<bam1_t*>();

        this->start(iSeqID, iStart, iEnd, (void*) pReturn);

        return pReturn;

    }

    size_t extractRegions(std::vector<bam1_t*>* pReads, std::vector<GenomicRegion*>* pMatchedRegions, std::vector<GenomicRegion*>* pSplicedRegions)
    {
        std::vector<char> vCIGARs;
        vCIGARs.push_back('M');
        vCIGARs.push_back('N');
        vCIGARs.push_back('X');

        return this->extractRegions(pReads, &vCIGARs, pMatchedRegions, pSplicedRegions);
    }

    size_t extractRegions(std::vector<bam1_t*>* pReads, std::vector<char>* pCIGARs, std::vector<GenomicRegion*>* pMatchedRegions, std::vector<GenomicRegion*>* pSplicedRegions)
    {

        size_t iValidReads = 0;

        // prepare things
        pMatchedRegions->reserve( pReads->size() );
        pSplicedRegions->reserve( pReads->size() );

        // actually extract regions

        const uint32_t iReadLength = pReads->size();
        for (uint32_t i = 0; i < iReadLength; ++i )
        {
            // this filters ambiguous reads

            bam1_t* pRead = (*pReads)[i];
            uint32_t iLength = pRead->core.l_qseq;

            bool bRemoveRead = ! (this->acceptReadMapping(pRead));

            if (bRemoveRead)
            {
                //bam_destroy1(pRead);
                //oIt = pReads->erase(oIt);
                continue;

            }

            std::vector<GenomicRegion*> vMatchedRegion = this->getRegions(pRead, pCIGARs);
            uint32_t iMatchedLength = GenomicRegionUtils::getLengthVec(&vMatchedRegion);
            bRemoveRead |= (float) iMatchedLength < 0.5f * (float) iLength; // must have at least 50% matched


            if (bRemoveRead)
            {
                GenomicRegionUtils::deleteVec(&vMatchedRegion);
                //bam_destroy1(pRead);
                //oIt = pReads->erase(oIt);
                continue;

            }

            ++iValidReads;

            const uint32_t iRegionsLength = vMatchedRegion.size();
            for (uint32_t j = 0; j < iRegionsLength; ++j)
            {
                ReadRegion* pRegion = (ReadRegion*)vMatchedRegion[j];

                // if region represents a split we do not want it in matchedregions
                if (pRegion->getCIGAR() == 'N')
                    pSplicedRegions->push_back((GenomicRegion*)pRegion);
                else
                    pMatchedRegions->push_back((GenomicRegion*)pRegion);

            }




        }

        return iValidReads;

    }

    float getCoverage(GenomicRegion* pElement)
    {
        std::vector<GenomicRegion*>* pMatchedRegions = new std::vector<GenomicRegion*>();

        std::vector<bam1_t*>* pReads = this->getReads(pElement);

        std::vector<bam1_t*>::iterator oIt;
        uint32_t iReadCount = 0;
        for (oIt = pReads->begin(); oIt != pReads->end(); )
        {

            ++iReadCount;

            // this filters ambiguous reads

            bam1_t* pRead = *oIt;
            uint32_t iLength = pRead->core.l_qseq;

            bool bRemoveRead = false;

            uint8_t *aux = bam_aux_get(pRead, "NH");
            if (aux){

                int32_t refvalue;
                refvalue = bam_aux2i(aux);
                bRemoveRead |= ( refvalue != 1);
            }

            std::vector<GenomicRegion*> vMatchedRegion = this->getAllRegions(pRead);
            uint32_t iMatchedLength = GenomicRegionUtils::getLengthVec(&vMatchedRegion);
            bRemoveRead |= (float) iMatchedLength < 0.5f * (float) iLength; // must have at least 50% matched


            if (bRemoveRead)
            {
                GenomicRegionUtils::deleteVecContent(&vMatchedRegion);
                bam_destroy1(pRead);

                oIt = pReads->erase(oIt);
                continue;

            }

            // is this read within read selection or not?
            if (GenomicRegionUtils::allContained((GenomicRegion*) pElement, &vMatchedRegion))
            {
                pMatchedRegions->insert(pMatchedRegions->end(), vMatchedRegion.begin(), vMatchedRegion.end());

            } else {
                bam_destroy1(pRead);
                GenomicRegionUtils::deleteVecContent(&vMatchedRegion);

                oIt = pReads->erase(oIt);
                continue;
            }

            ++oIt;

        }

        this->deleteReads(pReads);

        std::vector<uint32_t> oCoverages = GenomicRegionUtils::getCoverage(pMatchedRegions, (GenomicRegion*) pElement);
        float fCoverage = Maths<uint32_t>::getMean(&oCoverages);
        GenomicRegionUtils::deleteVec(pMatchedRegions);

        return fCoverage;

    }





};


#endif //PROJECT_XAMREADER_H
