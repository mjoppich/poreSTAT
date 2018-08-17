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

#ifndef PROJECT_XAMREADPROCESSOR_H
#define PROJECT_XAMREADPROCESSOR_H

#include <omp.h>
#include <vector>
#include <string>
#include <htslib/sam.h>
#include <vector>
#include <iostream>
#include <malloc.h>
#include <cstring>

#include "UtilsUnary.h"
#include "GenomicRegion.h"
#include "ReadRegion.h"
#include "macros.h"
#include "Sequence.h"


#if defined(__WIN32__) || defined(_WIN32) || defined(WIN32) || defined(__WINDOWS__) || defined(__TOS_WIN__)

#include <windows.h>

  inline void delay( unsigned long ms )
    {
    Sleep( ms );
    }

#else  /* presume POSIX */

#include <unistd.h>

inline void delay( unsigned long ms )
{
    usleep( ms * 1000 );
}

#endif


/**
 *
 * \brief Small wrapper for SAM/BAM/CRAM file handles in htslib to prevent mis-use
 *
 */
class XAMFile
{
public:

    /**
     *
     * @param sFileName the XAM-file to open
     * @return XAM-file wrapper for htslib
     */
    XAMFile(std::string& sFileName)
            : m_sFileName(sFileName.c_str(), sFileName.size())
    {
    }
    /**
     *
     * @param sFileName the XAM-file to open
     * @return XAM-file wrapper for htslib
     */
    XAMFile(const std::string* pFileName)
            : m_sFileName(pFileName->c_str(), pFileName->size())
    {

    }

    XAMFile( XAMFile* pOther )
        : m_sFileName(pOther->m_sFileName), m_sIdxFileName(pOther->m_sIdxFileName)
    {

    }

    /**
     *
     * @return true if file is currently open, false otherwise
     */
    bool is_open()
    {
        return (m_pXAMFile != NULL);
    }


    /**
     *
     * @return true if file has been opened successfully, false otherwise
     */
    bool open()
    {
        if (this->is_open())
        {
            this->close();
        }

        m_pXAMFile = sam_open(m_sFileName.c_str(), "r");

        if (m_pXAMFile == NULL)
            return false;

        m_pHeader = sam_hdr_read(m_pXAMFile);

        if (m_sIdxFileName.size() > 0)
        {

            this->addIndex(m_sIdxFileName);

        }

        return true;

    }

    /**
     *
     * @return true if file was closed, false otherwise (e.g. file was open)
     */
    bool close()
    {
        if (m_pXAMFile != NULL)
        {
            sam_close( m_pXAMFile );
            m_pXAMFile = NULL;

            bam_hdr_destroy(m_pHeader);
            m_pHeader = NULL;

            return true;

        }

        return false;
    }

    /**
     *
     * @return bam header of opened file
     */
    bam_hdr_t* getHeader()
    {

        if (!this->is_open())
            this->open();

        return m_pHeader;
    }


    /**
     *
     * @return XAM file file handle. Do not close this file handle or bad things will happen
     */
    samFile* getFileHandle()
    {
        if (!this->is_open())
            this->open();

        return m_pXAMFile;
    }

    void addIndex(std::string& sIndexFile)
    {
        m_sIdxFileName = std::string(sIndexFile);
        m_pXAMIdx = bam_index_load(m_sIdxFileName.c_str());
    }


    hts_itr_t* getIterator( size_t iSeqID, size_t iStart, size_t iEnd)
    {
        if (!this->is_open())
            this->open();

        if (m_pXAMIdx == NULL)
            return NULL;

        hts_itr_t* pXAMIterator = bam_itr_queryi( m_pXAMIdx ,iSeqID ,iStart ,iEnd);

        return pXAMIterator;
    }


protected:

    const std::string m_sFileName;
    std::string m_sIdxFileName = "";

    samFile* m_pXAMFile = NULL;
    bam_hdr_t* m_pHeader = NULL;

    hts_idx_t* m_pXAMIdx = NULL;

};


class XAMReadProcessor {
public:

    XAMReadProcessor(const std::string& sBAMFile, const std::string& sBAMidxFile)
            : XAMReadProcessor(&sBAMFile, &sBAMidxFile)
    {

    }

    XAMReadProcessor(const std::string* pBAMFile, const std::string* pBAMidxFile)
            : m_pBAMFile( pBAMFile == NULL ? NULL : new std::string(*pBAMFile)),
              m_pBAMidxFile( pBAMidxFile == NULL ? NULL : new std::string(*pBAMidxFile) )
    {

        m_pXAMFile = new XAMFile(m_pBAMFile);

        if (m_pBAMidxFile != NULL)
        {
            std::string sBAMIdxFile = *m_pBAMidxFile;

            m_pXAMFile->addIndex( sBAMIdxFile );
            m_pXAMIdx = bam_index_load(m_pBAMidxFile->c_str());
        } else {
            m_pXAMIdx = NULL;
        }

        pthread_mutex_init(&m_oLock, NULL);


    }

    ~XAMReadProcessor()
    {

        SAFEDEL(m_pXAMFile)
        hts_idx_destroy(m_pXAMIdx);

    }

    virtual void startAll()
    {

        std::vector<size_t> vSeqIDs = this->getSeqIDs();

        for (size_t i = 0; i < vSeqIDs.size(); ++i)
        {
            size_t iSeqID = vSeqIDs[i];

            std::cerr<< iSeqID << std::endl;

            this->start(iSeqID, 0, this->getSeqLength(iSeqID), NULL);

        }

    }

    uint32_t getSeqID(std::string* pDescr)
    {

        for (int32_t i = 0; i < m_pXAMFile->getHeader()->n_targets; ++i)
        {

            if (pDescr->compare( m_pXAMFile->getHeader()->target_name[i] ) == 0)
                return i;

        }

        return -1;

    }

    std::vector<size_t> getSeqIDs()
    {

        std::vector<size_t> vReturn;

        if (!m_pXAMFile->is_open())
            m_pXAMFile->open();

        for (int32_t i = 0; i < m_pXAMFile->getHeader()->n_targets; ++i)
        {
            vReturn.push_back(i);
        }

        return vReturn;

    }

    std::vector<std::string>* getSeqNames()
    {

        std::vector<std::string>* pSeqNames = new std::vector<std::string>();

        for (int32_t i = 0; i < m_pXAMFile->getHeader()->n_targets; ++i)
        {

            std::string sSeqName(m_pXAMFile->getHeader()->target_name[i]);
            pSeqNames->push_back(sSeqName);

        }

        return pSeqNames;

    }

    /**
     *
     * @param pRead the bam1_t read pointer
     * @param iFlag the sam flag (defined in sam.h, e.g. BAM_FSUPPLEMENTARY, BAM_FPROPER_PAIR ... ). Combine flags to query multiple flags, BAM_FSUPPLEMENTARY|BAM_FPROPER_PAIR
     * @return true if flag is set in read.
     */
    bool hasFlag(bam1_t* pRead, uint32_t iFlag)
    {
        uint32_t iReturn = pRead->core.flag & iFlag;

        if (iReturn == iFlag)
            return true;

        return false;
    }

    /**
     *
     * @param pRead bam1_t pointer
     * @return phread quality sequence (PHRED33)
     */
    std::string getQualitySequence(bam1_t* pRead)
    {
        uint8_t* pQuality = bam_get_qual(pRead);

        std::string sReturn = "";
        sReturn.resize( pRead->core.l_qseq );

        for (uint32_t i = 0; i < sReturn.size(); ++i)
        {
            sReturn[i] = pQuality[i] + 33;
        }

        return sReturn;

    }
    /**
     *
     * @param iID id of header sequence (0, .... n_targets)
     * @return lengths of this sequence, e.g. chromosome
     */
    uint32_t getSeqLength(uint32_t iID)
    {

        if ((int32_t) iID > m_pXAMFile->getHeader()->n_targets)
            return -1;

        return m_pXAMFile->getHeader()->target_len[iID];

    }

    /**
     *
     * @param iID id of header sequence (0, .... n_targets)
     * @return sequence name as in sam header (e.g. chr1)
     */
    std::string getSeqName(uint32_t iID)
    {
        return std::string( m_pXAMFile->getHeader()->target_name[iID] );

    }

    virtual void startAllParallel(int iChunkSize, void* pData)
    {
        uint32_t iReads = 0;
        std::vector<bam1_t*>* pReads = new std::vector<bam1_t*>();

        htsFile* pXAMFile = m_pXAMFile->getFileHandle();
        bam_hdr_t* pHeader = m_pXAMFile->getHeader();

        std::cout << "Preparing SAM FILE" << std::endl;

        bam1_t* pRes = bam_init1();

        


        int iActiveTasks = 0;

#pragma omp parallel
        {

            int iThreads = omp_get_max_threads();

            std::cout << "Starting up OPENMP with " << iThreads << " threads. " << omp_get_thread_num() << std::endl;

#pragma omp single nowait
            {


                while ( sam_read1(pXAMFile, pHeader, pRes)  >= 0)
                {
                
                    bam1_t* pRead = bam_init1();
                    bam_copy1(pRead, pRes);
                    pReads->push_back( pRead );

                    if (pReads-> size() > iChunkSize)
                    {

                        {

                            // wait until a task has finished
                            while(iActiveTasks > iThreads + 1)
                            {
                                delay(1000);
                            }
                        }


#pragma omp atomic
                        ++iActiveTasks;

#pragma omp task shared(iActiveTasks)
                        {
                            #pragma omp critical
                            std::cout << "Currently tasking in task " << omp_get_thread_num() << std::endl;

                            for (size_t i = 0; i < pReads->size(); ++i)
                            {
                                this->process(pReads->at(i), 0, pData);
                            }

                            this->deleteReads(pReads);

#pragma omp atomic
                            --iActiveTasks;

                        }

                        pReads = new std::vector<bam1_t*>();


                    }


                    ++iReads;
                }

#pragma omp taskwait
            }
        }



        bam_destroy1(pRes);
        sam_close(pXAMFile);
        //pthread_mutex_unlock(&m_oLock);

    }



    virtual std::vector<bam1_t*>* startParallel(uint32_t iChunkSize, uint32_t iSeqID, uint32_t iStart, uint32_t iEnd, std::function<void ()> oAfterTask = [] () {} ,void* pData=NULL)
    {

        std::cout << "Start Parallel launched" << std::endl;


        //pthread_mutex_lock(&m_oLock);

        samFile* pXAMFile;
        hts_idx_t* pXAMIdx;

        std::vector<bam1_t*>* pReturn = new std::vector<bam1_t*>();

//        kstring_t ks = { 0, 0, NULL };

        hts_itr_t* pXAMIterator = NULL;
        bam1_t* pRes;

        pXAMFile = sam_open(m_pBAMFile->c_str(), "r");
        pXAMIdx = bam_index_load(m_pBAMidxFile->c_str());

        std::cout << "Preparing BAM FILE" << std::endl;

        pXAMIterator = sam_itr_queryi(pXAMIdx, iSeqID, iStart, iEnd);
        pRes = bam_init1();



        uint32_t iReads = 0;

        std::vector<bam1_t*>* pReads = new std::vector<bam1_t*>();

        int iActiveTasks = 0;

        {

            int iThreads = omp_get_max_threads();

            std::cout << "Starting up OPENMP with " << iThreads << " threads. " << omp_get_thread_num() << std::endl;

#pragma omp single nowait
            {


                while (bam_itr_next(pXAMFile, pXAMIterator, pRes) > 0)
                {

                    bam1_t* pRead = bam_init1();
                    bam_copy1(pRead, pRes);
                    pReads->push_back( pRead );

                    if (pReads-> size() > iChunkSize)
                    {

                        {

                            // wait until a task has finished
                            while(iActiveTasks > iThreads + 1)
                            {
                                delay(1000);
                            }
                        }


#pragma omp atomic
                        ++iActiveTasks;

                        std::cout << "started chunk " << iActiveTasks << std::endl;


#pragma omp task shared(iActiveTasks)
                        {

                            std::cout << "in chunk " << iActiveTasks << std::endl;

                            for (size_t i = 0; i < pReads->size(); ++i)
                            {
                                this->process(pReads->at(i), iSeqID, pData);
                            }

                            this->deleteReads(pReads);

                            oAfterTask();

#pragma omp atomic
                            --iActiveTasks;

                            std::cout << "finished chunk " << iActiveTasks << std::endl;



                        }

                        pReads = new std::vector<bam1_t*>();


                    }


                    ++iReads;
                }

#pragma omp taskwait
            }
        }



        bam_destroy1(pRes);
        hts_itr_destroy(pXAMIterator);

        sam_close(pXAMFile);
        hts_idx_destroy(pXAMIdx);
        //pthread_mutex_unlock(&m_oLock);


        return pReturn;
    }


    virtual std::vector<bam1_t*>* start(uint32_t iSeqID, uint32_t iStart, uint32_t iEnd, void* pData=NULL)
    {
        //pthread_mutex_lock(&m_oLock);

        samFile* pXAMFile;
        hts_idx_t* pXAMIdx;

        std::cout << "Before opening file" << std::endl;

        pXAMFile = sam_open(m_pBAMFile->c_str(), "r");
        pXAMIdx = bam_index_load(m_pBAMidxFile->c_str());

        std::vector<bam1_t*>* pReturn = new std::vector<bam1_t*>();

//        kstring_t ks = { 0, 0, NULL };

        hts_itr_t* pXAMIterator = bam_itr_queryi( pXAMIdx ,iSeqID ,iStart ,iEnd);
        bam1_t* pRes = bam_init1();

        uint32_t iReads = 0;

        std::cout << "Starting iterations" << std::endl;


        while (bam_itr_next(pXAMFile, pXAMIterator, pRes) > 0)
        {

            this->process(pRes, iSeqID, pData);

            ++iReads;
        }

        bam_destroy1(pRes);
        hts_itr_destroy(pXAMIterator);

        sam_close(pXAMFile);
        hts_idx_destroy(pXAMIdx);
        //pthread_mutex_unlock(&m_oLock);


        return pReturn;
    }

    GenomicRegion* getSpannedRegion(bam1_t* pRead)
    {
        uint32_t iLeftStart = pRead->core.pos;
        uint32_t iCIGARitems = pRead->core.n_cigar;
        uint32_t* pCIGARs = bam_get_cigar(pRead);

//        uint32_t iMatch = 0;
        uint32_t iPosition = getReadStart(pRead);

        for (size_t i = 0; i < iCIGARitems; ++i)
        {
            uint32_t iCIGARop = bam_cigar_op(pCIGARs[i]);
//            char cCIGARop = BAM_CIGAR_STR[ iCIGARop ];
            uint32_t iCIGARlen = bam_cigar_oplen(pCIGARs[i]);

            if (iCIGARop == 1)
                continue;

            if (iCIGARop == 0) // MATCH M
            {
                iPosition += iCIGARlen;
            }

            if (iCIGARop == 1) // INSERTION TO REFERENCE I
                iPosition = iPosition + 0;

            if (iCIGARop == 2) // DELETION FROM REFERENCE D
                iPosition = iPosition + iCIGARlen;

            if (iCIGARop == 3) // SKIPPED REGION IN REFERENCE N
                iPosition = iPosition + iCIGARlen;

            if (iCIGARop == 4) // SOFT CLIPPING S
                iPosition = iPosition + 0;

            if (iCIGARop == 5) // HARD CLIPPING H
                iPosition = iPosition + iCIGARlen;

            if (iCIGARop == 6) // PADDING P
                iPosition = iPosition + iCIGARlen;

            if (iCIGARop == 7) // SEQUENCE MATCH =
                iPosition = iPosition + iCIGARlen;

            if (iCIGARop == 8) // SEQUENCE MISMATCH X
                iPosition = iPosition + iCIGARlen;

        }

        return new GenomicRegion(iLeftStart, iPosition-1);
    }

    std::vector<GenomicRegion*> getMatchedRegions(bam1_t* pRead)
    {
        return this->getRegions(pRead, 'M');
    }

    std::vector<GenomicRegion*> getRegions(bam1_t* pRead, char cRegion)
    {
        std::vector<char> vCIGARs;
        vCIGARs.push_back(cRegion);

        return this->getRegions(pRead, &vCIGARs);
    }

    std::vector<GenomicRegion*> getAllRegions(bam1_t* pRead)
    {

        std::vector<char> vCIGARs;
        vCIGARs.push_back('M');
        vCIGARs.push_back('N');
        vCIGARs.push_back('X');

        return this->getRegions(pRead, &vCIGARs);

    }

    ReadRegion* handleCIGARop(std::vector<GenomicRegion*>* pRegions, ReadRegion* pCurrentRegion, std::vector<char>* pOperations, char cOperation, uint32_t iPosition, uint32_t iCIGARlen)
    {
        bool bContainedCIGAR = (pOperations != NULL ) ? UtilsUnary<char>::contains(pOperations, cOperation) : true;

        if ( bContainedCIGAR )
        {
            if (pCurrentRegion == NULL) {
                pCurrentRegion = new ReadRegion(iPosition, iPosition + iCIGARlen - 1, cOperation, iCIGARlen);
                //pCurrentRegion->setAttribute("CIGAR", cOperation);
            } else {

                if ( pCurrentRegion->getCIGAR() != cOperation )
                {
                    // end current region
                    pCurrentRegion->setEnd(iPosition-1);
                    pRegions->push_back( pCurrentRegion );

                    // start next region
                    pCurrentRegion = new ReadRegion(iPosition, iPosition + iCIGARlen - 1, cOperation, iCIGARlen);
                    //pCurrentRegion->setAttribute("CIGAR", cOperation);
                } else {

                    // still same CIGAR (for whatever reason)
                    pCurrentRegion->setCIGARLenght( pCurrentRegion->getCIGARLength() + iCIGARlen );

                }

            }

        }


        if (!bContainedCIGAR)
        {
            if (pCurrentRegion != NULL)
            {

                pCurrentRegion->setEnd(iPosition-1);
                pRegions->push_back( (GenomicRegion*) pCurrentRegion );

                pCurrentRegion = NULL;
            }

        }

        return pCurrentRegion;
    }

    inline uint32_t getReadStart(bam1_t* pRead)
    {
        return pRead->core.pos + 1;
    }

    std::string getReadName(bam1_t* pRead)
    {

        char* pReadName = bam_get_qname(pRead);

        return std::string(pReadName, pRead->core.l_qname-1);
    }

    std::string getReadSequence(bam1_t* pRead)
    {
        return getReadSequence( bam_get_seq(pRead), 0, pRead->core.l_qseq );
    }

    std::string getReadSequence(uint8_t* pSequence, uint32_t iPositionStart, uint32_t iLength)
    {
        std::string sReturn = "";

        for (size_t i = iPositionStart; i < iPositionStart + iLength; ++i)
        {
            uint8_t iSeqAtI = bam_seqi(pSequence, i);

            // 1 for A, 2 for C, 4 for G, 8 for T and 15 for N
            char cBaseAtI = 'N';
            switch (iSeqAtI)
            {

                case 1: cBaseAtI =  'A'; break;
                case 2: cBaseAtI =  'C'; break;
                case 4: cBaseAtI =  'G'; break;
                case 8: cBaseAtI =  'T'; break;
                case 15: cBaseAtI = 'N'; break;

            }
            sReturn.append(1, cBaseAtI);
        }


        return sReturn;


    }

    std::vector<GenomicRegion*> getRegions(bam1_t* pRead, std::vector<char>* pOperations, bool bGetReadSequence = false)
    {
        uint32_t iCIGARitems = pRead->core.n_cigar;
        uint32_t* pCIGARs = bam_get_cigar(pRead);

        std::vector<GenomicRegion*> vRegions;// = new std::vector<GenomicRegion>();

//        uint32_t iMatch = 0;
        uint32_t iRefPosition = getReadStart(pRead);
        uint32_t iQueryPosition = 0;
        int32_t iReadLength = pRead->core.l_qseq;

        ReadRegion* pCurrentRegion = NULL;

        bool bRev = bam_is_rev(pRead);

        uint8_t* pSequence = bam_get_seq(pRead);

        for (size_t i = 0; i < iCIGARitems; ++i)
        {
            uint32_t iCIGARop = bam_cigar_op(pCIGARs[i]);
//            char cCIGARop = BAM_CIGAR_STR[ iCIGARop ];
            uint32_t iCIGARlen = bam_cigar_oplen(pCIGARs[i]);


            /*
            if ((iCIGARop == 0) && (iCIGARlen > 50))
            {

                size_t iTotal = 0;
                for (size_t j = 0; j < i; ++j)
                {

                    uint32_t qiCIGARop = bam_cigar_op(pCIGARs[j]);
                    char qcCIGARop = BAM_CIGAR_STR[ qiCIGARop ];
                    uint32_t qiCIGARlen = bam_cigar_oplen(pCIGARs[j]);

                    std::cout << qiCIGARop << " " << qcCIGARop << " " << qiCIGARlen << std::endl;
                    iTotal += qiCIGARlen;

                }


                std::cout << "asd" << iTotal << std::endl;
            }

            */

            if (iCIGARop == 0) // MATCH M
            {

                pCurrentRegion = this->handleCIGARop(&vRegions, pCurrentRegion, pOperations, 'M', iRefPosition, iCIGARlen);
                iRefPosition += iCIGARlen;
                iQueryPosition += iCIGARlen;
            }

            if (iCIGARop == 1) // INSERTION TO REFERENCE I
            {

                pCurrentRegion = this->handleCIGARop(&vRegions, pCurrentRegion, pOperations, 'I', iRefPosition, iCIGARlen);
                iRefPosition = iRefPosition + 0;
                iQueryPosition += iCIGARlen;
            }

            if (iCIGARop == 2) { // DELETION FROM REFERENCE D

                pCurrentRegion = this->handleCIGARop(&vRegions, pCurrentRegion, pOperations, 'D', iRefPosition, iCIGARlen);
                iRefPosition += iCIGARlen;
            }

            if (iCIGARop == 3) {// SKIPPED REGION IN REFERENCE N

                pCurrentRegion = this->handleCIGARop(&vRegions, pCurrentRegion, pOperations, 'N', iRefPosition, iCIGARlen);
                iRefPosition += iCIGARlen;
            }

            if (iCIGARop == 4)
            {
                // SOFT CLIPPING S
                pCurrentRegion = this->handleCIGARop(&vRegions, pCurrentRegion, pOperations, 'S', iRefPosition, iCIGARlen);
                //iRefPosition += iCIGARlen;
                iQueryPosition += iCIGARlen;
            }


            if (iCIGARop == 5) {
                // HARD CLIPPING H
                pCurrentRegion = this->handleCIGARop(&vRegions, pCurrentRegion, pOperations, 'H', iRefPosition, iCIGARlen);
                //iRefPosition += iCIGARlen;
            }

            if (iCIGARop == 6) // PADDING P
                iRefPosition += iCIGARlen;

            if (iCIGARop == 7) // SEQUENCE MATCH =
            {
                iRefPosition += iCIGARlen;
                iQueryPosition += iCIGARlen;
            }


            if (iCIGARop == 8) // SEQUENCE MISMATCH X
            {
                pCurrentRegion = this->handleCIGARop(&vRegions, pCurrentRegion, pOperations, 'X', iRefPosition, iCIGARlen);
                iRefPosition += iCIGARlen;
                iQueryPosition += iCIGARlen;
            }


            if (bGetReadSequence)
            {

                if ( (iCIGARop == 1) || (iCIGARop == 6) )
                {
                    continue;
                }

                // now retrieve sequence

                std::string sSequence = "";

                if (!bRev)
                {
                    sSequence = this->getReadSequence(pSequence, iQueryPosition - iCIGARlen, iCIGARlen);

                } else {

//                    uint32_t iStartPos = iQueryPosition - iCIGARlen;

                    sSequence = this->getReadSequence(pSequence, iQueryPosition - iCIGARlen, iCIGARlen);

                    Sequence oSeq(sSequence);
                    //oSeq.reverseComplement();

                    sSequence = oSeq;

                }

                pCurrentRegion->setSequence(sSequence);

                std::string sQueryName = this->getReadID(pRead);
                pCurrentRegion->setQueryName(sQueryName);

            }


        }

        if (pCurrentRegion != NULL)
            vRegions.push_back( (GenomicRegion*) pCurrentRegion);

        for (size_t i = 0; i < vRegions.size(); ++i)
        {
            ReadRegion* pReadRegion = (ReadRegion*)vRegions.at(i);

            if (m_bTurnReadDirection)
            {
                if (bRev) // not reverse
                    pReadRegion->setStrand('+');
                else // reverse
                    pReadRegion->setStrand('-');

            } else {

                if (!bRev) // not reverse
                    pReadRegion->setStrand('+');
                else // reverse
                    pReadRegion->setStrand('-');

            }




        }

        if (iQueryPosition != (uint32_t) iReadLength)
        {
            std::cerr << iQueryPosition << " " << iReadLength << std::endl;
        }

        return vRegions;
    }


    std::string getReadID(bam1_t* pRead)
    {
        char* pQueryName = bam_get_qname(pRead);
        std::string sQueryName(pQueryName);

        return sQueryName;
    }

    void setTurnReadDirection(bool bTurnDirection)
    {
        m_bTurnReadDirection = bTurnDirection;
    }

    bool acceptReadMapping(bam1_t* pRead)
    {

        bool bAccept = true;

        uint8_t *aux = bam_aux_get(pRead, "NH");
        if (aux){

            int32_t refvalue;
            refvalue = bam_aux2i(aux);
            bAccept = ( refvalue == 1);
        }

        return bAccept;
    }

    void deleteReads(std::vector<bam1_t*>* pReads)
    {
        const uint32_t iReadsLength = pReads->size();
        for (uint32_t i = 0; i < iReadsLength; ++i)
            bam_destroy1((*pReads)[i]);

        delete pReads;
    }


    bool readIsMapped(bam1_t* pRead)
    {
        return ((pRead->core.flag & BAM_FUNMAP) == 0);
    }

    bool readIsRead1(bam1_t* pRead)
    {
        return ((pRead->core.flag & BAM_FREAD1) != 0);
    }

    bool readIsRead2(bam1_t* pRead)
    {
        return ((pRead->core.flag & BAM_FREAD2) != 0);
    }

    bool readIsSecondary(bam1_t* pRead)
    {
        return ((pRead->core.flag & BAM_FSECONDARY) != 0);
    }

    bool readIsReverse(bam1_t* pRead)
    {
        return ((pRead->core.flag & BAM_FREVERSE) != 0);
    }

    bool readMateReverse(bam1_t* pRead)
    {
        return ((pRead->core.flag & BAM_FMREVERSE) != 0);
    }

    bool readIsPrimary(bam1_t* pRead)
    {
        return ((pRead->core.flag & BAM_FSECONDARY) == 0) && ((pRead->core.flag & BAM_FSUPPLEMENTARY) == 0);
    }

    bool readMateMapped(bam1_t* pRead)
    {
        return ((pRead->core.flag & BAM_FMUNMAP) == 0);
    }

    bool readMapped(bam1_t* pRead)
    {
        return ((pRead->core.flag & BAM_FUNMAP) == 0);
    }

    bool readHasMate(bam1_t* pRead)
    {

        int iNextSeqID = pRead->core.mtid;
        int iNextSeqPos = pRead->core.mpos;

        if (iNextSeqID == -1)
            return false;

        // TODO can this happen?
        if (iNextSeqPos == -1)
            return false;

        return true;
    }

protected:
    virtual void process(bam1_t* pRead, uint32_t iSeqID, void* pData) = 0;

    pthread_mutex_t m_oLock;
    const std::string* m_pBAMFile;
    const std::string* m_pBAMidxFile;

    XAMFile* m_pXAMFile;

    hts_idx_t* m_pXAMIdx;
    bool m_bTurnReadDirection = true;

};


#endif //PROJECT_XAMREADPROCESSOR_H
