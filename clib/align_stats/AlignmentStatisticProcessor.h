#include <vector>
#include <map>
#include "../utils/FASTAreader.h"
#include <assert.h>     /* assert */

#include "../utils/XAMReadProcessor.h"
#include "../utils/GenomicRegion.h"

#include "ReadStats.h"


class AlignmentStatisticProcessor : public XAMReadProcessor {
public:


    AlignmentStatisticProcessor(std::string* pBAMFile, std::string* pBAMidxFile, FASTAreader* pFASTAReader)
    : XAMReadProcessor(pBAMFile, pBAMidxFile), pFASTAReader(pFASTAReader)
    {

        this->vCIGARs.push_back('M');
        this->vCIGARs.push_back('S');
        this->vCIGARs.push_back('H');
        this->vCIGARs.push_back('I');
        this->vCIGARs.push_back('D');
        this->vCIGARs.push_back('N');
    
        this->pMsPerRead = new std::vector<uint32_t>();
        this->pSeqNames = this->getSeqNames();
    }

    virtual ~AlignmentStatisticProcessor()
    {
        delete this->pMsPerRead;
        delete this->pSeqNames;
    }


    virtual void startAllParallel(size_t iChunkSize, void* pData, size_t iMaxReads=-1)
    {
        uint32_t iReads = 0;
        std::vector<bam1_t*>* pReads = new std::vector<bam1_t*>();

        std::vector<PythonReadStats>* pResultStats = (std::vector<PythonReadStats>*) pData;


        htsFile* pXAMFile = m_pXAMFile->getFileHandle();
        bam_hdr_t* pHeader = m_pXAMFile->getHeader();

        std::cout << "Preparing SAM FILE" << std::endl;

        bam1_t* pRes = bam_init1();


        m_iNumThreads = omp_get_max_threads();

        int iActiveTasks = 0;

#pragma omp parallel num_threads(m_iNumThreads) shared(iActiveTasks)
        {

#pragma omp master
            {

                int iThreads = omp_get_num_threads();

                std::cout << "Starting up OPENMP with " << omp_get_num_threads() << " threads. Hello from master thread " << omp_get_thread_num() << std::endl;

                std::cerr << "Got a sequence file" << std::endl;


                while ( sam_read1(pXAMFile, pHeader, pRes)  >= 0)
                {

                    bam1_t* pRead = bam_init1();
                    bam_copy1(pRead, pRes);
                    pReads->push_back( pRead );

                    if ((iMaxReads != -1) && (iReads > iMaxReads))
                    {
                        break;
                    }

                    if (pReads->size() > iChunkSize)
                    {

                        iReads += pReads->size();
#pragma omp critical
                        {
                            std::cout << "Processed reads: " << iReads << std::endl;
                        }
                        {

                            // wait until a task has finished
                            while(iActiveTasks > iThreads + 1)
                            {
                                //std::cout << "Waiting for task to finish " << omp_get_thread_num() << " because " << iActiveTasks << std::endl;
                                delay(100);
                            }
                        }


#pragma omp critical
                        {
                            ++iActiveTasks;
                        }


#pragma omp task shared(iActiveTasks)
                        {
//#pragma omp critical
                            std::cout << "Currently tasking in task " << omp_get_thread_num() << std::endl;
                            std::vector<PythonReadStats> oLocalRes;
                            oLocalRes.reserve(pReads->size());

                            for (size_t i = 0; i < pReads->size(); ++i)
                            {
                                this->process(pReads->at(i), 0, &oLocalRes);
                            }

                            this->deleteReads(pReads);

#pragma omp critical
                            {

                                pResultStats->insert(pResultStats->end(), oLocalRes.begin(), oLocalRes.end());
                                --iActiveTasks;
                            }


                        }

                        pReads = new std::vector<bam1_t*>();


                    }

                }

                if (pReads->size() > 0)
                {

#pragma omp task shared(iActiveTasks)
                        {
#pragma omp critical
                            std::cout << "Currently tasking in task " << omp_get_thread_num() << std::endl;
                            std::vector<PythonReadStats> oLocalRes;
                            oLocalRes.reserve(pReads->size());

                            for (size_t i = 0; i < pReads->size(); ++i)
                            {
                                this->process(pReads->at(i), 0, &oLocalRes);
                            }

                            this->deleteReads(pReads);

#pragma omp critical
                            {

                                pResultStats->insert(pResultStats->end(), oLocalRes.begin(), oLocalRes.end());

                                --iActiveTasks;
                            }


                        }

                }

#pragma omp taskwait
            }
        }



        bam_destroy1(pRes);
        sam_close(pXAMFile);
        //pthread_mutex_unlock(&m_oLock);

    }



    int seenReads = 0;
    size_t seenRegions = 0;
    std::vector<char> vCIGARs;
    std::vector<uint32_t>* pMsPerRead;
    std::vector<std::string>* pSeqNames;

    FASTAreader* pFASTAReader;
    
protected:



    uint32_t getAlignQuality(bam1_t* pRead)
    {
        return pRead->core.qual;
    }

    std::vector<float> getSeqQuality(bam1_t* pRead)
    {
        std::vector<uint8_t> allQuals = this->getIQualitySequence(pRead);
        std::sort(allQuals.begin(), allQuals.end());

        std::vector<float> vRet;

        vRet.push_back((float) allQuals[0]);

        float fMedian;

        if (allQuals.size() % 2 == 0)
        {
            size_t iMedPos = allQuals.size() / 2;
            fMedian = allQuals[iMedPos] + allQuals[iMedPos-1];
            fMedian /= 2.0;

        } else {
            size_t iMedPos = (size_t) std::ceil(allQuals.size() / 2.0);
            fMedian = allQuals[iMedPos];
        }

        vRet.push_back(fMedian);
        vRet.push_back((float) allQuals[allQuals.size()-1]);

        return vRet;

    }

    uint32_t countGCBases(std::string& sseq)
    {
        uint32_t iGC =0;

        for (size_t i = 0; i < sseq.size(); ++i)
        {
            if ((sseq[i] == 'G') || (sseq[i]=='C'))
            {
                iGC += 1;
            }

        }

        return iGC;
    }


    virtual void process(bam1_t* pRead, uint32_t iSeqID, void* pData)
    {

        AlignedReadStats* pReadStats = new AlignedReadStats();
        std::string sReadName = this->getReadName(pRead);
        pReadStats->sReadID = sReadName;
        //std::cout << "sreadname " << sReadName << std::endl;
        pReadStats->aligned = false;

        if (this->readIsMapped(pRead) && !(this->readIsSecondary(pRead)))
        {
            this->seenReads += 1;
            pReadStats->aligned = true;

            std::string sAlignedSeqName = this->getSeqName(pRead->core.tid);

            //std::cout << sReadName << " " << sAlignedSeqName << std::endl;

            pReadStats->iAlignQual = this->getAlignQuality(pRead);

            std::vector<float> vSeqQuals = this->getSeqQuality(pRead);
            pReadStats->fSeqQual_min = vSeqQuals[0];
            pReadStats->fSeqQual_median = vSeqQuals[1];
            pReadStats->fSeqQual_max = vSeqQuals[2];

            //std::cout << sReadName << " " << sAlignedSeqName << " Qualitites done" << std::endl;

            std::vector<GenomicRegion*> allAlignedRegions = this->getRegions(pRead, &(this->vCIGARs));

            ReadRegion* pStartRegion = (ReadRegion*)allAlignedRegions[0];
            ReadRegion* pEndRegion = (ReadRegion*)allAlignedRegions[allAlignedRegions.size()-1];

            uint32_t iRefStart = -1;
            uint32_t iRefEnd = -1;

            iRefStart = pStartRegion->getStart();
            iRefEnd = pEndRegion->getEnd();

            uint32_t iLongestMatched = 0;
            uint32_t iReadGCBases = 0;
            uint32_t iRefGCBases = 0;
            uint32_t iMatchedRegionsBases=0;
            uint32_t iReadMatched = 0;

            uint32_t iReadMatches = 0;
            uint32_t iReadMismatches = 0;
            uint32_t iReadM = 0;

            /*
             * TODO this might break if read does not start with M
             */
            uint32_t iReadStart = 0;
            uint32_t iReadEnd = -1;
            iReadStart = pStartRegion->pReadRegion->getStart();
            iReadEnd = pEndRegion->pReadRegion->getEnd();

            pReadStats->iReadLength = iReadEnd-iReadStart;
            pReadStats->iRefLength = iRefEnd-iRefStart;

            uint32_t iRetRes = allAlignedRegions.size();
            pReadStats->seqCoverages.reserve(allAlignedRegions.size());

            if (iRetRes)
            {

                //std::cout << sReadName << " " << sAlignedSeqName << " iRetRes" << std::endl;

                std::string sReadSeq = this->getReadSequence(pRead);
                //std::cout << sReadName << " " << sAlignedSeqName << std::endl;

                for (uint32_t i = 0; i < iRetRes; ++i)
                {
                   
                    ReadRegion* pMRegion = (ReadRegion*)allAlignedRegions.at(i);

                    pReadStats->cigar2len.add( pMRegion->getCIGAR(), pMRegion->getCIGARLength() );

                    if (pMRegion->getCIGAR() == 'M')
                    {

                        iReadM += pMRegion->getCIGARLength();

                        std::string sRefRegion;
                        std::string sReadRegion;

                        sRefRegion = this->pFASTAReader->retrieveSequence(pMRegion, sAlignedSeqName);
                        sReadRegion = sReadSeq.substr(pMRegion->pReadRegion->getStart(), pMRegion->pReadRegion->getLength()-1);

                        /*
                        std::cerr << "Region " << pMRegion->toString() << std::endl;
                        std::cerr << "ChrName" << sAlignedSeqName << std::endl;
                        std::cerr << "Ref    " << sRefRegion << std::endl;
                        std::cerr << "Read   " << sReadRegion << std::endl;
                        */
                        

                        char* pSeqC = (char*) malloc(sizeof(char) * (sAlignedSeqName.size()+1));
                        sAlignedSeqName.copy(pSeqC, sAlignedSeqName.size());
                        pSeqC[sAlignedSeqName.size()] = '\0';

                        //std::cout << "coverage seq " << pSeqC << " " << sAlignedSeqName << std::endl;
                        // add coverage
                        pReadStats->seqCoverages.push_back(SeqCoverage(pRead->core.tid, pMRegion->getStart(), pMRegion->getEnd()));

                        iLongestMatched = std::max(pMRegion->pReadRegion->getLength(), iLongestMatched);

                        iRefGCBases += this->countGCBases(sRefRegion);
                        iReadGCBases += this->countGCBases(sReadRegion);
                        iMatchedRegionsBases += sRefRegion.size();

                        uint32_t iMatches = 0;
                        uint32_t iMisMatches = 0;

                        uint32_t equals = 0;
                        uint32_t mismatched = 0;

                        for (uint32_t j = 0; j < sReadRegion.size(); ++j)
                        {
                            char cRead = sReadRegion[j];
                            // because some references use atgc aswell to express uncertainty ... :|
                            char cRef = toupper(sRefRegion[j]);

                            if ((cRead <= 64) || (cRead >= 91))
                            {
                                std::cerr << "error character in read " << sReadName << " " << cRead << " -> " << cRef << std::endl;
                            }

                            if ((cRef <= 64) || (cRef >= 91))
                            {
                                std::cerr << "error character in ref " << sReadName << " " << cRead << " -> " << cRef << std::endl;
                            }


                            if (cRead == cRef)
                            {
                                ++equals;
                                ++iMatches;
                                ++iReadMatched;

                                if (mismatched > 0)
                                {
                                    mismatched = 0;
                                }

                            } else {
                                ++iMisMatches;

                                std::pair<char, char> subst(cRead, cRef);
                                pReadStats->ntSubst.add( subst, 1 );

                                if (equals > 0)
                                {

                                    if (j >= 5)
                                    {
                                        std::string sKmerBeforeMM = sReadRegion.substr(j-5, 5);
                                        pReadStats->mm2kmer.add(pReadStats->exactMatchCIGAR, sKmerBeforeMM); // TODO left this out
                                    }

                                    if (equals >= 21)
                                    {
                                        std::string sExactRead = sReadRegion.substr(j-equals, equals);
                                        std::vector<std::string> sKmers = this->calcKmers(sExactRead, 21);

                                        for (size_t k = 0; k < sKmers.size(); ++k)
                                        {
                                            if (sKmers.at(k).size() != 21)
                                            {
                                                std::cerr << sKmers.at(k) << " error length 32" << std::endl;
                                            }
                                            pReadStats->perfKmers.add(sKmers.at(k), 1);
                                        }
                                    }
                                }
                            }

                        }

                        //assert(iMatches+iMisMatches == pMRegion->getCIGARLength());

                        pReadStats->cigar2len.add(pReadStats->exactMatchCIGAR, iMatches);
                        pReadStats->cigar2len.add(pReadStats->mismatchCIGAR, iMisMatches);

                        iReadMatches += iMatches;
                        iReadMismatches += iMisMatches;

                        //std::cout << "to python" << pReadStats->cigar2len.keyString() << std::endl;


                        /*
                        //#pragma omp critical
                        {
                            std::cout << "Matches " << iMatches << " Mismatches " << iMisMatches << " Cigar Len"
                                      << pMRegion->getCIGARLength() << std::endl;
                        }

                        */

                    }



                }


            }

            pReadStats->iLongestMatched = iLongestMatched;
            pReadStats->fReadGCContent = (float)iReadGCBases / (float)iMatchedRegionsBases;
            pReadStats->fRefGCContent = (float)iRefGCBases / (float)iMatchedRegionsBases;

            pReadStats->fReadIdentity = (float) iReadMatched / (float) pReadStats->iReadLength;
            pReadStats->fRefIdentity = (float) iReadMatched / (float) pReadStats->iRefLength;

            /*
            if (pReadStats->sReadID == "ce27fd97-4547-4096-8d11-29b090889111_Basecall_1D_template")
            {
                std::cout << pReadStats->sReadID << " M " << iReadM << std::endl;
                std::cout << pReadStats->sReadID << " E " << iReadMatches << std::endl;
                std::cout << pReadStats->sReadID << " Z " << iReadMismatches << std::endl;

                uint32_t icount = 0;

                std::vector<uint32_t>*pvec = pReadStats->cigar2len.getVector('M');
                for (size_t j = 0; j < pvec->size(); ++j)
                {
                    icount += pvec->at(j);
                }

                std::cout << icount << std::endl;

                icount = 0;

                pvec = pReadStats->cigar2len.getVector('E');
                for (size_t j = 0; j < pvec->size(); ++j)
                {
                    icount += pvec->at(j);
                }

                std::cout << icount << std::endl;

                icount = 0;

                pvec = pReadStats->cigar2len.getVector('Z');
                for (size_t j = 0; j < pvec->size(); ++j)
                {
                    icount += pvec->at(j);
                }

                std::cout << icount << std::endl;
            }
            */

            PythonReadStats oStatsObj = pReadStats->toPython();
            std::vector<PythonReadStats>* pReads = (std::vector<PythonReadStats>*) pData;
            pReads->push_back(oStatsObj);

        } else {

            PythonReadStats oStatsObj = pReadStats->toPython();
            std::vector<PythonReadStats>* pReads = (std::vector<PythonReadStats>*) pData;
            pReads->push_back(oStatsObj);
        }

        delete pReadStats;
    }



private:

    std::vector<std::string> calcKmers(std::string& sSeq, uint32_t iK)
    {
        std::vector<std::string> oRet;

        if (iK > sSeq.size())
        {
            return oRet;
        }

        oRet.reserve(sSeq.size()-iK);

        for (size_t i = 0; i < sSeq.size()-iK; ++i)
        {
            oRet.push_back(sSeq.substr(i, iK));
        }

        return oRet;
    }


};