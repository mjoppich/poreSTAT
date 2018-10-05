#include <vector>
#include <map>
#include "../utils/FASTAreader.h"


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
        pReadStats->aligned = false;

        if (this->readIsMapped(pRead))
        {
            this->seenReads += 1;
            pReadStats->aligned = true;

            AlignedReadStats* pReadStats = new AlignedReadStats();
            
            std::string sAlignedSeqName = this->getSeqName(pRead->core.tid);

            pReadStats->iAlignQual = this->getAlignQuality(pRead);

            std::vector<float> vSeqQuals = this->getSeqQuality(pRead);
            pReadStats->fSeqQual_min = vSeqQuals[0];
            pReadStats->fSeqQual_median = vSeqQuals[1];
            pReadStats->fSeqQual_max = vSeqQuals[2];


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

            if (iRetRes)
            {
                std::string sReadSeq = this->getReadSequence(pRead);

                for (uint32_t i = 0; i < iRetRes; ++i)
                {
                    
                    ReadRegion* pMRegion = (ReadRegion*)allAlignedRegions.at(i);

                    pReadStats->cigar2len.add( pMRegion->getCIGAR(), pMRegion->getCIGARLength() );

                    if (pMRegion->getCIGAR() == 'M')
                    {

                        std::string sRefRegion;
                        std::string sReadRegion;

                        #pragma omp critical
                        {



                        std::cout << pMRegion->getStart() << " " << pMRegion->getLength() << " " << sAlignedSeqName << std::endl;
                        std::cout << pMRegion->pReadRegion->getStart() << " " << pMRegion->pReadRegion->getLength() << " " << this->getReadID(pRead) << std::endl;


                        sRefRegion = this->pFASTAReader->retrieveSequence(pMRegion, sAlignedSeqName);
                        sReadRegion = sReadSeq.substr(pMRegion->pReadRegion->getStart(), pMRegion->pReadRegion->getLength()-1);

                        std::cout << sRefRegion << std::endl;
                        std::cout << sReadRegion << std::endl;


                        };

                        // add coverage
                        pReadStats->seqCoverages.push_back(SeqCoverage(sAlignedSeqName.c_str(), pMRegion->getStart(), pMRegion->getEnd()));

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
                            if (sReadRegion[j] == sRefRegion[j])
                            {
                                ++equals;
                                ++iMatches;
                                ++iReadMatched;

                                if (mismatched > 0)
                                {
                                    pReadStats->cigar2len.add( pReadStats->mismatchCIGAR, mismatched );
                                    mismatched = 0;
                                }

                            } else {
                                ++iMisMatches;

                                std::pair<char, char> subst(sReadRegion[j], sRefRegion[j]);
                                pReadStats->ntSubst.add( subst, 1 );

                                if (equals > 0)
                                {
                                    pReadStats->cigar2len.add( pReadStats->exactMatchCIGAR, equals );

                                    if (j >= 5)
                                    {
                                        std::string sKmerBeforeMM = sReadRegion.substr(j-5, 5);
                                        pReadStats->mm2kmer.add(pReadStats->exactMatchCIGAR, sKmerBeforeMM);
                                    }

                                    if (equals >= 21)
                                    {
                                        std::string sExactRead = sReadRegion.substr(j-equals, equals);
                                        std::vector<std::string> sKmers = this->calcKmers(sExactRead, 21);

                                        for (size_t k = 0; k < sKmers.size(); ++k)
                                        {
                                            pReadStats->perfKmers.add(sKmers.at(k), 1);
                                        }
                                    }
                                }
                            }

                        }



                    }



                }


            }

            pReadStats->iLongestMatched = iLongestMatched;
            pReadStats->fReadGCContent = (float)iReadGCBases / (float)iMatchedRegionsBases;
            pReadStats->fRefGCContent = (float)iRefGCBases / (float)iMatchedRegionsBases;

            pReadStats->fReadIdentity = (float) iReadMatched / (float) pReadStats->iReadLength;
            pReadStats->fRefIdentity = (float) iReadMatched / (float) pReadStats->iRefLength;

            #pragma omp critical
            {
                std::vector<PythonReadStats>* pReads = (std::vector<PythonReadStats>*) pData;
                pReads->push_back(pReadStats->toPython());

            };

        }
    }



private:

    std::vector<std::string> calcKmers(std::string& sSeq, uint32_t iK)
    {
        std::vector<std::string> oRet;

        if (iK > sSeq.size())
        {
            return oRet;
        }

        for (size_t i = 0; i < sSeq.size()-iK; ++i)
        {
            oRet.push_back(sSeq.substr(i, iK));
        }

        return oRet;
    }


};