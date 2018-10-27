//
// Created by mjopp on 05/10/2018.
//

#ifndef FOO_READSTATS_H
#define FOO_READSTATS_H

#include <vector>
#include <map>
#include <inttypes.h>
#include <string.h>

template <typename K, typename V>
class DefDict
{
public:
    DefDict()
    {
    }

    size_t length(K key)
    {
        if (!this->hasKey(key))
        {
            return 0;
        }

        return this->mMap.at(key).size();
    }

    std::string keyString()
    {
        std::string sret = "";

        typename std::map<K, std::vector<V> >::iterator it = this->mMap.begin();

        while (it != this->mMap.end())
        {
            sret.push_back(it->first);
        }

        return sret;
    }

    std::vector<V>* getVector(K key)
    {
        if (!this->hasKey(key))
        {
            return NULL;
        }

        return &(this->mMap.at(key));
    }

    void add(K key, V val)
    {
        if (!this->hasKey(key))
        {
            this->mMap.insert(std::pair<K, std::vector<V> >(key, std::vector<V>()));
        }

        this->mMap.at(key).push_back(val);
    }

    bool hasKey(K key)
    {
        typename std::map<K, std::vector<V> >::iterator it = this->mMap.find(key);

        if (it != this->mMap.end())
        {
            return true;
        }
        return false;
    }

    typename std::map<K, std::vector<V> >::iterator begin()
    {
        return this->mMap.begin();
    }

    typename std::map<K, std::vector<V> >::iterator end()
    {
        return this->mMap.end();
    }

    size_t size()
    {
        return this->mMap.size();
    }

    ~DefDict()
    {
    }

protected:

    std::map<K, std::vector<V> > mMap;

};

template <typename K>
class Counter
{
public:
    Counter()
    {
    }

    void add(K key, uint32_t iCount)
    {
        if (!this->hasKey(key))
        {
            this->mMap.insert(std::pair<K, uint32_t >(key, 0));
        }

        typename std::map<K, uint32_t>::iterator it = this->mMap.find(key);

        if (it != this->mMap.end())
        {
            it->second += iCount;
        }

    }

    typename std::map<K,uint32_t >::iterator begin()
    {
        return this->mMap.begin();
    }

    typename std::map<K,uint32_t >::iterator end()
    {
        return this->mMap.end();
    }

    bool hasKey(K key)
    {
        typename std::map<K, uint32_t>::iterator it = this->mMap.find(key);

        if (it != this->mMap.end())
        {
            return true;
        }
        return false;
    }

    size_t size()
    {
        return this->mMap.size();
    }

    ~Counter()
    {
    }

protected:

    std::map<K, uint32_t > mMap;

};

struct CounterPair
{
    CounterPair(const char* pkey, uint32_t ival)
    {
        key = pkey;
        value = ival;
    }

    const char* key;
    uint32_t value;
};

struct SeqCoverage
{
    SeqCoverage(const char* pseq, uint32_t istart, uint32_t iend)
    {
        seq = pseq;
        start = istart;
        end = iend;
    }

    const char* seq;
    uint32_t start;
    uint32_t end;
};

struct NTSubstitution
{
    NTSubstitution(char cFrom, char cTo, uint32_t iCount)
    {
        from = cFrom;
        to = cTo;
        count = iCount;
    }

    char from;
    char to;
    uint32_t count;
};



struct PythonReadStats
{
    char* pReadID;
    bool aligned;

    uint32_t iReadLength;
    uint32_t iRefLength;

    uint32_t iAlignQual;
    float fSeqQual_min;
    float fSeqQual_median;
    float fSeqQual_max;

    uint32_t iLongestMatched;
    float fReadGCContent;
    float fRefGCContent;

    float fReadIdentity;
    float fRefIdentity;

    uint32_t PERF_KMER_COUNT;
    CounterPair* PERF_KMERS;

    uint32_t CIGAR_COUNT;
    char* CIGARS;
    uint32_t* CIGAR_VEC_LENGTHS;
    uint32_t** CIGAR_LENGTHS;

    uint32_t* MM2KMER_VEC_LENGTHS;
    char*** MM2KMER_LENGTHS;

    uint32_t COV_COUNT;
    SeqCoverage* COVERAGES;

    uint32_t NT_SUBST_COUNT;
    NTSubstitution* NT_SUBST;

};

#include <iostream>


class AlignedReadStats
{
public:
    AlignedReadStats()
    {
    }

    PythonReadStats toPython()
    {
        PythonReadStats oRet;

        oRet.pReadID = (char*) malloc(sizeof(char)*(this->sReadID.size()+1));
        strncpy(oRet.pReadID, this->sReadID.c_str(), this->sReadID.size());
        oRet.pReadID[this->sReadID.size()] = '\0';

        oRet.aligned = this->aligned;

        oRet.iReadLength = this->iReadLength;
        oRet.iRefLength = this->iRefLength;

        oRet.iAlignQual = this->iAlignQual;
        oRet.fSeqQual_min = this->fSeqQual_min;
        oRet.fSeqQual_median = this->fSeqQual_median;
        oRet.fSeqQual_max = this->fSeqQual_max;

        oRet.iLongestMatched = this->iLongestMatched;
        oRet.fReadGCContent = this->fReadGCContent;
        oRet.fRefGCContent = this->fRefGCContent;

        oRet. fReadIdentity = this->fReadIdentity;
        oRet. fRefIdentity = this->fRefIdentity;

        oRet.PERF_KMER_COUNT = (uint32_t) this->perfKmers.size();
        oRet.PERF_KMERS = (CounterPair*) malloc(sizeof(CounterPair) * oRet.PERF_KMER_COUNT);

        std::map<std::string, uint32_t>::iterator oPCit = this->perfKmers.begin();
        size_t idx = 0;
        for (; oPCit != this->perfKmers.end(); ++oPCit)
        {
            oRet.PERF_KMERS[idx++] = CounterPair( oPCit->first.c_str(), oPCit->second );
        }


        oRet.CIGAR_COUNT = (uint32_t) this->cigar2len.size();
        oRet.CIGARS = (char*) malloc(sizeof(char) * oRet.CIGAR_COUNT);

        oRet.CIGAR_VEC_LENGTHS = (uint32_t*) malloc(sizeof(uint32_t) * oRet.CIGAR_COUNT);
        oRet.CIGAR_LENGTHS = (uint32_t**) malloc(sizeof(uint32_t*) * oRet.CIGAR_COUNT);

        std::map<char, std::vector<uint32_t> >::iterator oCLit = this->cigar2len.begin();
        size_t iArrIdx = 0;

        for (; oCLit != this->cigar2len.end(); ++oCLit)
        {
            char cCIGAR = oCLit->first;
                oRet.CIGARS[iArrIdx] = cCIGAR;
            oRet.CIGAR_VEC_LENGTHS[iArrIdx] = (uint32_t) oCLit->second.size();
            oRet.CIGAR_LENGTHS[iArrIdx] = (uint32_t*) malloc(sizeof(uint32_t)*oCLit->second.size());

            /*
            if (this->sReadID == "ce27fd97-4547-4096-8d11-29b090889111_Basecall_1D_template")
            {
                if ((cCIGAR == 'M') || (cCIGAR == 'Z') || (cCIGAR == 'E'))
                {
                    uint64_t iCount = 0;
                    for (size_t i = 0; i < oCLit->second.size(); ++i)
                    {
                        iCount += oCLit->second[i];
                    }

                    std::cout << this->sReadID << " " << cCIGAR << " " << iCount << std::endl;
                }
            }
            */

            for (size_t i = 0; i < oCLit->second.size(); ++i)
            {
                oRet.CIGAR_LENGTHS[iArrIdx][i] = oCLit->second[i];
            }

            ++iArrIdx;
        }


        oRet.MM2KMER_VEC_LENGTHS = (uint32_t*) malloc(sizeof(uint32_t) * oRet.CIGAR_COUNT);
        oRet.MM2KMER_LENGTHS = (char***) malloc(sizeof(char**) * oRet.CIGAR_COUNT);
        for (size_t i=0; i < oRet.CIGAR_COUNT; ++i)
        {
            char cCIGAR = oRet.CIGARS[i];

            uint32_t iElemCount = (uint32_t) this->mm2kmer.length(cCIGAR);
            oRet.MM2KMER_VEC_LENGTHS[i] = iElemCount;
            oRet.MM2KMER_LENGTHS[i] = (char**) malloc(sizeof(char*)*iElemCount);

            std::vector<std::string>* pKmers = this->mm2kmer.getVector(cCIGAR);

            if (pKmers != NULL)
            {
                for (size_t j = 0; j < iElemCount; ++j)
                {
                    oRet.MM2KMER_LENGTHS[i][j] = strdup(pKmers->at(j).c_str());
                }
            }
        }


        oRet.COV_COUNT = (uint32_t) this->seqCoverages.size();
        oRet.COVERAGES = (SeqCoverage*) malloc(sizeof(SeqCoverage) * oRet.COV_COUNT);

        for (size_t i = 0; i != this->seqCoverages.size(); ++i)
        {
            oRet.COVERAGES[i] = this->seqCoverages[i];
        }

        oRet.NT_SUBST_COUNT = (uint32_t) this->ntSubst.size();
        oRet.NT_SUBST = (NTSubstitution*) malloc(sizeof(NTSubstitution) * oRet.NT_SUBST_COUNT);
        std::map<std::pair<char, char>, uint32_t>::iterator oNTit = this->ntSubst.begin();
        idx = 0;
        for (; oNTit != this->ntSubst.end(); ++oNTit)
        {
            oRet.NT_SUBST[idx++] = NTSubstitution( oNTit->first.first, oNTit->first.second, oNTit->second );
        }


        return oRet;
    }

    ~AlignedReadStats()
    {
        //delete this->cigar2len;
    }

    std::string sReadID;
    bool aligned;

    uint32_t iReadLength;
    uint32_t iRefLength;

    uint32_t iAlignQual;
    float fSeqQual_min;
    float fSeqQual_median;
    float fSeqQual_max;

    uint32_t iLongestMatched;
    float fReadGCContent;
    float fRefGCContent;

    float fReadIdentity;
    float fRefIdentity;


    DefDict<char, uint32_t> cigar2len;
    DefDict<char, std::string> mm2kmer;
    Counter<std::string> perfKmers;

    Counter<std::pair<char, char>> ntSubst;
    std::vector<SeqCoverage> seqCoverages;

    const char mismatchCIGAR = 'Z';
    const char exactMatchCIGAR = 'E';

};

#endif //FOO_READSTATS_H
