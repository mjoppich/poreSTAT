//
// Created by joppich on 12/1/15.
//

#ifndef PROJECT_SEQUENCEUTILS_H
#define PROJECT_SEQUENCEUTILS_H

#include "Sequence.h"
#include <inttypes.h>
#include <limits>
#include <iostream>
#include <vector>
#include <numeric> 

class SequenceUtils {
public:
    static uint8_t pack(char cBase)
    {
        switch (cBase)
        {
            case 'A': case 'a': return 0;

            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;

            case 'T': case 't': return 3;
            case 'U': case 'u': return 3;


            default: return rand() % 4; // keine hochwertige Zufallszahl, aber besser als ein Wert größer als 3
        }
    }

    static char unpack(uint8_t iBase)
    {
        switch (iBase)
        {
            case 0: return 'A';
            case 1: return 'C';
            case 2: return 'G';
            case 3: return 'T';

            default: return '.';

        }

        return 'x';
    }

    static uint8_t* packNedSequence(std::string *pString, std::string *pId, uint8_t* pOut = NULL)
    {
        uint32_t iSize = pString->size();

        size_t iMallocSize = (sizeof(uint8_t)*iSize + 3) / 4;

        bool allocCalled = false;
        if (pOut == NULL) {
            pOut = (uint8_t*)calloc(iMallocSize , sizeof(uint8_t));
            allocCalled = true;
        }

        uint8_t* pData = (uint8_t*)pString->c_str();

        uint8_t cVal = 0;
        for (size_t i = 0; i < iSize; ++i)
        {

            switch (pData[i])
            {
                case 'A': case 'a':
                    // 			std::cout << "A";
                    cVal = 0;
                    break;
                case 'C': case 'c':
                    // 			std::cout << "C";
                    cVal = 1;
                    break;

                case 'G': case 'g':
                    // 			std::cout << "G";
                    cVal = 2;
                    break;

                case 'T': case 't':  case 'U':  case 'u':
                    // 			std::cout << "T";
                    cVal = 3;
                    break;

                case 'N':  case 'n':
                    cVal = (i * pId->at( i % pId->size() ) ) % 4; // insert a random base in the hope that will be corrected
                    break;

                default:
                    pData = NULL;
                    i = iSize;
                    std::cerr << "non handled character" << std::endl;
                    break;
            }

            uint8_t iI8 = (uint8_t) i;
            cVal <<= (2 * (3 - (iI8 & 0x3)));
            pOut[i / 4] |= cVal;
        }

        if (pData == NULL)
        {
            if (allocCalled) {
                free(pOut);
            }
            return NULL;
        }

        return pOut;
    }

    static char* pAminoAcids;

    static char* cCodonToAA()
    {

        char* pAA = (char*) malloc(sizeof(char) * 64);


        // ALL A
        // A
        pAA[0] = 'K'; // A
        pAA[1] = 'N'; // C
        pAA[2] = 'K'; // G
        pAA[3] = 'N'; // T
        //C
        pAA[4] = 'T'; // A
        pAA[5] = 'T'; // C
        pAA[6] = 'T'; // G
        pAA[7] = 'T'; // T
        //G
        pAA[8] = 'R'; // A
        pAA[9] = 'S'; // C
        pAA[10] = 'R'; // G
        pAA[11] = 'S'; // T
        //T
        pAA[12] = 'I'; // A
        pAA[13] = 'I'; // C
        pAA[14] = 'M'; // G // AUG/ATG START CODON
        pAA[15] = 'I'; // T


        // ALL C
        // A
        pAA[16] = 'Q'; // A
        pAA[17] = 'H'; // C
        pAA[18] = 'Q'; // G
        pAA[19] = 'H'; // T
        // C
        pAA[20] = 'P'; // A
        pAA[21] = 'P'; // C
        pAA[22] = 'P'; // G
        pAA[23] = 'P'; // T
        // G
        pAA[24] = 'R'; // A
        pAA[25] = 'R'; // C
        pAA[26] = 'R'; // G
        pAA[27] = 'R'; // T
        // T
        pAA[28] = 'L'; // A
        pAA[29] = 'L'; // C
        pAA[30] = 'L'; // G
        pAA[31] = 'L'; // T


        // ALL G
        // A
        pAA[32] = 'E'; // A
        pAA[33] = 'D'; // C
        pAA[34] = 'E'; // G
        pAA[35] = 'D'; // T
        // C
        pAA[36] = 'A'; // A
        pAA[37] = 'A'; // C
        pAA[38] = 'A'; // G
        pAA[39] = 'A'; // T
        // G
        pAA[40] = 'G'; // A
        pAA[41] = 'G'; // C
        pAA[42] = 'G'; // G
        pAA[43] = 'G'; // T
        // T
        pAA[44] = 'V'; // A
        pAA[45] = 'V'; // C
        pAA[46] = 'V'; // G
        pAA[47] = 'V'; // T

        // ALL T/U
        // A
        pAA[48] = '.'; // A
        pAA[49] = 'Y'; // C
        pAA[50] = '.'; // G
        pAA[51] = 'Y'; // T
        // C
        pAA[52] = 'S'; // A
        pAA[53] = 'S'; // C
        pAA[54] = 'S'; // G
        pAA[55] = 'S'; // T
        // G
        pAA[56] = '.'; // A
        pAA[57] = 'C'; // C
        pAA[58] = 'W'; // G
        pAA[59] = 'C'; // T
        // T
        pAA[60] = 'L'; // A
        pAA[61] = 'F'; // C
        pAA[62] = 'L'; // G
        pAA[63] = 'F'; // T

        return pAA;

    }

    static Sequence translate(Sequence& oSequence, bool bChecked = true)
    {

        Sequence oReturnSeq;

        if ((bChecked) && (oSequence.size() % 3 != 0))
        {
            throw "Invalid Sequence length for translation";
        }

        size_t iAASeqLength = oSequence.size() / 3;
        oReturnSeq.resize(iAASeqLength, '.');

        for (size_t i = 0; i < oSequence.size()-2; i += 3)
        {

            char cA1 = oSequence.at(i);
            char cA2 = oSequence.at(i+1);
            char cA3 = oSequence.at(i+2);

            char cAA = SequenceUtils::translate( cA1, cA2, cA3 );

            oReturnSeq[i/3] = cAA;
        }


        return oReturnSeq;
    }

    static char translate(char c1, char c2, char c3)
    {

        uint8_t iIndex = 16 * SequenceUtils::pack(c1) + 4 * SequenceUtils::pack(c2) + SequenceUtils::pack(c3);

        return SequenceUtils::pAminoAcids[iIndex];
    }

    static bool isStartCodon(char (&cCodon)[3])
    {
        return SequenceUtils::isStartCodon(cCodon[0], cCodon[1], cCodon[2]);
    }

    static bool isStartCodon(char c1, char c2, char c3)
    {
        uint8_t iIndex = 16 * SequenceUtils::pack(c1) + 4 * SequenceUtils::pack(c2) + SequenceUtils::pack(c3);

        return (iIndex == 14);
    }

    static bool isStopCodon(char (&cCodon)[3])
    {
        return SequenceUtils::isStopCodon(cCodon[0], cCodon[1], cCodon[2]);
    }

    static bool isStopCodon(char c1, char c2, char c3)
    {
        uint8_t iIndex = 16 * SequenceUtils::pack(c1) + 4 * SequenceUtils::pack(c2) + SequenceUtils::pack(c3);

        return (SequenceUtils::pAminoAcids[iIndex] == '.');
    }

    static bool isStartCodon(Sequence* pSequence, uint32_t iShift = 0)
    {

        if (pSequence->length() < iShift+3)
            return false;

        return isStartCodon(pSequence->at(iShift+0), pSequence->at(iShift+1), pSequence->at(iShift+2));

    }

    static bool isStopCodon(Sequence* pSequence, uint32_t iShift = 0)
    {

        if (pSequence->length() < iShift+3)
            return false;

        return isStopCodon(pSequence->at(iShift+0), pSequence->at(iShift+1), pSequence->at(iShift+2));

    }

    static int32_t KOZAKdiff(std::string sSequence1)
    {
        if (sSequence1.length() != 6)
            return std::numeric_limits<std::int32_t>::max();

        std::string sKozak = "AAAAAA"; // (A|U)A(A|C)(A|C)AA
        int32_t iScore = 6;
        int32_t iMaxScore = iScore;


        for (uint8_t i = 0; i < sSequence1.length(); ++i)
        {

            char cValue = sSequence1.at(i);
            if (i == 0)
            {
                if (cValue == 'T')
                    continue;
            }

            if ((i == 2) || (i == 3))
            {
                if (cValue == 'C')
                    continue;
            }

            if (cValue == sKozak.at(i))
                continue;

            --iScore;
        }


        return iMaxScore-iScore;

    }



    /**
     *
     * \brief assumes that s1 and s2 are of same length. summing the array then would return the hamming distance
     *
     * \param s1 string1
     * \param s2 string2
     * \return array of size (min(s1,s2)) containing the distance for each position
     */
    static uint8_t* hamming_distance_array(const std::string &s1, const std::string &s2)
    {

//        int iDistance = 0;
        size_t iMin = std::min(s1.size(), s2.size());

        uint8_t* pReturn = (uint8_t*) calloc( sizeof(uint8_t), iMin);

        for (size_t i = 0; i < iMin; ++i)
        {

            if (s1.at( i ) != s2.at( i ))
            {
                pReturn[i] = 1;
            } else {
                pReturn[i] = 0;
            }

        }

        return pReturn;
    }

    static int hamming_distance(const std::string &s1, const std::string &s2)
    {

        int iDistance = 0;
        size_t iMin = std::min(s1.size(), s2.size());

        for (size_t i = 0; i < iMin; ++i)
        {

            if (s1.at( i ) != s2.at( i ))
            {
                iDistance += 1;
            }

        }


        if (iMin < s1.size())
            iDistance += (s1.size() - iMin);
        else
            iDistance += (s2.size() - iMin);

        return iDistance;
    }


    // from https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C.2B.2B
    static int levenshtein_distance(const std::string &s1, const std::string &s2)
    {
        // To change the type this function manipulates and returns, change
        // the return type and the types of the two variables below.
        int s1len = s1.size();
        int s2len = s2.size();

        auto column_start = (decltype(s1len))1;

        auto column = new decltype(s1len)[s1len + 1];
        std::iota(column + column_start, column + s1len + 1, column_start);

        for (auto x = column_start; x <= s2len; x++) {
            column[0] = x;
            auto last_diagonal = x - column_start;
            for (auto y = column_start; y <= s1len; y++) {
                auto old_diagonal = column[y];
                auto possibilities = {
                        column[y] + 1,
                        column[y - 1] + 1,
                        last_diagonal + (s1[y - 1] == s2[x - 1]? 0 : 1)
                };
                column[y] = std::min(possibilities);
                last_diagonal = old_diagonal;
            }
        }
        auto result = column[s1len];
        delete[] column;
        return result;
    }

    static std::vector<Sequence> splitLength(std::string *pString, uint32_t iWidth) {

        std::vector<Sequence> vReturn;

        std::string sCopy(pString->c_str());

        while( sCopy.size() > iWidth)
        {
            std::string sStart = sCopy.substr(0, iWidth);
            vReturn.push_back(sStart);

            sCopy = sCopy.substr(iWidth, -1);

        }

        if (sCopy.size() > 0)
            vReturn.push_back(sCopy);

        return vReturn;
    }
};


#endif //PROJECT_SEQUENCEUTILS_H
