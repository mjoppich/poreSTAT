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

#ifndef PROJECT_UTILSUNARY_H
#define PROJECT_UTILSUNARY_H

#include <inttypes.h>
#include <stddef.h>
#include <vector>
#include <algorithm>

template<typename T1>
class UtilsUnary {
public:

    static std::vector<T1>* copyVec(std::vector<T1>* pOriginal)
    {
        if (pOriginal == NULL)
            return NULL;


        std::vector<T1>* pVec = new std::vector<T1>();

        pVec->reserve( pOriginal->size() );

        // TODO make this more efficient
        for (uint32_t i = 0; i < pOriginal->size(); ++i)
        {
            pVec->push_back( pOriginal->at(i) );
        }

        return pVec;
    }

    static bool contains(std::vector<T1>* pVec, T1 oElem)
    {
        return std::find(pVec->begin(), pVec->end(), oElem) != pVec->end();
    }

    static size_t find(std::vector<T1>* pVec, T1 oElem)
    {

        for (size_t i = 0; i < pVec->size(); ++i)
        {

            if (pVec->at(i) == oElem)
                return i;

        }


        return -1;

    }

    static std::vector<std::vector<T1> > getCombinations(std::vector<T1>* pCompleteSet, uint32_t iChoose)
    {

        std::vector<std::vector<T1> > vReturn;

        std::string bitmask(iChoose, 1); // K leading 1's
        bitmask.resize(pCompleteSet->size(), 0); // N-K trailing 0's

        // print integers and permute bitmask
        do {

            std::vector<T1> vCombination;

            for (int i = 0; i < iChoose; ++i) // [0..N-1] integers
            {
                if (bitmask[i])
                {
                    vCombination.push_back( pCompleteSet->at(i) );
                }

            }
            vReturn.push_back(vCombination);

        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));


        return vReturn;
    }

    static std::vector< std::vector<T1> > getAllCombinations(std::vector<T1>* pCompleteSet)
    {

        std::vector<std::vector<T1> > vReturn;

        // always insert empty set
        vReturn.push_back( std::vector<T1>() );

        const size_t iElems = pCompleteSet->size() + 1;
        for (size_t i = 1; i < iElems; ++i)
        {

            std::vector<std::vector<T1> > vReturned = UtilsUnary<T1>::getCombinations(pCompleteSet, i);
            vReturn.insert(vReturn.end(), vReturned.begin(), vReturned.end());

        }


        return vReturn;
    }

};


#endif //PROJECT_UTILSUNARY_H
