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

#ifndef PROJECT_VECTORIZEDHASHMAP_H
#define PROJECT_VECTORIZEDHASHMAP_H

#include <map>
#include <omp.h>
#include <vector>
#include <stdlib.h>
#include <functional>
#include <pthread.h>
#include <set>
#include <unordered_map>

template <typename K, typename V>
class VectorizedMap  {

public:

    VectorizedMap()
    {
        m_pMap = new std::unordered_map<K,std::vector<V>*>();

        pthread_mutex_init(&mLock, NULL);

    }

    void applyAll(std::function< void (K, V) > oFun)
    {
        typename std::unordered_map<K,std::vector<V>* >::iterator oIt = m_pMap->begin();

        for ( ; oIt != m_pMap->end(); ++oIt)
        {

            for (typename std::vector<V>::iterator oJt = oIt->second->begin(); oJt != oIt->second->end(); ++ oJt)
            {
                oFun( oIt->first, *oJt);
            }
        }
    }

    size_t vecInsert(K oKey, V oVal)
    {
        pthread_mutex_lock(&mLock);

        std::vector<V>* pElements = this->get(oKey);

        if (pElements == NULL)
        {
            pElements = new std::vector<V>();
            m_pMap->insert(m_pMap->end(), std::pair<K, std::vector<V>*>(oKey, pElements) );
        }

        pElements->push_back(oVal);
        pthread_mutex_unlock(&mLock);

        return pElements->size()-1;

    }

    typename std::unordered_map<K,std::vector<V>* >* getMap()
    {
        return m_pMap;
    }


    std::vector<V>* get(K oKey)
    {
        typename std::unordered_map<K, std::vector< V>* >::iterator oIt = m_pMap->find(oKey);

        std::vector<V>* pElements = NULL;
        if (oIt != m_pMap->end())
        {
            pElements = (*oIt).second;
        }

        return pElements;
    }


    std::vector<V>* pop(K oKey)
    {
        pthread_mutex_lock(&mLock);

        typename std::unordered_map<K, std::vector< V>* >::iterator oIt = m_pMap->find(oKey);

        std::vector<V>* pElements = NULL;
        if (oIt != m_pMap->end())
        {
            pElements = (*oIt).second;
            (*oIt).second = new std::vector<V>();
        }

        pthread_mutex_unlock(&mLock);


        return pElements;
    }


    typename std::unordered_map<K, std::vector< V>* >::iterator begin()
    {
        return m_pMap->begin();
    }

    typename std::unordered_map<K, std::vector< V>* >::iterator end()
    {
        return m_pMap->end();
    }

    bool contains(K oKey)
    {

        typename std::unordered_map<K, std::vector< V>* >::iterator oIt = m_pMap->find(oKey);

        return (oIt != m_pMap->end());

    }

    void pruneHashMap(std::function<int(K,std::vector<V>*)> oFun)
    {
        typename std::unordered_map<K,std::vector<V>* >::iterator oIt = m_pMap->begin();

        for ( ; oIt != m_pMap->end(); )
        {
            oFun( oIt->first, oIt->second );

            if (oIt->second != NULL)
                delete (oIt->second);

            oIt = m_pMap->erase(oIt);
        }


    }


    ~VectorizedMap()
    {
        delete m_pMap;
    }

    std::set<K> keySet()
    {

        std::set<K> vReturn;

        typename std::unordered_map<K,std::vector<V>* >::iterator oIt = m_pMap->begin();

        for ( ; oIt != m_pMap->end(); ++oIt)
        {

            vReturn.insert( oIt->first );

        }

        return vReturn;

    }

    void reserve(size_t iSize)
    {
        m_pMap->reserve(iSize);
    }

protected:
    typename std::unordered_map<K,std::vector<V>* >* m_pMap;

    pthread_mutex_t mLock;


};


#endif //PROJECT_VECTORIZEDHASHMAP_H
