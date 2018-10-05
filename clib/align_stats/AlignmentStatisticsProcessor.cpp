#include "AlignmentStatisticProcessor.h"

/*
template <typename K, typename V> void DefDict<K,V>::add(K key, V val)
{
    if (!this->hasKey(key))
    {
        this->pMap->insert(std::pair<K, std::vector<V> >(key, std::vector<V>()));
    }

    this->pMap->at(key).push_back(val);
}


template <typename K, typename V> bool DefDict<K,V>::hasKey(K key)
{
    return false;
    if (this->pMap == NULL)
    {
        std::cerr << "NULL PMAP" << std::endl;
    }

    typename std::map<K, std::vector<V> >::iterator it = this->pMap->find(key);

    if (it != this->pMap->end())
    {
        return true;
    }
    return false;
}
*/