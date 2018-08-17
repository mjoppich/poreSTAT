//
// Created by joppich on 4/12/16.
//

#ifndef PROJECT_STRINGUTILS_H
#define PROJECT_STRINGUTILS_H

#include <string>
#include <vector>
#include <sstream>

class StringUtils
{
public:

    /**
     *
     * PREFIX AND SUFFIX
     *
     *
     */

    static std::string prefix(std::string& sString, size_t iLength)
    {
        return sString.substr(0, iLength);
    }

    static std::string suffix(std::string& sString, size_t iLength)
    {
        size_t iStringLength = sString.size();

        return sString.substr( iStringLength-iLength , iLength);
    }


    static bool startsWith(std::string& sString, std::string& sIsPrefix)
    {

        std::string sPrefix = StringUtils::prefix(sString, sIsPrefix.size());

        return (sPrefix.compare(sIsPrefix) == 0);


    }

    static bool endsWith(std::string& sString, std::string& sIsSuffix)
    {

        std::string sSuffix = StringUtils::suffix(sString, sIsSuffix.size());

        return (sSuffix.compare(sIsSuffix) == 0);

    }

    /**
     *
     *
     * SPLIT VECTOR
     *
     *
     */

    static std::vector<std::string> &split(const std::string &sString, char cDelim, std::vector<std::string>* pElems) {
        std::stringstream sStringStream( sString );
        std::string sItem;
        while (std::getline(sStringStream, sItem, cDelim)) {
            pElems->push_back(sItem);
        }
        return *pElems;
    }

    static std::vector<std::string> split(const std::string &sString, char cDelim) {
        std::vector<std::string> vElems;
        StringUtils::split(sString, cDelim, &vElems);
        return vElems;
    }

    static std::vector<std::string> split(std::string *pString, char cDelim) {
        std::vector<std::string> vElems;

        if (pString != NULL)
            StringUtils::split(*pString, cDelim, &vElems);

        return vElems;
    }

    static std::vector<std::string> split(std::string* pString, char cDelim, std::vector<std::string>* pElems)
    {
        StringUtils::split(*pString, cDelim, pElems);
        return *pElems;
    }


    /*
     *
     * JOIN VECTOR
     *
     */

    static std::string join(std::vector<std::string>* pVector, std::string sDelim)
    {

        std::string sReturn;

        for (size_t i = 0; i < pVector->size(); ++i)
        {
            if (i > 0)
                sReturn.append(sDelim);

            sReturn.append(pVector->at(i));
        }

        return sReturn;
    }

};


#endif //PROJECT_STRINGUTILS_H
