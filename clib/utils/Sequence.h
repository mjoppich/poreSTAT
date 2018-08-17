//
// Created by joppich on 12/1/15.
//

#ifndef PROJECT_SEQUENCE_H
#define PROJECT_SEQUENCE_H

#include <string>
#include <algorithm>

class Sequence : public std::string {

public:

    Sequence(std::string::iterator oBegin, std::string::iterator oEnd)
            : std::string(oBegin, oEnd)
    {

    }

    Sequence(std::string& oString, size_t iPos, size_t iLen = std::string::npos)
            : std::string(oString, iPos, iLen)
    {

    }

    Sequence(std::string* pSeq)
     : std::string(*pSeq)
    {

    }

    Sequence(std::string sSeq)
            : std::string(sSeq)
    {

    }

    Sequence(const char* seq) :
        std::string(seq) {
    }

    Sequence()
        :std::string()
    {

    }

    static char complement(char cIn)
    {

        switch (cIn)
        {
            case 'A': return 'T';
            case 'T': return 'A';
            case 'U': return 'A';

            case 'C': return 'G';
            case 'G': return 'C';

        }

        // just return whatever you got
        return cIn;

    }

    const Sequence& self() {return *this;};

    void toUpper()
    {
        std::transform(this->begin(), this->end(), this->begin(), [](unsigned char c) { return std::toupper(c); });
    }

    void replaceN();


    void reverse()
    {

        std::reverse(this->begin(), this->end());

    }

    void complement()
    {

        this->toUpper();

        for (size_t i = 0; i < this->length(); ++i)
        {
            (*this)[i] = Sequence::complement((*this)[i]);
        }


    }


    void reverseComplement()
    {

        this->reverse();
        this->complement();

    }

    Sequence prefix(size_t iLength)
    {
        return this->substr(0, iLength);
    }

    Sequence suffix( size_t iLength )
    {
        size_t iStringLength = this->size();

        return this->substr( iStringLength-iLength , iLength);
    }


    bool startsWith(const std::string& sIsPrefix)
    {

        Sequence sPrefix = this->prefix(sIsPrefix.size());

        return (sPrefix.compare(sIsPrefix) == 0);


    }

    bool endsWith(const std::string& sIsSuffix)
    {

        Sequence sSuffix = this->suffix( sIsSuffix.size() );

        return (sSuffix.compare(sIsSuffix) == 0);

    }

};


#endif //PROJECT_SEQUENCE_H
