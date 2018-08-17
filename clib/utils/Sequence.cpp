//
// Created by joppich on 12/1/15.
//

#include "Sequence.h"
#include "SequenceUtils.h"

void Sequence::replaceN()
{

    for (size_t i = 0; i < this->size(); ++i)
    {

        if (this->self()[i] == 'N')
        {

            uint8_t iValue = ((i * this->size() * this->size()) % 4);

            this->replace(i, 1,1, SequenceUtils::unpack(iValue) );

        }

    }

}