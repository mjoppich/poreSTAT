/**
    Splitter Application Package - some toolkit to work with gtf/gff files
    Copyright (C) 2015  Markus Joppich

    The Splitter Application Package is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The Splitter Application Package is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "GenomicRegionUtils.h"


void GenomicRegionUtils::printCoverage(std::vector<uint32_t> *pCoverage) {


    std::stringstream oSS;

    for (size_t i = 0; i < pCoverage->size(); ++i)
    {

        if (pCoverage->at(i) > 0)
        {
            oSS << "1";
        } else {
            oSS << "0";
        }

    }

    std::cout << oSS.str() << std::endl;


}
