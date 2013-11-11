/*
    Copyright (C) 2011 Matthias Kretz <kretz@kde.org>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <Vc/cpuid.h>
#include <iostream>
#include <cstring>

const char *g_fun = 0;

#define printFeature(fun) \
    if (g_fun == 0 || 0 == std::strcmp(g_fun, #fun)) { \
        std::cout << #fun ": " << Vc::CpuId::fun() << '\n'; \
    }

int main(int argc, char **argv)
{
    if (argc == 2) {
        g_fun = argv[1];
    }
    Vc::CpuId::init();
    printFeature(cacheLineSize);
    printFeature(processorType);
    printFeature(processorFamily);
    printFeature(processorModel);
    printFeature(logicalProcessors);
    printFeature(isAmd);
    printFeature(isIntel);
    printFeature(hasSse3);
    printFeature(hasVmx);
    printFeature(hasSmx);
    printFeature(hasEist);
    printFeature(hasTm2);
    printFeature(hasSsse3);
    printFeature(hasPdcm);
    printFeature(hasSse41);
    printFeature(hasSse42);
    printFeature(hasAes);
    printFeature(hasOsxsave);
    printFeature(hasAvx);
    printFeature(hasFpu);
    printFeature(hasVme);
    printFeature(hasDe);
    printFeature(hasPse);
    printFeature(hasTsc);
    printFeature(hasMsr);
    printFeature(hasPae);
    printFeature(hasMtrr);
    printFeature(hasCmov);
    printFeature(hasClfsh);
    printFeature(hasAcpi);
    printFeature(hasMmx);
    printFeature(hasSse);
    printFeature(hasSse2);
    printFeature(hasHtt);
    printFeature(hasSse4a);
    printFeature(hasMisAlignSse);
    printFeature(hasAmdPrefetch);
    printFeature(hasXop);
    printFeature(hasFma4);
    printFeature(hasRdtscp);
    printFeature(has3DNow);
    printFeature(has3DNowExt);
    printFeature(L1Instruction);
    printFeature(L1Data);
    printFeature(L2Data);
    printFeature(L3Data);
    printFeature(L1Associativity);
    printFeature(L2Associativity);
    printFeature(L3Associativity);
    printFeature(L1InstructionLineSize);
    printFeature(L1DataLineSize);
    printFeature(L2DataLineSize);
    printFeature(L3DataLineSize);
    printFeature(prefetch);
    return 0;
}
