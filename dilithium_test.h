/** The header file for the tests of CRYSTALS-Dilithium. 
Adapted by Julius Hekkala from the public-domain reference
implementation of Dilithium by the CRYSTALS team 
(https://github.com/pq-crystals/dilithium) */


#ifndef DILITHIUM_TEST_H
#define DILITHIUM_TEST_H

#include "cryptlib.h"

NAMESPACE_BEGIN(CryptoPP)
NAMESPACE_BEGIN(Test)

int RunDilithiumTests();


NAMESPACE_END
NAMESPACE_END

#endif