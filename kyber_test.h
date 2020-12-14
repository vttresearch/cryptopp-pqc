/** The header file for the tests of CRYSTALS-Kyber. 
Adapted from the reference
implementation of kyber by the CRYSTALS team 
(https://github.com/pq-crystals/kyber) */


#ifndef KYBER_TEST_H
#define KYBER_TEST_H

#include "cryptlib.h"

NAMESPACE_BEGIN(CryptoPP)
NAMESPACE_BEGIN(Test)

int RunKyberTests();


NAMESPACE_END
NAMESPACE_END

#endif