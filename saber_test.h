/** The header file for the tests of the post-quantum KEM algorithm SABER. 
Adapted from the reference
implementation of SABER by the SABER team 
(https://github.com/KULeuven-COSIC/SABER) */


#ifndef SABER_TEST_H
#define SABER_TEST_H

#include "cryptlib.h"

NAMESPACE_BEGIN(CryptoPP)
NAMESPACE_BEGIN(Test)

int RunSaberTests();


NAMESPACE_END
NAMESPACE_END

#endif