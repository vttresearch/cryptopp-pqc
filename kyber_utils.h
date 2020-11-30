/*
 * Utility functions for Kyber. Adapted from the reference
 * implementation of Kyber by the CRYSTALS team
 * (https://github.com/pq-crystals/kyber)
 */

#ifndef CRYPTOPP_KYBER_UTILS_H
#define CRYPTOPP_KYBER_UTILS_H

#include "cryptlib.h"

NAMESPACE_BEGIN(CryptoPP)


class RandomBytes {
public:
    void randombytes(uint8_t *out, size_t outlen);
};



NAMESPACE_END

#endif