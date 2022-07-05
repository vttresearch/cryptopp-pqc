/* The header file for PQC FIPS 202. Adapted from the reference
FIPS 202 implementation of kyber by the CRYSTALS team 
(https://github.com/pq-crystals/kyber fips202.c and fips202.h) by Julius Hekkala
ORIGINAL LICENSE OF FIPS202.c
 * Based on the public domain implementation in crypto_hash/keccakc512/simple/ from
 * http://bench.cr.yp.to/supercop.html by Ronny Van Keer and the public domain "TweetFips202"
 * implementation from https://twitter.com/tweetfips202 by Gilles Van Assche, Daniel J. Bernstein,
 * and Peter Schwabe */ 

#ifndef CRYPTOPP_PQC_FIPS202
#define CRYPTOPP_PQC_FIPS202

#include "cryptlib.h"
#include "secblock.h"

NAMESPACE_BEGIN(CryptoPP)

/**
 * Class for FIPS202 functionality of Keccak needed by PQC algorithms.
 * Created so that the identical code is not in multiple places.
 */
class PQCFips202
{

public:
    //Constructor
    PQCFips202() {KeccakInit();};

    void Shake128Init();    
    void Shake128Absorb(const byte *input, size_t inputLen);
    void Shake128Finalize();
    void Shake128Squeeze(byte *output, size_t outputLen);
    void Shake128AbsorbOnce(const byte *input, size_t inputLen);
    void Shake128SqueezeBlocks(byte *output, size_t nBlocks);

    void Shake256Init();
    void Shake256Absorb(const byte *input, size_t inputLen);
    void Shake256Finalize();
    void Shake256Squeeze(byte *output, size_t outputLen);
    void Shake256AbsorbOnce(const byte *input, size_t inputLen);
    void Shake256SqueezeBlocks(byte *output, size_t nBlocks);
        
    void Shake128(byte *output, size_t outputLen, const byte *input, size_t inputLen);
    void Shake256(byte *output, size_t outputLen, const byte *input, size_t inputLen);
    void Sha3_256(byte h[32], const byte *input, size_t inputLen);
    void Sha3_512(byte h[64], const byte *input, size_t inputLen);

    

    CRYPTOPP_CONSTANT(SHAKE128_RATE=168);
    CRYPTOPP_CONSTANT(SHAKE256_RATE=136);
    CRYPTOPP_CONSTANT(SHA3_256_RATE=136);
    CRYPTOPP_CONSTANT(SHA3_512_RATE=72);
private:

    FixedSizeSecBlock<word64, 25> mState;
   //word64 mState[25];
    word32 mPos;

    void KeccakInit();
    void KeccakAbsorb(word32 rate, const byte *input, size_t inputLen);
    void KeccakFinalize(word32 rate, byte p);
    void KeccakSqueeze(byte *output, size_t outputLen, word32 rate);
    void KeccakAbsorbOnce(word32 rate, const byte *input, size_t inputLen, byte p);
    void KeccakSqueezeBlocks(byte *output, size_t nBlocks, word32 rate);





};



NAMESPACE_END
#endif
