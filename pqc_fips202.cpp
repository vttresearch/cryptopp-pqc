/* The C++ file for PQC FIPS 202. Adapted from the reference
implementation of kyber by the CRYSTALS team 
(https://github.com/pq-crystals/kyber fips202.c and fips202.h)
ORIGINAL LICENSE OF fips202.c
 * Based on the public domain implementation in crypto_hash/keccakc512/simple/ from
 * http://bench.cr.yp.to/supercop.html by Ronny Van Keer and the public domain "TweetFips202"
 * implementation from https://twitter.com/tweetfips202 by Gilles Van Assche, Daniel J. Bernstein,
 * and Peter Schwabe */ 

#include "pch.h"
#include "pqc_fips202.h"
#include <iostream>

NAMESPACE_BEGIN(CryptoPP)

// The Keccak core function
extern void KeccakF1600(word64 *state);

/*
* Load 8 bytes into word64 in little-endian order
*
*               - const byte *x: pointer to input byte array
*
* Returns the loaded 64-bit unsigned integer
*/
static word64 Load64(const byte x[8]) {
  word32 i;
  word64 r = 0;

  for(i=0;i<8;i++)
    r |= (word64)x[i] << 8*i;

  return r;
}

/*
* Store a 64-bit integer to array of 8 bytes in little-endian order
*
*              - byte *x: pointer to the output byte array (allocated)
*              - byte u: input 64-bit unsigned integer
*/
static void Store64(byte x[8], word64 u) {
  word32 i;

  for(i=0;i<8;i++)
    x[i] = u >> 8*i;
}


/*
* Initializes the Keccak state.
*
*              - keccakState *state: pointer to Keccak state
*/
void PQCFips202::KeccakInit()
{
  memset(mState, 0, mState.SizeInBytes());
  mPos = 0;
}

/*
* Absorb step of Keccak; incremental.
*
*              - word32 rate:   rate in bytes (e.g., 168 for SHAKE128)
*              - const byte *input: pointer to input to be absorbed into s
*              - size_t inputLen:      length of input in bytes
*
*/
void PQCFips202::KeccakAbsorb(word32 rate, const byte *input, size_t inputLen) {
  word32 i;

  while(mPos+inputLen >= rate) {
    for(i=mPos;i<rate;i++) {
      mState.data()[i/8] ^= (word64)*input++ << 8*(i%8);   
    }
      
    inputLen -= rate-mPos;
    KeccakF1600(mState);
    mPos = 0;
  }
  for(i=mPos;i<mPos+inputLen;i++) {
    mState.data()[i/8] ^= (word64)*input++ << 8*(i%8);
  }
  mPos = i;
}

/*
* Finalize absorb step.
*              - word32 rate: rate in bytes (e.g., 168 for SHAKE128)
*
*              - byte p: domain separation byte
*/
void PQCFips202::KeccakFinalize(word32 rate, byte p)
{
  mState.data()[mPos/8] ^= (word64)p << 8*(mPos%8);
  mState.data()[rate/8-1] ^= 1ULL << 63;
}


/*              
* Squeeze step of Keccak. Squeezes arbitratrily many bytes.
* Modifies the state. Can be called multiple times to keep
* squeezing, i.e., is incremental.
*
*              - byte *output: pointer to output
*              - size_t outputLen: number of bytes to be squeezed (written to out)
*              - word32 rate: rate in bytes (e.g., 168 for SHAKE128)
*
*
*/
void PQCFips202::KeccakSqueeze(byte *output, size_t outputLen, word32 rate)
{
  word32 i;

  while(outputLen) {
    if(mPos == rate) {
      KeccakF1600(mState);
      mPos = 0;
    }
    for(i=mPos;i < rate && i < mPos+outputLen; i++)
      *output++ = mState.data()[i/8] >> 8*(i%8);
    outputLen -= i-mPos;
    mPos = i;
  }
}

/*
* Absorb step of Keccak;
* non-incremental, starts by zeroeing the state.
*
*              - word32 rate: rate in bytes (e.g., 168 for SHAKE128)
*              - const byte *input: pointer to input to be absorbed into state
*              - size_t inputLen: length of input in bytes
*              - byte p: domain-separation byte for different Keccak-derived functions
*/
void PQCFips202::KeccakAbsorbOnce(word32 rate, const byte *input, size_t inputLen, byte p)
{
  word32 i;

  memset(mState, 0, mState.SizeInBytes());

  while(inputLen >= rate) {
    for(i=0;i<rate/8;i++)
      mState.data()[i] ^= Load64(input+8*i);
    input += rate;
    inputLen -= rate;
    KeccakF1600(mState);
  }
  for(i=0;i<inputLen;i++)
    mState.data()[i/8] ^= (word64)input[i] << 8*(i%8);
  mState.data()[i/8] ^= (word64)p << 8*(i%8);
  mState.data()[(rate-1)/8] ^= 1ULL << 63;

}

/*
* Squeeze step of Keccak. Squeezes full blocks of r bytes each.
* Modifies the state. Can be called multiple times to keep
* squeezing, i.e., is incremental. Assumes zero bytes of current
* block have already been squeezed.
*
*              - byte *output: pointer to output blocks
*              - size_t nBlocks: number of blocks to be squeezed (written to output)
*              - word32 rate: rate in bytes (e.g., 168 for SHAKE128)
**************************************************/
void PQCFips202::KeccakSqueezeBlocks(byte *output, size_t nBlocks, word32 rate)
{
  word32 i;

  while(nBlocks) {
    KeccakF1600(mState);
   for(i=0;i<rate/8;i++)
      Store64(output+8*i, mState.data()[i]);
    output += rate;
    nBlocks -= 1;
  }
}

/*
* Initilizes Keccak state for use as SHAKE128 XOF
*

*/
void PQCFips202::Shake128Init()
{
    KeccakInit();
}

/*
* Absorb step of the SHAKE128 XOF; incremental.
*
*              - const byte *input: pointer to input to be absorbed into state
*              - size_t inputLen: length of input in bytes
*/
void PQCFips202::Shake128Absorb(const byte *input, size_t inputLen)
{
  KeccakAbsorb(SHAKE128_RATE, input, inputLen);
}

/*
* Finalize absorb step of the SHAKE128 XOF.
*
*/
void PQCFips202::Shake128Finalize()
{
  KeccakFinalize(SHAKE128_RATE, 0x1F);
  mPos = SHAKE128_RATE;
}


/*
* Squeeze step of SHAKE128 XOF. Squeezes arbitraily many
* bytes. Can be called multiple times to keep squeezing.
*
*              - byte *output: pointer to output blocks
*              - size_t outputLen : number of bytes to be squeezed (written to output)
*/
void PQCFips202::Shake128Squeeze(byte *output, size_t outputLen)
{
  KeccakSqueeze(output, outputLen, SHAKE128_RATE);
}


/*
* Initialize, absorb into and finalize SHAKE128 XOF; non-incremental.
*
*              - const byte *input: pointer to input to be absorbed into state
*              - size_t inputLen: length of input in bytes
*/
void PQCFips202::Shake128AbsorbOnce(const byte *input, size_t inputLen)
{
  KeccakAbsorbOnce(SHAKE128_RATE, input, inputLen, 0x1F);
  mPos = SHAKE128_RATE;
}

/*
* Squeeze step of SHAKE128 XOF. Squeezes full blocks of
* SHAKE128_RATE bytes each. Can be called multiple times
* to keep squeezing. Assumes new block has not yet been
* started (mPos = SHAKE128_RATE).
*
*              - byte *output: pointer to output blocks
*              - size_t nBlocks: number of blocks to be squeezed (written to output)      
*/
void PQCFips202::Shake128SqueezeBlocks(byte *output, size_t nBlocks)
{
  KeccakSqueezeBlocks(output, nBlocks, SHAKE128_RATE);
}

/*
* Initilizes Keccak state for use as SHAKE256 XOF
*
*/
void PQCFips202::Shake256Init()
{
    KeccakInit();
}


/*
* Absorb step of the SHAKE256 XOF; incremental.
*
*              - const byte *input: pointer to input to be absorbed into state
*              - size_t inputLen: length of input in bytes
*/
void PQCFips202::Shake256Absorb(const byte *input, size_t inputLen)
{
  KeccakAbsorb(SHAKE256_RATE, input, inputLen);
}

/*
* Finalize absorb step of the SHAKE256 XOF.
*
*/
void PQCFips202::Shake256Finalize()
{
  KeccakFinalize(SHAKE256_RATE, 0x1F);
  mPos = SHAKE256_RATE;
}

/*
* Squeeze step of SHAKE256 XOF. Squeezes arbitraily many
* bytes. Can be called multiple times to keep squeezing.
*
*              - byte *output: pointer to output blocks
*              - size_t outputLen : number of bytes to be squeezed (written to output)
*/
void PQCFips202::Shake256Squeeze(byte *output, size_t outputLen)
{
  KeccakSqueeze(output, outputLen, SHAKE256_RATE);
}

/*
* Initialize, absorb into and finalize SHAKE256 XOF; non-incremental.
*
*              - const byte *input: pointer to input to be absorbed into state
*              - size_t inputLen: length of input in bytes
*/
void PQCFips202::Shake256AbsorbOnce(const byte *input, size_t inputLen)
{
  KeccakAbsorbOnce(SHAKE256_RATE, input, inputLen, 0x1F);
  mPos = SHAKE256_RATE;
}

/*
* Squeeze step of SHAKE256 XOF. Squeezes full blocks of
* SHAKE256_RATE bytes each. Can be called multiple times
* to keep squeezing. Assumes new block has not yet been
* started (mPos = SHAKE256_RATE).
*
*              - byte *output: pointer to output blocks
*              - size_t nBlocks: number of blocks to be squeezed (written to output)
*/
void PQCFips202::Shake256SqueezeBlocks(byte *output, size_t nBlocks)
{
  KeccakSqueezeBlocks(output, nBlocks, SHAKE256_RATE);
}

/*
* SHAKE128 XOF with non-incremental API
*
*              - byte *output: pointer to output
*              - size_t outputLen: requested output length in bytes
*              - const byte *input: pointer to input
*              - size_t inputLen: length of input in bytes
*/
void PQCFips202::Shake128(byte *output, size_t outputLen, const byte *input, size_t inputLen)
{
  size_t nBlocks;

  Shake128AbsorbOnce(input, inputLen);
  nBlocks = outputLen/SHAKE128_RATE;
  Shake128SqueezeBlocks(output, nBlocks);
  outputLen -= nBlocks*SHAKE128_RATE;
  output += nBlocks*SHAKE128_RATE;
  Shake128Squeeze(output, outputLen);
}

/*
* SHAKE256 XOF with non-incremental API
*
*              - byte *output: pointer to output
*              - size_t outputLen: requested output length in bytes
*              - const byte *input: pointer to input
*              - size_t inputLen: length of input in bytes
*/
void PQCFips202::Shake256(byte *output, size_t outputLen, const byte *input, size_t inputLen)
{
  size_t nBlocks;

  Shake256AbsorbOnce(input, inputLen);
  nBlocks = outputLen/SHAKE256_RATE;
  Shake256SqueezeBlocks(output, nBlocks);
  outputLen -= nBlocks*SHAKE256_RATE;
  output += nBlocks*SHAKE256_RATE;
  Shake256Squeeze(output, outputLen);
}

/*
* SHA3-256 with non-incremental API
*
*              - byte *h: pointer to output (32 bytes)
*              - const byte *input: pointer to input
*              - size_t inputLen: length of input in bytes
*/
void PQCFips202::Sha3_256(byte h[32], const byte *input, size_t inputLen)
{
  word32 i;

  KeccakAbsorbOnce(SHA3_256_RATE, input, inputLen, 0x06);
  KeccakF1600(mState);
  for(i=0;i<4;i++)
    Store64(h+8*i,mState.data()[i]);
}

/*
* SHA3-512 with non-incremental API
*
*              - byte *h: pointer to output (64 bytes)
*              - const byte *input: pointer to input
*              - size_t inputLen: length of input in bytes
*/
void PQCFips202::Sha3_512(byte h[64], const byte *input, size_t inputLen)
{
  word32 i;

  KeccakAbsorbOnce(SHA3_512_RATE, input, inputLen, 0x06);
  KeccakF1600(mState);
  for(i=0;i<8;i++)
    Store64(h+8*i,mState.data()[i]);

}




NAMESPACE_END