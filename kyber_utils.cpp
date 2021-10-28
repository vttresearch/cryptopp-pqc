/*
 * Miscellaneous utility functions for Kyber. Adapted from the reference
 * implementation of Kyber by the CRYSTALS team
 * (https://github.com/pq-crystals/kyber)
 */


#include "pch.h"
#include "kyber.h"
#include "osrng.h"


NAMESPACE_BEGIN(CryptoPP)


//Random bytes implementation for kyber
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::RandomBytes(byte *out, size_t outLen) {
  AutoSeededRandomPool rng;
  rng.GenerateBlock(out, outLen);
}


//Montgomery reduction; given a 32-bit integer a, computes
//16-bit integer congruent to a * R^-1 mod q, where R=2^16
//returns integer in {-q+1,...,q-1} congruent to a * R^-1 modulo q.
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
sword16 Kyber<T_K, T_C, T_C2, T_ETA1>::MontgomeryReduce(sword32 a)
{
  sword16 t;

  t = (sword16)a*QINV;
  t = (a - (sword32)t*Q) >> 16;
  return t;
}


//Barrett reduction; given a 16-bit integer a, computes
//centered representative congruent to a mod q in {-(q-1)/2,...,(q-1)/2}
//returns integer in {-(q-1)/2,...,(q-1)/2} congruent to a modulo q.
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
sword16 Kyber<T_K, T_C, T_C2, T_ETA1>::BarrettReduce(sword16 a) {
  sword16 t;
  const sword16 v = ((1 << 26) + Q/2)/Q;
  t  = ((sword32)v*a + (1<<25)) >> 26;
  t *= Q;
  return a - t;
}


//load 4 bytes into a 32-bit integer
//in little-endian order
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
word32 Kyber<T_K, T_C, T_C2, T_ETA1>::Load32LittleEndian(const byte x[4])
{
  word32 r;
  r  = (word32)x[0];
  r |= (word32)x[1] << 8;
  r |= (word32)x[2] << 16;
  r |= (word32)x[3] << 24;
  return r;
}

//load 3 bytes into a 32-bit integer
//in little-endian order.
//This function is only needed for Kyber-512
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
word32 Kyber<T_K, T_C, T_C2, T_ETA1>::Load24LittleEndian(const byte x[3])
{
  if (ETA1 != 3) {
    return -1;
  }
  word32 r;
  r  = (word32)x[0];
  r |= (word32)x[1] << 8;
  r |= (word32)x[2] << 16;
  return r;
}

//Given an array of uniformly random bytes, compute
// polynomial with coefficients distributed according to
//a centered binomial distribution with parameter ETA1 or ETA2
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::Cbd(poly *r, const byte *buf, const byte eta)
{
  word32 i,j;
  word32 t,d;
  sword16 a,b;
  if (eta == 2) {

    for(i=0;i<N/8;i++) {
      t  = Load32LittleEndian(buf+4*i);
      d  = t & 0x55555555;
      d += (t>>1) & 0x55555555;

      for(j=0;j<8;j++) {
        a = (d >> (4*j+0)) & 0x3;
        b = (d >> (4*j+2)) & 0x3;
        r->coeffs[8*i+j] = a - b;
      }
    }
  } else if (eta == 3) {
    for(i=0;i<N/4;i++) {
      t  = Load24LittleEndian(buf+3*i);
      d  = t & 0x00249249;
      d += (t>>1) & 0x00249249;
      d += (t>>2) & 0x00249249;

      for(j=0;j<4;j++) {
        a = (d >> (6*j+0)) & 0x7;
        b = (d >> (6*j+3)) & 0x7;
        r->coeffs[4*i+j] = a - b;
      }
    }
  }
}

// Copy len bytes from x to r if b is 1;
// don't modify x if b is 0. Requires b to be in {0,1};
// assumes two's complement representation of negative integers.
// Runs in constant time.

template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::Cmov(byte *r, const byte *x, size_t len, byte b)
{
  size_t i;

  b = -b;
  for(i=0;i<len;i++)
    r[i] ^= b & (r[i] ^ x[i]);
}


/*
* Absorb step of the SHAKE128 specialized for the Kyber context.
*
*              - const byte *seed: pointer to KYBER_SYMBYTES input to be absorbed into state
*              - byte x: additional byte of input
*              - byte y: additional byte of input
*/
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::KyberShake128Absorb(const byte *seed, byte x, byte y)
{
  byte extSeed[SYMBYTES+2];

  std::memcpy(extSeed, seed, SYMBYTES);
  extSeed[SYMBYTES+0] = x;
  extSeed[SYMBYTES+1] = y;

  Shake128AbsorbOnce(extSeed, SYMBYTES + 2);
}

/*
* Usage of SHAKE256 as a PRF, concatenates secret and public input
* and then generates outputLen bytes of SHAKE256 output
*
*              - byte *output: pointer to output
*              - size_t outputLen: number of requested output bytes
*              - const byte *key: pointer to the key (of length SYMBYTES)
*              - byte nonce: single-byte nonce (public PRF input)
*/
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::KyberShake256PRF(byte *output, size_t outputLen, const byte *key, byte nonce)
{
  byte extKey[SYMBYTES+1];

  std::memcpy(extKey, key, SYMBYTES);
  extKey[SYMBYTES] = nonce;

  Shake256(output, outputLen, extKey, SYMBYTES + 1);
}


template class Kyber<2, 320, 128, 3>;
template class Kyber<3, 320, 128, 2>;
template class Kyber<4, 352, 160, 2>;

NAMESPACE_END