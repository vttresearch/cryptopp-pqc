/*
 * Miscellaneous utility functions for Dilithium. Adapted by Julius Hekkala
 * from the public-domain reference
 * implementation of Dilithium by the CRYSTALS team
 * (https://github.com/pq-crystals/dilithium)
 */

#include "pch.h"
#include "dilithium.h"
#include "osrng.h"


NAMESPACE_BEGIN(CryptoPP)


/**
 * Get random bytes in a byte array the size of outLen.
 * Uses Crypto++'s AutoSeededRandomPool (osrng.h).
 */ 
void Dilithium::RandomBytes(byte *out, size_t outLen) {
  AutoSeededRandomPool rng;
  rng.GenerateBlock(out, outLen);
}

/* For finite field element a, compute a0, a1 such that
*  a mod^+ Q = a1*2^mD + a0 with -2^{mD-1} < a0 <= 2^{mD-1}.
*  Assumes a to be standard representative.
*  Arguments:  - sword32 a: input element
*              - sword32 *a0: pointer to output element a0
*
* Returns a1.
*/
sword32 Dilithium::Power2Round(sword32 *a0, sword32 a)  {
  sword32 a1;

  a1 = (a + (1 << (mD-1)) - 1) >> mD;
  *a0 = a - (a1 << mD);
  return a1;
}

/*  For finite field element a, compute high and low bits a0, a1 such
*   that a mod^+ mQ = a1*mAlpha + a0 with -mAlpha/2 < a0 <= mAlpha/2 except
*   if a1 = (mQ-1)/mAlpha where we set a1 = 0 and
*   -ALPHA/2 <= a0 = a^+ mod Q - Q < 0. Assumes a to be standard
*   representative.
*   Arguments: - sword32 a: input element
*              - sword32 *a0: pointer to output element a0
*
* Returns a1.
*/
sword32 Dilithium::Decompose(sword32 *a0, sword32 a) {
  sword32 a1;

  a1  = (a + 127) >> 7;
  if (mGamma2 == (mQ-1)/32) {
    a1  = (a1*1025 + (1 << 21)) >> 22;
    a1 &= 15;
  }
  
  else if (mGamma2 == (mQ-1)/88) {
    a1  = (a1*11275 + (1 << 23)) >> 24;
    a1 ^= ((43 - a1) >> 31) & a1;
  }

  *a0  = a - a1*2*mGamma2;
  *a0 -= (((mQ-1)/2 - *a0) >> 31) & mQ;
  return a1;
}

/* Compute hint bit indicating whether the low bits of the
*  input element overflow into the high bits. 
*  Arguments:   - sword32 a0: low bits of input element
*               - sword32 a1: high bits of input element
*
* Returns 1 if overflow.
*/
word32 Dilithium::MakeHint(sword32 a0, sword32 a1) {
  if(a0 > mGamma2 || a0 < - mGamma2 || (a0 == - mGamma2 && a1 != 0))
    return 1;

  return 0;
}

/* Correct high bits according to hint.
* Arguments:   - sword32 a: input element
*              - word32 hint: hint bit
*
* Returns corrected high bits.
*/
sword32 Dilithium::UseHint(sword32 a, word32 hint) {
  sword32 a0, a1;

  a1 = Decompose(&a0, a);
  if(hint == 0)
    return a1;

  if (mGamma2 == (mQ-1)/32) {
    if(a0 > 0)
      return (a1 + 1) & 15;
    else
      return (a1 - 1) & 15;
  }
  
  else if (mGamma2 == (mQ-1)/88) {
    if(a0 > 0)
      return (a1 == 43) ?  0 : a1 + 1;
    else
      return (a1 ==  0) ? 43 : a1 - 1;
  }
}

/*
* Bit-pack public key pk = (rho, t1).
* Arguments:   - byte *pk: pointer to output byte array
*              - const byte *rho: pointer to byte array containing rho
*              - const polyvec *t1: pointer to vector t1
*/
void Dilithium::PackPk(byte *pk, const byte *rho, const polyvec *t1)
{
  word32 i;

  for(i = 0; i < mSeedBytes; ++i)
    pk[i] = rho[i];
  pk += mSeedBytes;

  for(i = 0; i < mK; ++i)
    Polyt1Pack(pk + i*mPolyt1PackedBytes, &t1->at(i));
}

/*
* Unpack public key pk = (rho, t1).
* Arguments:   - const byte *rho: pointer to output byte array for rho
*              - const polyvec *t1: pointer to output vector t1
*              - byte *pk: byte array containing bit-packed pk
*/
void Dilithium::UnpackPk(byte *rho, polyvec *t1, const byte *pk)
{
  unsigned int i;

  for(i = 0; i < mSeedBytes; ++i)
    rho[i] = pk[i];
  pk += mSeedBytes;

  for(i = 0; i < mK; ++i)
    Polyt1Unpack(&t1->at(i), pk + i*mPolyt1PackedBytes);
}

/*
* Bit-pack secret key sk = (rho, tr, key, t0, s1, s2).
*
* Arguments:   - byte *sk: pointer to output byte array
*              - const byte *rho: pointer to byte array containing rho
*              - const byte *tr: byte array containing tr
*              - const byte *key: pointer to byte array containing key
*              - const polyvec *t0: pointer to vector t0
*              - const polyvec *s1: pointer to vector s1
*              - const polyvec *s2: pointer to vector s2
*/
void Dilithium::PackSk(byte *sk, const byte *rho, const byte *tr, const byte *key, const polyvec *t0,
             const polyvec *s1, const polyvec *s2)
{
  word32 i;

  for(i = 0; i < mSeedBytes; ++i)
    sk[i] = rho[i];
  sk += mSeedBytes;

  for(i = 0; i < mSeedBytes; ++i)
    sk[i] = key[i];
  sk += mSeedBytes;

  for(i = 0; i < mSeedBytes; ++i)
    sk[i] = tr[i];
  sk += mSeedBytes;

  for(i = 0; i < mL; ++i)
    PolyEtaPack(sk + i*mPolyEtaPackedBytes, &s1->at(i));
  sk += mL*mPolyEtaPackedBytes;

  for(i = 0; i < mK; ++i)
    PolyEtaPack(sk + i*mPolyEtaPackedBytes, &s2->at(i));
  sk += mK*mPolyEtaPackedBytes;

  for(i = 0; i < mK; ++i)
    Polyt0Pack(sk + i*mPolyt0PackedBytes, &t0->at(i));
}

/* Unpack secret key sk = (rho, tr, key, t0, s1, s2).
* Arguments:   - const byte *rho: pointer to output byte array for rho
*              - const byte *tr: pointer to output byte array for tr
*              - const byte *key: pointer to output byte array for key
*              - const polyvec *t0: pointer to output vector t0
*              - const polyvec *s1: pointer to output vector s1
*              - const polyvec *s2: pointer to output vector s2
*              - byte *sk: byte array containing bit-packed sk
*/
void Dilithium::UnpackSk(byte *rho, byte *tr, byte *key, polyvec *t0,  
                          polyvec *s1, polyvec *s2, const byte *sk)
{
  word32 i;

  for(i = 0; i < mSeedBytes; ++i)
    rho[i] = sk[i];
  sk += mSeedBytes;

  for(i = 0; i < mSeedBytes; ++i)
    key[i] = sk[i];
  sk += mSeedBytes;

  for(i = 0; i < mSeedBytes; ++i)
    tr[i] = sk[i];
  sk += mSeedBytes;
  
  for(i=0; i < mL; ++i)
    PolyEtaUnpack(&s1->at(i), sk + i*mPolyEtaPackedBytes);
  sk += mL*mPolyEtaPackedBytes;

  for(i=0; i < mK; ++i)
    PolyEtaUnpack(&s2->at(i), sk + i*mPolyEtaPackedBytes);
  sk += mK*mPolyEtaPackedBytes;
  
  for(i=0; i < mK; ++i)
    Polyt0Unpack(&t0->at(i), sk + i*mPolyt0PackedBytes);
}

/* Bit-pack signature sig = (c, z, h).
* Arguments:   - byte *sig: pointer to output byte array
*              - const byte *c: pointer to challenge hash 
*              - const polyvec *z: pointer to vector z
*              - const polyvec *h: pointer to hint vector h

*/
void Dilithium::PackSig(byte *sig, const byte *c, const polyvec *z, const polyvec *h)
{
  word32 i, j, k;

  for(i=0; i < mSeedBytes; ++i)
    sig[i] = c[i];
  sig += mSeedBytes;

  for(i = 0; i < mL; ++i)
    PolyzPack(sig + i*mPolyzPackedBytes, &z->at(i));
  sig += mL*mPolyzPackedBytes;

  /* Encode h */
  for(i = 0; i < mOmega + mK; ++i)
    sig[i] = 0;

  k = 0;
  for(i = 0; i < mK; ++i) {
    for(j = 0; j < mN; ++j)
      if(h->at(i).at(j) != 0)
        sig[k++] = j;

    sig[mOmega + i] = k;
  }
}


/* Unpack signature sig = (c, h, z).
* Arguments: 
*              - byte *c: pointer to output challenge hash
                - polyvec *z: pointer to output vector z
*              - polyvec *h: pointer to output hint vector h

*              - const byte *sig: pointer to byte array containing
*                bit-packed signature
*
* Returns 1 in case of malformed signature; otherwise 0.
*/
int Dilithium::UnpackSig(byte *c, polyvec *z, polyvec *h, const byte *sig)
{
  word32 i, j, k;
  
  for(i = 0; i < mSeedBytes; ++i)
    c[i] = sig[i];
  sig += mSeedBytes;

  for(i = 0; i < mL; ++i)
    PolyzUnpack(&z->at(i), sig + i*mPolyzPackedBytes);
  sig += mL*mPolyzPackedBytes;

  /* Decode h */
  k = 0;
  for(i = 0; i < mK; ++i) {
    for(j = 0; j < mN; ++j)
      h->at(i).at(j) = 0;

    if(sig[mOmega + i] < k || sig[mOmega + i] > mOmega)
      return 1;

    for(j = k; j < sig[mOmega + i]; ++j) {
      /* Coefficients are ordered for strong unforgeability */
      if(j > k && sig[j] <= sig[j-1]) return 1;
      h->at(i).at(sig[j]) = 1;
    }

    k = sig[mOmega + i];
  }

  /* Extra indices are zero for strong unforgeability */
  for(j = k; j < mOmega; ++j)
    if(sig[j])
      return 1;

  return 0;
}


/*
* For finite field element a with mQ*(-2)^31 <= a <= mQ*2^31,
* compute r \equiv a*2^{-32} (mod Q) such that -mQ <= r < mQ.
* Arguments:   - sword64: finite field element a
*
* Returns r.
*/
sword32 Dilithium::MontgomeryReduce(sword64 a) {
  sword32 t;

  t = (sword64)(sword32)a*mQInv;
  t = (a - (sword64)t*mQ) >> 32;

  return t;
}

/* 
* For finite field element a with a <= 2^{31} - 2^{22} - 1, compute r \equiv a (mod Q)
* such that -6283009 <= r <= 6283007.
* Arguments:   - sword32: finite field element a
*
* Returns r.
*/
sword32 Dilithium::Reduce32(sword32 a) {
  sword32 t;

  t = (a + (1 << 22)) >> 23;
  t = a - t*mQ;
  return t;
}

/* 
*  For finite field element a, compute standard
* representative r = a mod^+ Q.
* Arguments:   - sword32: finite field element a
*
* Returns r.
*/
sword32 Dilithium::Freeze(sword32 a) {
  a = Reduce32(a);
  a = CAddQ(a);
  return a;
}

/*
* Add Q if input coefficient is negative.
*
* Arguments:   - sword32: finite field element a
*
* Returns r.
*/
sword32 Dilithium::CAddQ(sword32 a) {
  a += (a >> 31) & mQ;
  return a;
}



/*
 * Initiate Shake-128 Stream.
 * Arguments: - *keccakState 
 * 
 * 
 */
void Dilithium::Shake128StreamInit(const byte *seed, word16 nonce)
{
  byte t[2];
  t[0] = nonce;
  t[1] = nonce >> 8;

  Shake128Init();
  Shake128Absorb(seed, mSeedBytes);
  Shake128Absorb(t, 2);
  Shake128Finalize();
}

void Dilithium::Shake256StreamInit(const byte *seed, word16 nonce)
{
  byte t[2];
  t[0] = nonce;
  t[1] = nonce >> 8;

  Shake256Init();
  Shake256Absorb(seed, mCrhBytes);
  Shake256Absorb(t, 2);
  Shake256Finalize();
}

NAMESPACE_END