/*
 * Polynomial functions for Dilithium. Adapted by Julius Hekkala
 * from the public-domain reference
 * implementation of Dilithium by the CRYSTALS team
 * (https://github.com/pq-crystals/dilithium)
 */

#include "pch.h"
#include "dilithium.h"



NAMESPACE_BEGIN(CryptoPP)

/*
* Inplace reduction of all coefficients of polynomial to
* representative in [0,2*mQ[.
*
* Arguments:   - poly *a: pointer to input/output polynomial
*/
void Dilithium::PolyReduce(poly *a) {
  word32 i;

  for(i = 0; i < mN; ++i)
    a->at(i) = Reduce32(a->at(i));
}

/*
* For all coefficients of in/out polynomial subtract mQ if
* coefficient is bigger than mQ.
*
* Arguments:   - poly *a: pointer to input/output polynomial
*/
void Dilithium::PolyCsubQ(poly *a) {
  word32 i;
  

  for(i = 0; i < mN; ++i)
    a->at(i) = CSubQ(a->at(i));

}

/*
* Inplace reduction of all coefficients of polynomial to
* standard representatives.
*
* Arguments:   - poly *a: pointer to input/output polynomial
*/
void Dilithium::PolyFreeze(poly *a) {
  word32 i;


  for(i = 0; i < mN; ++i)
    a->at(i) = Freeze(a->at(i));

}

/*
* Add polynomials. No modular reduction is performed.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first summand
*              - const poly *b: pointer to second summand
*/
void Dilithium::PolyAdd(poly *c, const poly *a, const poly *b)  {
  word32 i;

  for(i = 0; i < mN; ++i)
    c->at(i) = a->at(i) + b->at(i);

}

/*
* Subtract polynomials. Assumes coefficients of second input
* polynomial to be less than 2*mQ. No modular reduction is
* performed.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial to be
*                               subtraced from first input polynomial
*/
void Dilithium::PolySub(poly *c, const poly *a, const poly *b) {
  word32 i;

  for(i = 0; i < mN; ++i)
    c->at(i) = a->at(i) + 2*mQ - b->at(i);

}

/*
* Multiply polynomial by 2^D without modular reduction. Assumes
* input coefficients to be less than 2^{32-D}.
*
* Arguments:   - poly *a: pointer to input/output polynomial
*/
void Dilithium::PolyShiftL(poly *a) {
  word32 i;

  for(i = 0; i < mN; ++i)
    a->at(i) <<= mD;

}

/*
* Inplace forward NTT. Output coefficients can be up to
*              16*mQ larger than input coefficients.
*
* Arguments:   - poly *a: pointer to input/output polynomial
*/
void Dilithium::PolyNtt(poly *a) {

  Ntt(a->data());

}

/*
* Inplace inverse NTT and multiplication by 2^{32}.
* Input coefficients need to be less than 2*mQ.
* Output coefficients are less than 2*mQ.
*
* Arguments:   - poly *a: pointer to input/output polynomial
*/
void Dilithium::PolyInvNttToMont(poly *a) {

  InvNttToMont(a->data());

}

/*
* Pointwise multiplication of polynomials in NTT domain
* representation and multiplication of resulting polynomial
* by 2^{-32}. Output coefficients are less than 2*mQ if input
* coefficient are less than 22*mQ.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
*/
void Dilithium::PolyPointwiseMontgomery(poly *c, const poly *a, const poly *b) {
  word32 i;

  for(i = 0; i < mN; ++i)
    c->at(i) = MontgomeryReduce((word64)a->at(i) * b->at(i));

}

/*
* For all coefficients c of the input polynomial,
* compute c0, c1 such that c mod mQ = c1*2^D + c0
* with -2^{D-1} < c0 <= 2^{D-1}. Assumes coefficients to be
* standard representatives.
*
* Arguments:   - poly *a1: pointer to output polynomial with coefficients c1
*              - poly *a0: pointer to output polynomial with coefficients mQ + c0
*              - const poly *v: pointer to input polynomial
*/
void Dilithium::PolyPower2Round(poly *a1, poly *a0, const poly *a) {
  word32 i;

  for(i = 0; i < mN; ++i)
    a1->at(i) = Power2Round(a->at(i), &a0->at(i));

}

/*
* For all coefficients c of the input polynomial,
* compute high and low bits c0, c1 such c mod mQ = c1*mAlpha + c0
* with -mAlpha/2 < c0 <= mAlpha/2 except c1 = (mQ-1)/mAlpha where we
* set c1 = 0 and -mAlpha/2 <= c0 = c mod mQ - mQ < 0.
* Assumes coefficients to be standard representatives.
*
* Arguments:   - poly *a1: pointer to output polynomial with coefficients c1
*              - poly *a0: pointer to output polynomial with coefficients mQ + c0
*              - const poly *c: pointer to input polynomial
*/
void Dilithium::PolyDecompose(poly *a1, poly *a0, const poly *a) {
  word32 i;

  for(i = 0; i < mN; ++i)
    a1->at(i) = Decompose(a->at(i), &a0->at(i));

}



/*
* Compute hint polynomial. The coefficients of which indicate
* whether the low bits of the corresponding coefficient of
* the input polynomial overflow into the high bits.
*
* Arguments:   - poly *h: pointer to output hint polynomial
*              - const poly *a0: pointer to low part of input polynomial
*              - const poly *a1: pointer to high part of input polynomial
*
* Returns number of 1 bits.
*/
word32 Dilithium::PolyMakeHint(poly *h, const poly *a0, const poly *a1) {
  word32 i, s = 0;

  for(i = 0; i < mN; ++i) {
    h->at(i) = MakeHint(a0->at(i), a1->at(i));
    s += h->at(i);
  }

  return s;
}

/*
* Use hint polynomial to correct the high bits of a polynomial.
*
* Arguments:   - poly *b: pointer to output polynomial with corrected high bits
*              - const poly *a: pointer to input polynomial
*              - const poly *h: pointer to input hint polynomial
*/
void Dilithium::PolyUseHint(poly *b, const poly *a, const poly *h) {
  word32 i;

  for(i = 0; i < mN; ++i)
    b->at(i) = UseHint(a->at(i), h->at(i));

}

/*
* Check infinity norm of polynomial against given bound.
* Assumes input coefficients to be standard representatives.
*
* Arguments:   - const poly *a: pointer to polynomial
*              - word32 B: norm bound
*
* Returns 0 if norm is strictly smaller than B and 1 otherwise.
*/
int Dilithium::PolyChkNorm(const poly *a, word32 B) {
  word32 i;
  word32 t;

  /* It is ok to leak which coefficient violates the bound since
     the probability for each coefficient is independent of secret
     data but we must not leak the sign of the centralized representative. */
  for(i = 0; i < mN; ++i) {
    /* Absolute value of centralized representative */
    t = (mQ-1)/2 - a->at(i);
    t ^= (sword32)t >> 31;
    t = (mQ-1)/2 - t;

    if(t >= B) {
      return 1;
    }
  }

  return 0;
}

/*
* Sample uniformly random coefficients in [0, mQ-1] by
* performing rejection sampling using array of random bytes.
*
* Arguments:   - word32 *a: pointer to output array (allocated)
*              - word32 len: number of coefficients to be sampled
*              - const byte *buf: array of random bytes
*              - word32 bufLen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
*/
word32 Dilithium::RejUniform(word32 *a, word32 len, const byte *buf, word32 bufLen)
{
  word32 ctr, pos;
  word32 t;

  ctr = pos = 0;
  while(ctr < len && pos + 3 <= bufLen) {
    t  = buf[pos++];
    t |= (word32)buf[pos++] << 8;
    t |= (word32)buf[pos++] << 16;
    t &= 0x7FFFFF;

    if(t < mQ)
      a[ctr++] = t;
  }

  return ctr;
}

/*
* Sample polynomial with uniformly random coefficients
* in [0,mQ-1] by performing rejection sampling using the
* output stream of SHAKE256(seed|nonce) or AES256CTR(seed,nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const byte seed[]: byte array with seed of length SEEDBYTES
*              - word16 nonce: 2-byte nonce
*/
#define POLY_UNIFORM_NBLOCKS ((768+168-1)/168)
void Dilithium::PolyUniform(poly *a, const byte *seed, word16 nonce)
{
  word32 i, ctr, off;
  word32 bufLen = POLY_UNIFORM_NBLOCKS*SHAKE128_RATE;
  byte buf[POLY_UNIFORM_NBLOCKS*SHAKE128_RATE + 2];
  keccakState state;

  Shake128StreamInit(&state, seed, nonce);
  Shake128SqueezeBlocks(buf, POLY_UNIFORM_NBLOCKS, &state);

  ctr = RejUniform(a->data(), mN, buf, bufLen);

  while(ctr < mN) {
    off = bufLen % 3;
    for(i = 0; i < off; ++i)
      buf[i] = buf[bufLen - off + i];

    bufLen = SHAKE128_RATE + off;
    Shake128SqueezeBlocks(buf + off, 1, &state);
    ctr += RejUniform(a->data() + ctr, mN - ctr, buf, bufLen);
  }
}

/*
* Sample uniformly random coefficients in [-mEta, mEta] by
* performing rejection sampling using array of random bytes.
*
* Arguments:   - word32 *a: pointer to output array (allocated)
*              - word32 len: number of coefficients to be sampled
*              - const byte *buf: array of random bytes
*              - word32 bufLen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
*/
word32 Dilithium::RejEta(word32 *a, word32 len, const byte *buf, word32 bufLen)
{
  word32 ctr, pos;
  word32 t0, t1;

  ctr = pos = 0;
  while(ctr < len && pos < bufLen) {
  if (mEta <= 3) {
    t0 = buf[pos] & 0x07;
    t1 = buf[pos++] >> 5;
  }
  else {
    t0 = buf[pos] & 0x0F;
    t1 = buf[pos++] >> 4;
  }
  if(t0 <= 2*mEta)
    a[ctr++] = mQ + mEta - t0;
  if(t1 <= 2*mEta && ctr < len)
    a[ctr++] = mQ + mEta - t1;
  }

  return ctr;
}

/*
* Sample polynomial with uniformly random coefficients
* in [-mEta,mEta] by performing rejection sampling using the
* output stream from SHAKE256(seed|nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const byte *seed: pointer to byte array with seed of length SEEDBYTES
*              - word16 nonce: 2-byte nonce
**************************************************/
#define POLY_UNIFORM_ETA_NBLOCKS ((192 + 168 - 1) \
                                  /168)
void Dilithium::PolyUniformEta(poly *a, const byte *seed, word16 nonce)
{
  word32 ctr;
  word32 bufLen = POLY_UNIFORM_ETA_NBLOCKS*SHAKE128_RATE;
  byte buf[POLY_UNIFORM_ETA_NBLOCKS*SHAKE128_RATE];
  keccakState state;

  Shake128StreamInit(&state, seed, nonce);
  Shake128SqueezeBlocks(buf, POLY_UNIFORM_ETA_NBLOCKS, &state);

  ctr = RejEta(a->data(), mN, buf, bufLen);

  while(ctr < mN) {
    Shake128SqueezeBlocks(buf, 1, &state);
    ctr += RejEta(a->data() + ctr, mN - ctr, buf, SHAKE128_RATE);
  }
}

/*
* Sample uniformly random coefficients
* in [-(mGamma1 - 1), mGamma1 - 1] by performing rejection sampling
* using array of random bytes.
*
* Arguments:   - word32 *a: pointer to output array (allocated)
*              - word32 len: number of coefficients to be sampled
*              - const byte *buf: array of random bytes
*              - word32 buflen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
*/
word32 Dilithium::RejGamma1m1(word32 *a, word32 len, const byte *buf, word32 buflen)
{
  word32 ctr, pos;
  word32 t0, t1;

  ctr = pos = 0;
  while(ctr < len && pos + 5 <= buflen) {
    t0  = buf[pos];
    t0 |= (word32)buf[pos + 1] << 8;
    t0 |= (word32)buf[pos + 2] << 16;
    t0 &= 0xFFFFF;

    t1  = buf[pos + 2] >> 4;
    t1 |= (word32)buf[pos + 3] << 4;
    t1 |= (word32)buf[pos + 4] << 12;

    pos += 5;

    if(t0 <= 2*mGamma1 - 2)
      a[ctr++] = mQ + mGamma1 - 1 - t0;
    if(t1 <= 2*mGamma1 - 2 && ctr < len)
      a[ctr++] = mQ + mGamma1 - 1 - t1;
  }

  return ctr;
}

/*
* Sample polynomial with uniformly random coefficients
* in [-(mGamma1 - 1), mGamma1 - 1] by performing rejection
* sampling on output stream of SHAKE256(seed|nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const byte *seed: pointer to byte array with seed of length mCrhBytes
*              - word16 nonce: 16-bit nonce
*/
#define POLY_UNIFORM_GAMMA1M1_NBLOCKS ((640 + 136 - 1) \
                                       /136)
void Dilithium::PolyUniformGamma1m1(poly *a, const byte *seed, word16 nonce)
{
  word32 i, ctr, off;
  word32 bufLen = POLY_UNIFORM_GAMMA1M1_NBLOCKS*SHAKE256_RATE;
  byte buf[POLY_UNIFORM_GAMMA1M1_NBLOCKS*SHAKE256_RATE + 4];
  keccakState state;

  Shake256StreamInit(&state, seed, nonce);
  Shake256SqueezeBlocks(buf, POLY_UNIFORM_GAMMA1M1_NBLOCKS, &state);

  ctr = RejGamma1m1(a->data(), mN, buf, bufLen);

  while(ctr < mN) {
    off = bufLen % 5;
    for(i = 0; i < off; ++i)
      buf[i] = buf[bufLen - off + i];

    bufLen = SHAKE256_RATE + off;
    Shake256SqueezeBlocks(buf + off, 1, &state);
    ctr += RejGamma1m1(a->data() + ctr, mN - ctr, buf, bufLen);
  }
}

/*
* Bit-pack polynomial with coefficients in [-mEta,mEta].
* Input coefficients are assumed to lie in [mQ-mEta,mQ+mEta].
*
* Arguments:   - byte *r: pointer to output byte array with at least mPolyEtaPackedBytes bytes
*              - const poly *a: pointer to input polynomial
*/
void Dilithium::PolyEtaPack(byte *r, const poly *a) {
  word32 i;
  byte t[8];

  if (2*mEta <= 7) {
    for(i = 0; i < mN/8; ++i) {
      t[0] = mQ + mEta - a->at(8*i+0);
      t[1] = mQ + mEta - a->at(8*i+1);
      t[2] = mQ + mEta - a->at(8*i+2);
      t[3] = mQ + mEta - a->at(8*i+3);
      t[4] = mQ + mEta - a->at(8*i+4);
      t[5] = mQ + mEta - a->at(8*i+5);
      t[6] = mQ + mEta - a->at(8*i+6);
      t[7] = mQ + mEta - a->at(8*i+7);

      r[3*i+0]  = (t[0] >> 0) | (t[1] << 3) | (t[2] << 6);
      r[3*i+1]  = (t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7);
      r[3*i+2]  = (t[5] >> 1) | (t[6] << 2) | (t[7] << 5);
    }
  }
    
  else {
    for(i = 0; i < mN/2; ++i) {
      t[0] = mQ + mEta - a->at(2*i+0);
      t[1] = mQ + mEta - a->at(2*i+1);
      r[i] = t[0] | (t[1] << 4);
    }
  }

}

/*
* Unpack polynomial with coefficients in [-mEta,mEta].
* Output coefficients lie in [mQ-mEta,mQ+mEta].
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const byte *a: byte array with bit-packed polynomial
*/
void Dilithium::PolyEtaUnpack(poly *r, const byte *a) {
  word32 i;

  if (mEta <= 3) {
    for(i = 0; i < mN/8; ++i) {
      r->at(8*i+0) = a[3*i+0] & 0x07;
      r->at(8*i+1) = (a[3*i+0] >> 3) & 0x07;
      r->at(8*i+2) = ((a[3*i+0] >> 6) | (a[3*i+1] << 2)) & 0x07;
      r->at(8*i+3) = (a[3*i+1] >> 1) & 0x07;
      r->at(8*i+4) = (a[3*i+1] >> 4) & 0x07;
      r->at(8*i+5) = ((a[3*i+1] >> 7) | (a[3*i+2] << 1)) & 0x07;
      r->at(8*i+6) = (a[3*i+2] >> 2) & 0x07;
      r->at(8*i+7) = (a[3*i+2] >> 5) & 0x07;

      r->at(8*i+0) = mQ + mEta - r->at(8*i+0);
      r->at(8*i+1) = mQ + mEta - r->at(8*i+1);
      r->at(8*i+2) = mQ + mEta - r->at(8*i+2);
      r->at(8*i+3) = mQ + mEta - r->at(8*i+3);
      r->at(8*i+4) = mQ + mEta - r->at(8*i+4);
      r->at(8*i+5) = mQ + mEta - r->at(8*i+5);
      r->at(8*i+6) = mQ + mEta - r->at(8*i+6);
      r->at(8*i+7) = mQ + mEta - r->at(8*i+7);
    }
  } else {
    for(i = 0; i < mN/2; ++i) {
      r->at(2*i+0) = a[i] & 0x0F;
      r->at(2*i+1) = a[i] >> 4;
      r->at(2*i+0) = mQ + mEta - r->at(2*i+0);
      r->at(2*i+1) = mQ + mEta - r->at(2*i+1);
    }
  }
}

/*
* Bit-pack polynomial t1 with coefficients fitting in 9 bits.
* Input coefficients are assumed to be standard representatives.
*
* Arguments:   - byte *r: pointer to output byte array with at least mPolyt1PackedBytes bytes
*              - const poly *a: pointer to input polynomial
*/
void Dilithium::Polyt1Pack(byte *r, const poly *a) {
  word32 i;

  for(i = 0; i < mN/8; ++i) {
    r[9*i+0] = (a->at(8*i+0) >> 0);
    r[9*i+1] = (a->at(8*i+0) >> 8) | (a->at(8*i+1) << 1);
    r[9*i+2] = (a->at(8*i+1) >> 7) | (a->at(8*i+2) << 2);
    r[9*i+3] = (a->at(8*i+2) >> 6) | (a->at(8*i+3) << 3);
    r[9*i+4] = (a->at(8*i+3) >> 5) | (a->at(8*i+4) << 4);
    r[9*i+5] = (a->at(8*i+4) >> 4) | (a->at(8*i+5) << 5);
    r[9*i+6] = (a->at(8*i+5) >> 3) | (a->at(8*i+6) << 6);
    r[9*i+7] = (a->at(8*i+6) >> 2) | (a->at(8*i+7) << 7);
    r[9*i+8] = (a->at(8*i+7) >> 1);
  }
}

/*
* Unpack polynomial t1 with 9-bit coefficients.
* Output coefficients are standard representatives.
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const byte *a: byte array with bit-packed polynomial
*/
void Dilithium::Polyt1Unpack(poly *r, const byte *a) {
  word32 i;

  for(i = 0; i < mN/8; ++i) {
    r->at(8*i+0) = ((a[9*i+0] >> 0) | ((word32)a[9*i+1] << 8)) & 0x1FF;
    r->at(8*i+1) = ((a[9*i+1] >> 1) | ((word32)a[9*i+2] << 7)) & 0x1FF;
    r->at(8*i+2) = ((a[9*i+2] >> 2) | ((word32)a[9*i+3] << 6)) & 0x1FF;
    r->at(8*i+3) = ((a[9*i+3] >> 3) | ((word32)a[9*i+4] << 5)) & 0x1FF;
    r->at(8*i+4) = ((a[9*i+4] >> 4) | ((word32)a[9*i+5] << 4)) & 0x1FF;
    r->at(8*i+5) = ((a[9*i+5] >> 5) | ((word32)a[9*i+6] << 3)) & 0x1FF;
    r->at(8*i+6) = ((a[9*i+6] >> 6) | ((word32)a[9*i+7] << 2)) & 0x1FF;
    r->at(8*i+7) = ((a[9*i+7] >> 7) | ((word32)a[9*i+8] << 1)) & 0x1FF;
  }
}

/*
* Bit-pack polynomial t0 with coefficients in ]-2^{D-1}, 2^{D-1}].
* Input coefficients are assumed to lie in ]mQ-2^{D-1}, mQ+2^{D-1}].
*
* Arguments:   - byte *r: pointer to output byte array with at least mPolyt0PackedBytes bytes
*              - const poly *a: pointer to input polynomial
*/
void Dilithium::Polyt0Pack(byte *r, const poly *a) {
  word32 i;
  word32 t[4];

  for(i = 0; i < mN/4; ++i) {
    t[0] = mQ + (1U << (mD-1)) - a->at(4*i+0);
    t[1] = mQ + (1U << (mD-1)) - a->at(4*i+1);
    t[2] = mQ + (1U << (mD-1)) - a->at(4*i+2);
    t[3] = mQ + (1U << (mD-1)) - a->at(4*i+3);

    r[7*i+0]  =  t[0];
    r[7*i+1]  =  t[0] >> 8;
    r[7*i+1] |=  t[1] << 6;
    r[7*i+2]  =  t[1] >> 2;
    r[7*i+3]  =  t[1] >> 10;
    r[7*i+3] |=  t[2] << 4;
    r[7*i+4]  =  t[2] >> 4;
    r[7*i+5]  =  t[2] >> 12;
    r[7*i+5] |=  t[3] << 2;
    r[7*i+6]  =  t[3] >> 6;
  }
}

/*
* Unpack polynomial t0 with coefficients in ]-2^{D-1}, 2^{D-1}].
* Output coefficients lie in ]mQ-2^{D-1},mQ+2^{D-1}].
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const byte *a: byte array with bit-packed polynomial
*/
void Dilithium::Polyt0Unpack(poly *r, const byte *a) {
  word32 i;

  for(i = 0; i < mN/4; ++i) {
    r->at(4*i+0)  = a[7*i+0];
    r->at(4*i+0) |= (word32)a[7*i+1] << 8;
    r->at(4*i+0) &= 0x3FFF;

    r->at(4*i+1)  = a[7*i+1] >> 6;
    r->at(4*i+1) |= (word32)a[7*i+2] << 2;
    r->at(4*i+1) |= (word32)a[7*i+3] << 10;
    r->at(4*i+1) &= 0x3FFF;

    r->at(4*i+2)  = a[7*i+3] >> 4;
    r->at(4*i+2) |= (word32)a[7*i+4] << 4;
    r->at(4*i+2) |= (word32)a[7*i+5] << 12;
    r->at(4*i+2) &= 0x3FFF;

    r->at(4*i+3)  = a[7*i+5] >> 2;
    r->at(4*i+3) |= (word32)a[7*i+6] << 6;

    r->at(4*i+0) = mQ + (1U << (mD-1)) - r->at(4*i+0);
    r->at(4*i+1) = mQ + (1U << (mD-1)) - r->at(4*i+1);
    r->at(4*i+2) = mQ + (1U << (mD-1)) - r->at(4*i+2);
    r->at(4*i+3) = mQ + (1U << (mD-1)) - r->at(4*i+3);
  }
}

/*
* Bit-pack polynomial z with coefficients
* in [-(mGamma1 - 1), mGamma1 - 1].
* Input coefficients are assumed to be standard representatives.
*
* Arguments:   - byte *r: pointer to output byte array with at least mPolyzPackedBytes bytes
*              - const poly *a: pointer to input polynomial
*/
void Dilithium::PolyzPack(byte *r, const poly *a) {
  word32 i;
  word32 t[2];

  for(i = 0; i < mN/2; ++i) {
    /* Map to {0,...,2*mGamma1 - 2} */
    t[0] = mGamma1 - 1 - a->at(2*i+0);
    t[0] += ((sword32)t[0] >> 31) & mQ;
    t[1] = mGamma1 - 1 - a->at(2*i+1);
    t[1] += ((sword32)t[1] >> 31) & mQ;

    r[5*i+0]  = t[0];
    r[5*i+1]  = t[0] >> 8;
    r[5*i+2]  = t[0] >> 16;
    r[5*i+2] |= t[1] << 4;
    r[5*i+3]  = t[1] >> 4;
    r[5*i+4]  = t[1] >> 12;
  }
}

/*
* Unpack polynomial z with coefficients
* in [-(mGamma1 - 1), mGamma1 - 1].
* Output coefficients are standard representatives.
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const byte *a: byte array with bit-packed polynomial
*/
void Dilithium::PolyzUnpack(poly *r, const byte *a) {
  word32 i;

  for(i = 0; i < mN/2; ++i) {
    r->at(2*i+0)  = a[5*i+0];
    r->at(2*i+0) |= (word32)a[5*i+1] << 8;
    r->at(2*i+0) |= (word32)a[5*i+2] << 16;
    r->at(2*i+0) &= 0xFFFFF;

    r->at(2*i+1)  = a[5*i+2] >> 4;
    r->at(2*i+1) |= (word32)a[5*i+3] << 4;
    r->at(2*i+1) |= (word32)a[5*i+4] << 12;

    r->at(2*i+0) = mGamma1 - 1 - r->at(2*i+0);
    r->at(2*i+0) += ((sword32)r->at(2*i+0) >> 31) & mQ;
    r->at(2*i+1) = mGamma1 - 1 - r->at(2*i+1);
    r->at(2*i+1) += ((sword32)r->at(2*i+1) >> 31) & mQ;
  }
}

/*
* Bit-pack polynomial w1 with coefficients in [0, 15].
* Input coefficients are assumed to be standard representatives.
*
* Arguments:   - byte *r: pointer to output byte array with at least
*                            mPolyw1PackedBytes bytes
*              - const poly *a: pointer to input polynomial
*/
void Dilithium::Polyw1Pack(byte *r, const poly *a) {
  word32 i;

  for(i = 0; i < mN/2; ++i)
    r[i] = a->at(2*i+0) | (a->at(2*i+1) << 4);
}

/*
* Implementation of ExpandA. Generates matrix A with uniformly
* random coefficients a_{i,j} by performing rejection
* sampling on the output stream of SHAKE128(rho|j|i).
*
* Arguments:   - std::vector<polyvec> *mat: pointer to output matrix
*              - const byte *rho: pointer to byte array containing seed rho
**************************************************/
void Dilithium::ExpandMat(std::vector<polyvec> *mat, const byte *rho) {
  word32 i, j;
  for(i = 0; i < mK; ++i)
    for(j = 0; j < mL; ++j)
      PolyUniform(&mat->at(i).at(j), rho, (i << 8) + j);
}



/* Reduce coefficients of polynomials in vector of length mL or mK
*  to standard representatives.
*
* Arguments:   - polyvec *v: pointer to input/output vector
*              - byte vectorLen: length of the input/output vector
*/
void Dilithium::PolyvecFreeze(polyvec *v, byte vectorLen) {
  word32 i;

  for(i = 0; i < vectorLen; ++i)
    PolyFreeze(&v->at(i));
}

/* Add vectors of polynomials of length mL or mK.
*  No modular reduction is performed.
*
* Arguments:   - polyvec *w: pointer to output vector
*              - const polyvec *u: pointer to first summand
*              - const polyvec *v: pointer to second summand
*              - byte vectorLen: length of the vectors
*/
void Dilithium::PolyvecAdd(polyvec *w, const polyvec *u, const polyvec *v, byte vectorLen) {
  word32 i;

  for(i = 0; i < vectorLen; ++i)
    PolyAdd(&w->at(i), &u->at(i), &v->at(i));
}

/*
* Forward NTT of all polynomials in vector of length mL or mK. Output
* coefficients can be up to 16*mQ larger than input coefficients.
*
* Arguments:   - polyvec *v: pointer to input/output vector
*              - byte vectorLen: lenght of the input/output vector
*/
void Dilithium::PolyvecNtt(polyvec *v, byte vectorLen) {
  word32 i;

  for(i = 0; i < vectorLen; ++i)
    PolyNtt(&v->at(i));
}

/*
* Pointwise multiply vectors of polynomials of length mL, multiply
* resulting vector by 2^{-32} and add (accumulate) polynomials
* in it. Input/output vectors are in NTT domain representation.
* Input coefficients are assumed to be less than 22*mQ. Output
* coeffcient are less than 2*mL*mQ.
*
* Arguments:   - poly *w: output polynomial
*              - const polyvec *u: pointer to first input vector
*              - const polyvec *v: pointer to second input vector
*              - byte vectorLen: lenght of the vectors
*/
void Dilithium::PolyvecPointwiseAccMontgomery(poly *w, const polyvec *u, const polyvec *v, byte vectorLen)
{
  word32 i;
  poly t;

  PolyPointwiseMontgomery(w, &u->at(0), &v->at(0));

  for(i = 1; i < vectorLen; ++i) {
    PolyPointwiseMontgomery(&t, &u->at(i), &v->at(i));
    PolyAdd(w, w, &t);
  }
}

/*
* Check infinity norm of polynomials in vector of length mL or mK.
* Assumes input coefficients to be standard representatives.
*
* Arguments:   - const polyvec *v: pointer to vector
*              - word32 B: norm bound
*              - byte vectorLen: lenght of vector
*
* Returns 0 if norm of all polynomials is strictly smaller than B and 1
* otherwise.
*/
int Dilithium::PolyvecChkNorm(const polyvec *v, word32 bound, byte vectorLen)  {
  word32 i;

  for(i = 0; i < vectorLen; ++i)
    if(PolyChkNorm(&v->at(i), bound))
      return 1;

  return 0;
}



/*
* Reduce coefficients of polynomials in vector of length mK
* to representatives in [0,2*mQ[.
*
* Arguments:   - polyvec *v: pointer to input/output vector
*              - byte vectorLen: length of vector
*/
void Dilithium::PolyvecReduce(polyvec *v, byte vectorLen) {
  word32 i;

  for(i = 0; i < vectorLen; ++i)
    PolyReduce(&v->at(i));
}

/*
* For all coefficients of polynomials in vector of length mK
* subtract mQ if coefficient is bigger than mQ.
*
* Arguments:   - polyvec *v: pointer to input/output vector
               - byte vectorLen: length of vector
*/
void Dilithium::PolyvecCsubQ(polyvec *v, byte vectorLen) {
  word32 i;

  for(i = 0; i < vectorLen; ++i)
    PolyCsubQ(&v->at(i));
}



/*
* Subtract vectors of polynomials of length mK.
* Assumes coefficients of polynomials in second input vector
* to be less than 2*mQ. No modular reduction is performed.
*
* Arguments:   - polyvec *w: pointer to output vector
*              - const polyvec *u: pointer to first input vector
*              - const polyvec *v: pointer to second input vector to be
*                                   subtracted from first input vector
*              - byte vectorLen: lenght of vectors
*/
void Dilithium::PolyvecSub(polyvec *w, const polyvec *u, const polyvec *v, byte vectorLen) {
  word32 i;

  for(i = 0; i < vectorLen; ++i)
    PolySub(&w->at(i), &u->at(i), &v->at(i));
}

/*
* Multiply vector of polynomials of Length mK by 2^mD without modular
* reduction. Assumes input coefficients to be less than 2^{32-mD}.
*
* Arguments:   - polyvec *v: pointer to input/output vector
*              - byte vectorLen: length of input/output vector
*/
void Dilithium::PolyvecShiftL(polyvec *v, byte vectorLen) {
  word32 i;

  for(i = 0; i < vectorLen; ++i)
    PolyShiftL(&v->at(i));
}


/*
* Inverse NTT and multiplication by 2^{32} of polynomials
* in vector of length mK. Input coefficients need to be less
* than 2*mQ.
*
* Arguments:   - polyvec *v: pointer to input/output vector
*              - byte vectorLen: length of input/output vector
*/
void Dilithium::PolyvecInvNttToMont(polyvec *v, byte vectorLen) {
  word32 i;

  for(i = 0; i < vectorLen; ++i)
    PolyInvNttToMont(&v->at(i));
}

/*
* For all coefficients a of polynomials in vector of length mK,
* compute a0, a1 such that a mod mQ = a1*2^mD + a0
* with -2^{mD-1} < a0 <= 2^{mD-1}. Assumes coefficients to be
* standard representatives.
*
* Arguments:   - polyvec *v1: pointer to output vector of polynomials with
*                              coefficients a1
*              - polyvec *v0: pointer to output vector of polynomials with
*                              coefficients mQ + a0
*              - const polyvec *v: pointer to input vector
*              - byte vectorLen: length of vectors
*/
void Dilithium::PolyvecPower2Round(polyvec *v1, polyvec *v0, const polyvec *v, byte vectorLen) {
  word32 i;

  for(i = 0; i < vectorLen; ++i)
    PolyPower2Round(&v1->at(i), &v0->at(i), &v->at(i));
}

/*
* For all coefficients a of polynomials in vector of length K,
* compute high and low bits a0, a1 such a mod mQ = a1*mAlpha + a0
* with -mAlpha/2 < a0 <= mAlpha/2 except a1 = (mQ-1)/mAlpha where we
* set a1 = 0 and -mAlpha/2 <= a0 = a mod mQ - mQ < 0.
* Assumes coefficients to be standard representatives.
*
* Arguments:   - polyvec *v1: pointer to output vector of polynomials with
*                              coefficients a1
*              - polyvec *v0: pointer to output vector of polynomials with
*                              coefficients mQ + a0
*              - const polyvec *v: pointer to input vector
*              - byte vectorLen: length of vectors
*/
void Dilithium::PolyvecDecompose(polyvec *v1, polyvec *v0, const polyvec *v, byte vectorLen) {
  word32 i;

  for(i = 0; i < vectorLen; ++i)
    PolyDecompose(&v1->at(i), &v0->at(i), &v->at(i));
}

/*
* Compute hint vector.
*
* Arguments:   - polyvec *h: pointer to output vector
*              - const polyvec *v0: pointer to low part of input vector
*              - const polyvec *v1: pointer to high part of input vector
*              - byte vectorLen: length of vectors
*
* Returns number of 1 bits.
*/
word32 Dilithium::PolyvecMakeHint(polyvec *h, const polyvec *v0, const polyvec *v1, byte vectorLen)
{
  word32 i, s = 0;

  for(i = 0; i < vectorLen; ++i)
    s += PolyMakeHint(&h->at(i), &v0->at(i), &v1->at(i));

  return s;
}

/*
* Use hint vector to correct the high bits of input vector.
*
* Arguments:   - polyvec *w: pointer to output vector of polynomials with
*                             corrected high bits
*              - const polyvec *u: pointer to input vector
*              - const polyvec *h: pointer to input hint vector
*              - byte vectorLen: length of vectors
*/
void Dilithium::PolyvecUseHint(polyvec *w, const polyvec *u, const polyvec *h, byte vectorLen) {
  word32 i;

  for(i = 0; i < vectorLen; ++i)
    PolyUseHint(&w->at(i), &u->at(i), &h->at(i));
}


NAMESPACE_END