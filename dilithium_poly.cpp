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
* representative in [-6283009,6283007].
*
* Arguments:   - poly *a: pointer to input/output polynomial
*/
void Dilithium::PolyReduce(poly *a) {
  word32 i;

  for(i = 0; i < mN; ++i)
    a->at(i) = Reduce32(a->at(i));
}


/*
*For all coefficients of in/out polynomial add mQ if
*              coefficient is negative.
*
* Arguments:   - poly *a: pointer to input/output polynomial
*/
void Dilithium::PolyCAddQ(poly *a) {
  word32 i;

  for(i = 0; i < mN; ++i)
    a->at(i) = CAddQ(a->at(i));
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
* Subtract polynomials. No modular reduction is
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
    c->at(i) = a->at(i) - b->at(i);

}

/*
* Multiply polynomial by 2^D without modular reduction. Assumes
* input coefficients to be less than 2^{31-mD} in absolute value.
*
* Arguments:   - poly *a: pointer to input/output polynomial
*/
void Dilithium::PolyShiftL(poly *a) {
  word32 i;

  for(i = 0; i < mN; ++i)
    a->at(i) <<= mD;

}

/*
* Inplace forward NTT. Output coefficients can grow by
*              8*mQ in absolute value.
*
* Arguments:   - poly *a: pointer to input/output polynomial
*/
void Dilithium::PolyNtt(poly *a) {

  Ntt(a->data());

}

void Dilithium::PolyvecInvNttToMont(polyvec *v, const byte vectorLen) {
  word32 i;

  for(i = 0; i < vectorLen; ++i)
    PolyInvNttToMont(&v->at(i));
}


void Dilithium::PolyvecPointwisePolyMontgomery(polyvec *r, const poly *a, const polyvec *v, const byte vectorLen) {
  word32 i;

  for(i = 0; i < vectorLen; ++i)
    PolyPointwiseMontgomery(&r->at(i), a, &v->at(i));
}

/*
* Inplace inverse NTT and multiplication by 2^{32}.
* Input coefficients need to be less than mQ in absolute value.
* Output coefficients are less than mQ in absolute value.
*
* Arguments:   - poly *a: pointer to input/output polynomial
*/
void Dilithium::PolyInvNttToMont(poly *a) {

  InvNttToMont(a->data());

}

/*
* Pointwise multiplication of polynomials in NTT domain
* representation and multiplication of resulting polynomial
* by 2^{-32}. 
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
*/
void Dilithium::PolyPointwiseMontgomery(poly *c, const poly *a, const poly *b) {
  word32 i;

  for(i = 0; i < mN; ++i)
    c->at(i) = MontgomeryReduce((sword64)a->at(i) * b->at(i));

}

/*
* For all coefficients c of the input polynomial,
* compute c0, c1 such that c mod mQ = c1*2^D + c0
* with -2^{D-1} < c0 <= 2^{D-1}. Assumes coefficients to be
* standard representatives.
*
* Arguments:   - poly *a1: pointer to output polynomial with coefficients c1
*              - poly *a0: pointer to output polynomial with coefficients c0
*              - const poly *a: pointer to input polynomial
*/
void Dilithium::PolyPower2Round(poly *a1, poly *a0, const poly *a) {
  word32 i;

  for(i = 0; i < mN; ++i)
    a1->at(i) = Power2Round(&a0->at(i), a->at(i));

}

/*
* For all coefficients c of the input polynomial,
* compute high and low bits c0, c1 such c mod mQ = c1*mAlpha + c0
* with -mAlpha/2 < c0 <= mAlpha/2 except c1 = (mQ-1)/mAlpha where we
* set c1 = 0 and -mAlpha/2 <= c0 = c mod mQ - mQ < 0.
* Assumes coefficients to be standard representatives.
*
* Arguments:   - poly *a1: pointer to output polynomial with coefficients c1
*              - poly *a0: pointer to output polynomial with coefficients c0
*              - const poly *a: pointer to input polynomial
*/
void Dilithium::PolyDecompose(poly *a1, poly *a0, const poly *a) {
  word32 i;

  for(i = 0; i < mN; ++i)
    a1->at(i) = Decompose(&a0->at(i), a->at(i));

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
* Assumes input coefficients were reduced by Reduce32.
*
* Arguments:   - const poly *a: pointer to polynomial
*              - sword32 B: norm bound
*
* Returns 0 if norm is strictly smaller than B <= (Q-1)/8 and 1 otherwise.
*/
int Dilithium::PolyChkNorm(const poly *a, sword32 B) {
  word32 i;
  sword32 t;
  if(B > (mQ-1)/8) {
    return 1;
  }

  /* It is ok to leak which coefficient violates the bound since
     the probability for each coefficient is independent of secret
     data but we must not leak the sign of the centralized representative. */
  for(i = 0; i < mN; ++i) {
    /* Absolute value of centralized representative */
    t = a->at(i) >> 31;
    t = a->at(i) - (t & 2*a->at(i));

    if(t >= B) {
      return 1;
    }
  }
  return 0;
}

/*
* Sample uniformly random coefficients in [0, mQ-1] by
* performing rejection sampling on array of random bytes.
*
* Arguments:   - sword32 *a: pointer to output array (allocated)
*              - word32 len: number of coefficients to be sampled
*              - const byte *buf: array of random bytes
*              - word32 bufLen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
*/
word32 Dilithium::RejUniform(sword32 *a, word32 len, const byte *buf, word32 bufLen)
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
* in [0,mQ-1] by performing rejection sampling on the
* output stream of SHAKE256(seed|nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const byte seed[]: byte array with seed of length SEEDBYTES
*              - word16 nonce: 2-byte nonce
*/
void Dilithium::PolyUniform(poly *a, const byte *seed, word16 nonce)
{
  word32 i, ctr, off;
  word32 mPolyUniformNBlocks = (768 + SHAKE128_RATE - 1)/SHAKE128_RATE; 
  word32 bufLen = mPolyUniformNBlocks*SHAKE128_RATE;
  byte buf[mPolyUniformNBlocks*SHAKE128_RATE + 2];

  Shake128StreamInit(seed, nonce);
  Shake128SqueezeBlocks(buf, mPolyUniformNBlocks);

  ctr = RejUniform(a->data(), mN, buf, bufLen);

  while(ctr < mN) {
    off = bufLen % 3;
    for(i = 0; i < off; ++i)
      buf[i] = buf[bufLen - off + i];

    
    Shake128SqueezeBlocks(buf + off, 1);
    bufLen = SHAKE128_RATE + off;
    ctr += RejUniform(a->data() + ctr, mN - ctr, buf, bufLen);
  }
}

/*
* Sample uniformly random coefficients in [-mEta, mEta] by
* performing rejection sampling on array of random bytes.
*
* Arguments:   - sword32 *a: pointer to output array (allocated)
*              - word32 len: number of coefficients to be sampled
*              - const byte *buf: array of random bytes
*              - word32 bufLen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
*/
word32 Dilithium::RejEta(sword32 *a, word32 len, const byte *buf, word32 bufLen)
{
  word32 ctr, pos;
  word32 t0, t1;

  ctr = pos = 0;
  while(ctr < len && pos < bufLen) {
    t0 = buf[pos] & 0x0F;
    t1 = buf[pos++] >> 4;
    if (mEta == 2) {
      if(t0 < 15) {
        t0 = t0 - (205*t0 >> 10)*5;
        a[ctr++] = 2 - t0;
      }
      if(t1 < 15 && ctr < len) {
        t1 = t1 - (205*t1 >> 10)*5;
        a[ctr++] = 2 - t1;
      }
    } else { //mEta == 4

      if(t0 < 9)
        a[ctr++] = 4 - t0;
      if(t1 < 9 && ctr < len)
        a[ctr++] = 4 - t1;
    }
  }
  return ctr;
}

/*
* Sample polynomial with uniformly random coefficients
* in [-mEta,mEta] by performing rejection sampling on the
* output stream from SHAKE256(seed|nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const byte *seed: pointer to byte array with seed of length mCrhBytes
*              - word16 nonce: 2-byte nonce
**************************************************/

void Dilithium::PolyUniformEta(poly *a, const byte *seed, word16 nonce)
{
  word32 ctr;
  word32 polyUniformEtaNBlocks;
  if (mEta == 2) {
    polyUniformEtaNBlocks = (136 + SHAKE256_RATE - 1)/SHAKE256_RATE;
  } else { //mEta == 4
    polyUniformEtaNBlocks = (227 + SHAKE256_RATE - 1)/SHAKE256_RATE;

  };

  word32 bufLen = polyUniformEtaNBlocks*SHAKE256_RATE;
  std::vector<byte> bufVector(bufLen);
  std::vector<byte> *buf = &bufVector;

  Shake256StreamInit(seed, nonce);
  Shake256SqueezeBlocks(buf->data(), polyUniformEtaNBlocks);

  ctr = RejEta(a->data(), mN, buf->data(), bufLen);

  while(ctr < mN) {
    Shake256SqueezeBlocks(buf->data(), 1);
    ctr += RejEta(a->data() + ctr, mN - ctr, buf->data(), SHAKE256_RATE);
  }
}


/*  
* Sample polynomial with uniformly random coefficients
* in [-(mGamma1 - 1), mGamma1] by performing rejection
* sampling on output stream of SHAKE256(seed|nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const byte *seed: pointer to byte array with seed of length mCrhBytes
*              - word16 nonce: 16-bit nonce
*/

void Dilithium::PolyUniformGamma1(poly *a, const byte *seed, word16 nonce)
{
  word16 polyUniformGamma1NBlocks = (mPolyzPackedBytes + SHAKE256_RATE - 1)/SHAKE256_RATE;
  std::vector<byte> bufVector(polyUniformGamma1NBlocks*SHAKE256_RATE);
  std::vector<byte> *buf = &bufVector;

  Shake256StreamInit(seed, nonce);
  Shake256SqueezeBlocks(buf->data(), polyUniformGamma1NBlocks);
  PolyzUnpack(a, buf->data());
}

/*
* Bit-pack polynomial with coefficients in [-mEta,mEta].
* Arguments:   - byte *r: pointer to output byte array with at least mPolyEtaPackedBytes bytes
*              - const poly *a: pointer to input polynomial
*/
void Dilithium::PolyEtaPack(byte *r, const poly *a) {
  word32 i;
  byte t[8];

  if (mEta == 2) {
    for(i = 0; i < mN/8; ++i) {
      t[0] = mEta - a->at(8*i+0);
      t[1] = mEta - a->at(8*i+1);
      t[2] = mEta - a->at(8*i+2);
      t[3] = mEta - a->at(8*i+3);
      t[4] = mEta - a->at(8*i+4);
      t[5] = mEta - a->at(8*i+5);
      t[6] = mEta - a->at(8*i+6);
      t[7] = mEta - a->at(8*i+7);

      r[3*i+0]  = (t[0] >> 0) | (t[1] << 3) | (t[2] << 6);
      r[3*i+1]  = (t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7);
      r[3*i+2]  = (t[5] >> 1) | (t[6] << 2) | (t[7] << 5);
    }
  }
    
  else { //mEta == 4
    for(i = 0; i < mN/2; ++i) {
      t[0] = mEta - a->at(2*i+0);
      t[1] = mEta - a->at(2*i+1);
      r[i] = t[0] | (t[1] << 4);
    }
  }

}

/*
* Unpack polynomial with coefficients in [-mEta,mEta].
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const byte *a: byte array with bit-packed polynomial
*/
void Dilithium::PolyEtaUnpack(poly *r, const byte *a) {
  word32 i;

  if (mEta == 2) {
    for(i = 0; i < mN/8; ++i) {
      r->at(8*i+0) = (a[3*i+0] >> 0) & 7;
      r->at(8*i+1) = (a[3*i+0] >> 3) & 7;
      r->at(8*i+2) = ((a[3*i+0] >> 6) | (a[3*i+1] << 2)) & 7;
      r->at(8*i+3) = (a[3*i+1] >> 1) & 7;
      r->at(8*i+4) = (a[3*i+1] >> 4) & 7;
      r->at(8*i+5) = ((a[3*i+1] >> 7) | (a[3*i+2] << 1)) & 7;
      r->at(8*i+6) = (a[3*i+2] >> 2) & 7;
      r->at(8*i+7) = (a[3*i+2] >> 5) & 7;

      r->at(8*i+0) = mEta - r->at(8*i+0);
      r->at(8*i+1) = mEta - r->at(8*i+1);
      r->at(8*i+2) = mEta - r->at(8*i+2);
      r->at(8*i+3) = mEta - r->at(8*i+3);
      r->at(8*i+4) = mEta - r->at(8*i+4);
      r->at(8*i+5) = mEta - r->at(8*i+5);
      r->at(8*i+6) = mEta - r->at(8*i+6);
      r->at(8*i+7) = mEta - r->at(8*i+7);
    }
  } else { //mEta == 4
    for(i = 0; i < mN/2; ++i) {
      r->at(2*i+0) = a[i] & 0x0F;
      r->at(2*i+1) = a[i] >> 4;
      r->at(2*i+0) = mEta - r->at(2*i+0);
      r->at(2*i+1) = mEta - r->at(2*i+1);
    }
  }
}

/*
* Bit-pack polynomial t1 with coefficients fitting in 10 bits.
* Input coefficients are assumed to be standard representatives.
*
* Arguments:   - byte *r: pointer to output byte array with at least mPolyt1PackedBytes bytes
*              - const poly *a: pointer to input polynomial
*/
void Dilithium::Polyt1Pack(byte *r, const poly *a) {
  word32 i;

  for(i = 0; i < mN/4; ++i) {
    r[5*i+0] = (a->at(4*i+0) >> 0);
    r[5*i+1] = (a->at(4*i+0) >> 8) | (a->at(4*i+1) << 2);
    r[5*i+2] = (a->at(4*i+1) >> 6) | (a->at(4*i+2) << 4);
    r[5*i+3] = (a->at(4*i+2) >> 4) | (a->at(4*i+3) << 6);
    r[5*i+4] = (a->at(4*i+3) >> 2);
  }
}

/*
* Unpack polynomial t1 with 10-bit coefficients.
* Output coefficients are standard representatives.
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const byte *a: byte array with bit-packed polynomial
*/
void Dilithium::Polyt1Unpack(poly *r, const byte *a) {
  word32 i;

  for(i = 0; i < mN/4; ++i) {
    r->at(4*i+0) = ((a[5*i+0] >> 0) | ((word32)a[5*i+1] << 8)) & 0x3FF;
    r->at(4*i+1) = ((a[5*i+1] >> 2) | ((word32)a[5*i+2] << 6)) & 0x3FF;
    r->at(4*i+2) = ((a[5*i+2] >> 4) | ((word32)a[5*i+3] << 4)) & 0x3FF;
    r->at(4*i+3) = ((a[5*i+3] >> 6) | ((word32)a[5*i+4] << 2)) & 0x3FF;
  }
}

/*
* Bit-pack polynomial t0 with coefficients in ]-2^{mD-1}, 2^{mD-1}].
*
* Arguments:   - byte *r: pointer to output byte array with at least mPolyt0PackedBytes bytes
*              - const poly *a: pointer to input polynomial
*/
void Dilithium::Polyt0Pack(byte *r, const poly *a) {
  word32 i;
  word32 t[8];

  for(i = 0; i < mN/8; ++i) {
    t[0] = (1 << (mD-1)) - a->at(8*i+0);
    t[1] = (1 << (mD-1)) - a->at(8*i+1);
    t[2] = (1 << (mD-1)) - a->at(8*i+2);
    t[3] = (1 << (mD-1)) - a->at(8*i+3);
    t[4] = (1 << (mD-1)) - a->at(8*i+4);
    t[5] = (1 << (mD-1)) - a->at(8*i+5);
    t[6] = (1 << (mD-1)) - a->at(8*i+6);
    t[7] = (1 << (mD-1)) - a->at(8*i+7);

    r[13*i+ 0]  =  t[0];
    r[13*i+ 1]  =  t[0] >>  8;
    r[13*i+ 1] |=  t[1] <<  5;
    r[13*i+ 2]  =  t[1] >>  3;
    r[13*i+ 3]  =  t[1] >> 11;
    r[13*i+ 3] |=  t[2] <<  2;
    r[13*i+ 4]  =  t[2] >>  6;
    r[13*i+ 4] |=  t[3] <<  7;
    r[13*i+ 5]  =  t[3] >>  1;
    r[13*i+ 6]  =  t[3] >>  9;
    r[13*i+ 6] |=  t[4] <<  4;
    r[13*i+ 7]  =  t[4] >>  4;
    r[13*i+ 8]  =  t[4] >> 12;
    r[13*i+ 8] |=  t[5] <<  1;
    r[13*i+ 9]  =  t[5] >>  7;
    r[13*i+ 9] |=  t[6] <<  6;
    r[13*i+10]  =  t[6] >>  2;
    r[13*i+11]  =  t[6] >> 10;
    r[13*i+11] |=  t[7] <<  3;
    r[13*i+12]  =  t[7] >>  5;
  }
}

/*
* Unpack polynomial t0 with coefficients in ]-2^{mD-1}, 2^{mD-1}].
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const byte *a: byte array with bit-packed polynomial
*/
void Dilithium::Polyt0Unpack(poly *r, const byte *a) {
  word32 i;

  for(i = 0; i < mN/8; ++i) {
    r->at(8*i+0)  = a[13*i+0];
    r->at(8*i+0) |= (word32)a[13*i+1] << 8;
    r->at(8*i+0) &= 0x1FFF;

    r->at(8*i+1)  = a[13*i+1] >> 5;
    r->at(8*i+1) |= (word32)a[13*i+2] << 3;
    r->at(8*i+1) |= (word32)a[13*i+3] << 11;
    r->at(8*i+1) &= 0x1FFF;

    r->at(8*i+2)  = a[13*i+3] >> 2;
    r->at(8*i+2) |= (word32)a[13*i+4] << 6;
    r->at(8*i+2) &= 0x1FFF;

    r->at(8*i+3)  = a[13*i+4] >> 7;
    r->at(8*i+3) |= (word32)a[13*i+5] << 1;
    r->at(8*i+3) |= (word32)a[13*i+6] << 9;
    r->at(8*i+3) &= 0x1FFF;

    r->at(8*i+4)  = a[13*i+6] >> 4;
    r->at(8*i+4) |= (word32)a[13*i+7] << 4;
    r->at(8*i+4) |= (word32)a[13*i+8] << 12;
    r->at(8*i+4) &= 0x1FFF;

    r->at(8*i+5)  = a[13*i+8] >> 1;
    r->at(8*i+5) |= (word32)a[13*i+9] << 7;
    r->at(8*i+5) &= 0x1FFF;

    r->at(8*i+6)  = a[13*i+9] >> 6;
    r->at(8*i+6) |= (word32)a[13*i+10] << 2;
    r->at(8*i+6) |= (word32)a[13*i+11] << 10;
    r->at(8*i+6) &= 0x1FFF;

    r->at(8*i+7)  = a[13*i+11] >> 3;
    r->at(8*i+7) |= (word32)a[13*i+12] << 5;
    r->at(8*i+7) &= 0x1FFF;

    r->at(8*i+0) = (1 << (mD-1)) - r->at(8*i+0);
    r->at(8*i+1) = (1 << (mD-1)) - r->at(8*i+1);
    r->at(8*i+2) = (1 << (mD-1)) - r->at(8*i+2);
    r->at(8*i+3) = (1 << (mD-1)) - r->at(8*i+3);
    r->at(8*i+4) = (1 << (mD-1)) - r->at(8*i+4);
    r->at(8*i+5) = (1 << (mD-1)) - r->at(8*i+5);
    r->at(8*i+6) = (1 << (mD-1)) - r->at(8*i+6);
    r->at(8*i+7) = (1 << (mD-1)) - r->at(8*i+7);
  }
}

/*
* Bit-pack polynomial z with coefficients
* in [-(mGamma1 - 1), mGamma1].
* Input coefficients are assumed to be standard representatives.
*
* Arguments:   - byte *r: pointer to output byte array with at least mPolyzPackedBytes bytes
*              - const poly *a: pointer to input polynomial
*/
void Dilithium::PolyzPack(byte *r, const poly *a) {
  word32 i;
  word32 t[4];

  if (mGamma1 == (1 << 17)) {
    for(i = 0; i < mN/4; ++i) {
      t[0] = mGamma1 - a->at(4*i+0);
      t[1] = mGamma1 - a->at(4*i+1);
      t[2] = mGamma1 - a->at(4*i+2);
      t[3] = mGamma1 - a->at(4*i+3);

      r[9*i+0]  = t[0];
      r[9*i+1]  = t[0] >> 8;
      r[9*i+2]  = t[0] >> 16;
      r[9*i+2] |= t[1] << 2;
      r[9*i+3]  = t[1] >> 6;
      r[9*i+4]  = t[1] >> 14;
      r[9*i+4] |= t[2] << 4;
      r[9*i+5]  = t[2] >> 4;
      r[9*i+6]  = t[2] >> 12;
      r[9*i+6] |= t[3] << 6;
      r[9*i+7]  = t[3] >> 2;
      r[9*i+8]  = t[3] >> 10;
    }
  }
  
  else if (mGamma1 == (1 << 19)) {
    for(i = 0; i < mN/2; ++i) {
      t[0] = mGamma1 - a->at(2*i+0);
      t[1] = mGamma1 - a->at(2*i+1);

      r[5*i+0]  = t[0];
      r[5*i+1]  = t[0] >> 8;
      r[5*i+2]  = t[0] >> 16;
      r[5*i+2] |= t[1] << 4;
      r[5*i+3]  = t[1] >> 4;
      r[5*i+4]  = t[1] >> 12;
    }

  }
}

/*
* Unpack polynomial z with coefficients
* in [-(mGamma1 - 1), mGamma1].
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const byte *a: byte array with bit-packed polynomial
*/
void Dilithium::PolyzUnpack(poly *r, const byte *a) {
  word32 i;

  if (mGamma1 == (1 << 17)) {
    for(i = 0; i < mN/4; ++i) {
      r->at(4*i+0)  = a[9*i+0];
      r->at(4*i+0) |= (word32)a[9*i+1] << 8;
      r->at(4*i+0) |= (word32)a[9*i+2] << 16;
      r->at(4*i+0) &= 0x3FFFF;

      r->at(4*i+1)  = a[9*i+2] >> 2;
      r->at(4*i+1) |= (word32)a[9*i+3] << 6;
      r->at(4*i+1) |= (word32)a[9*i+4] << 14;
      r->at(4*i+1) &= 0x3FFFF;

      r->at(4*i+2)  = a[9*i+4] >> 4;
      r->at(4*i+2) |= (word32)a[9*i+5] << 4;
      r->at(4*i+2) |= (word32)a[9*i+6] << 12;
      r->at(4*i+2) &= 0x3FFFF;

      r->at(4*i+3)  = a[9*i+6] >> 6;
      r->at(4*i+3) |= (word32)a[9*i+7] << 2;
      r->at(4*i+3) |= (word32)a[9*i+8] << 10;
      r->at(4*i+3) &= 0x3FFFF;

      r->at(4*i+0) = mGamma1 - r->at(4*i+0);
      r->at(4*i+1) = mGamma1 - r->at(4*i+1);
      r->at(4*i+2) = mGamma1 - r->at(4*i+2);
      r->at(4*i+3) = mGamma1 - r->at(4*i+3);
    }
  }
  
  else if (mGamma1 == (1 << 19)) {
    for(i = 0; i < mN/2; ++i) {
      r->at(2*i+0)  = a[5*i+0];
      r->at(2*i+0) |= (word32)a[5*i+1] << 8;
      r->at(2*i+0) |= (word32)a[5*i+2] << 16;
      r->at(2*i+0) &= 0xFFFFF;

      r->at(2*i+1)  = a[5*i+2] >> 4;
      r->at(2*i+1) |= (word32)a[5*i+3] << 4;
      r->at(2*i+1) |= (word32)a[5*i+4] << 12;
      r->at(2*i+0) &= 0xFFFFF;

      r->at(2*i+0) = mGamma1 - r->at(2*i+0);
      r->at(2*i+1) = mGamma1 - r->at(2*i+1);
    }
  }
}

void Dilithium::Polyvecw1Pack(byte *r, const polyvec *w1, const byte vectorLen) {
  word32 i;

  for(i = 0; i < vectorLen; ++i)
    Polyw1Pack(&r[i*mPolyw1PackedBytes], &w1->at(i));
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
  
  if (mGamma2 == (mQ-1)/88) {
    for(i = 0; i < mN/4; ++i) {
      r[3*i+0]  = a->at(4*i+0);
      r[3*i+0] |= a->at(4*i+1) << 6;
      r[3*i+1]  = a->at(4*i+1) >> 2;
      r[3*i+1] |= a->at(4*i+2) << 4;
      r[3*i+2]  = a->at(4*i+2) >> 4;
      r[3*i+2] |= a->at(4*i+3) << 2;
    }
  }
  else if (mGamma2 == (mQ-1)/32) {
     for(i = 0; i < mN/2; ++i)
      r[i] = a->at(2*i+0) | (a->at(2*i+1) << 4);

  }
 
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

void Dilithium::PolyvecUniformEta(polyvec *v, const byte *seed, word16 nonce, const byte vectorLen) {
  word32 i;

  for(i = 0; i < vectorLen; ++i)
    PolyUniformEta(&v->at(i), seed, nonce++);
}

void Dilithium::PolyvecUniformGamma1(polyvec *v, const byte *seed, word16 nonce) {
  word32 i;

  for(i = 0; i < mL; ++i)
    PolyUniformGamma1(&v->at(i), seed, mL*nonce + i);
}


void Dilithium::PolyvecMatrixPointwiseMontgomery(polyvec *t, const std::vector<polyvec> *mat, const polyvec *v) {
  word32 i;

  for(i = 0; i < mK; ++i)
    PolyvecPointwiseAccMontgomery(&t->at(i), &mat->at(i), v, mL);
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

/*************************************************
* Name:        polyveck_caddq
*
* Description: For all coefficients of polynomials in vector of length K
*              add Q if coefficient is negative.
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void Dilithium::PolyvecCAddQ(polyvec *v, byte vectorLen) {
  word32 i;

  for(i = 0; i < vectorLen; ++i)
    PolyCAddQ(&v->at(i));
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
* Assumes v to be reduced by PolyvecReduce().
*
* Arguments:   - const polyvec *v: pointer to vector
*              - sword32 B: norm bound
*              - byte vectorLen: lenght of vector
*
* Returns 0 if norm of all polynomials is strictly smaller than B <= (mQ - 1)/8 and 1
* otherwise.
*/
int Dilithium::PolyvecChkNorm(const polyvec *v, sword32 bound, byte vectorLen)  {
  word32 i;

  for(i = 0; i < vectorLen; ++i) {
    if(PolyChkNorm(&v->at(i), bound)) {
      return 1;
    }
  }
    

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
* Subtract vectors of polynomials of length mK.
* No modular reduction is performed.
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
* reduction. Assumes input coefficients to be less than 2^{31-mD}.
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
//void Dilithium::PolyvecInvNttToMont(polyvec *v, byte vectorLen) {
  //word32 i;

 // for(i = 0; i < vectorLen; ++i)
   // PolyInvNttToMont(&v->at(i));
//}

/*
* For all coefficients a of polynomials in vector of length mK,
* compute a0, a1 such that a mod^+ mQ = a1*2^mD + a0
* with -2^{mD-1} < a0 <= 2^{mD-1}. Assumes coefficients to be
* standard representatives.
*
* Arguments:   - polyvec *v1: pointer to output vector of polynomials with
*                              coefficients a1
*              - polyvec *v0: pointer to output vector of polynomials with
*                              coefficients a0
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
* compute high and low bits a0, a1 such a mod^+ mQ = a1*mAlpha + a0
* with -mAlpha/2 < a0 <= mAlpha/2 except a1 = (mQ-1)/mAlpha where we
* set a1 = 0 and -mAlpha/2 <= a0 = a mod mQ - mQ < 0.
* Assumes coefficients to be standard representatives.
*
* Arguments:   - polyvec *v1: pointer to output vector of polynomials with
*                              coefficients a1
*              - polyvec *v0: pointer to output vector of polynomials with
*                              coefficients a0
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