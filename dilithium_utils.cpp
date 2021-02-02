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

// The Keccak core function
extern void KeccakF1600(word64 *state);


/**
 * Get random bytes in a byte array the size of outLen.
 * Uses Crypto++'s AutoSeededRandomPool (osrng.h).
 */ 
void Dilithium::RandomBytes(byte *out, size_t outLen) {
  AutoSeededRandomPool rng;
  rng.GenerateBlock(out, outLen);
}


/* For finite field element a, compute a0, a1 such that
*  a mod Q = a1*2^D + a0 with -2^{D-1} < a0 <= 2^{D-1}.
*  Assumes a to be standard representative.
*  Arguments:  - word32 a: input element
*              - word32 *a0: pointer to output element Q + a0
*
* Returns a1.
*/
word32 Dilithium::Power2Round(word32 a, word32 *a0)  {
  sword32 t;

  /* Centralized remainder mod 2^D */
  t = a & ((1U << mD) - 1);
  t -= (1U << (mD-1)) + 1;
  t += (t >> 31) & (1U << mD);
  t -= (1U << (mD-1)) - 1;
  *a0 = mQ + t;
  a = (a - t) >> mD;
  return a;
}

/*  For finite field element a, compute high and low bits a0, a1 such
*   that a mod Q = a1*ALPHA + a0 with -ALPHA/2 < a0 <= ALPHA/2 except
*   if a1 = (Q-1)/ALPHA where we set a1 = 0 and
*   -ALPHA/2 <= a0 = a mod Q - Q < 0. Assumes a to be standard
*   representative.
*   Arguments: - word32 a: input element
*              - word32 *a0: pointer to output element Q + a0
*
* Returns a1.
*/
word32 Dilithium::Decompose(word32 a, word32 *a0) {
  sword32 t, u;

  /* Centralized remainder mod ALPHA */
  t = a & 0x7FFFF;
  t += (a >> 19) << 9;
  t -= mAlpha/2 + 1;
  t += (t >> 31) & mAlpha;
  t -= mAlpha/2 - 1;
  a -= t;

  /* Divide by ALPHA (possible to avoid) */
  u = a - 1;
  u >>= 31;
  a = (a >> 19) + 1;
  a -= u & 1;

  /* Border case */
  *a0 = mQ + t - (a >> 4);
  a &= 0xF;
  return a;
}

/* Compute hint bit indicating whether the low bits of the
*  input element overflow into the high bits. Inputs assumed to be
*  standard representatives.
*  Arguments:   - word32 a0: low bits of input element
*               - word32 a1: high bits of input element
*
* Returns 1 if high bits of a and b differ and 0 otherwise.
*/
word32 Dilithium::MakeHint(word32 a0, word32 a1) {
  if(a0 <= mGamma2 || a0 > mQ - mGamma2 || (a0 == mQ - mGamma2 && a1 == 0))
    return 0;

  return 1;
}

/* Correct high bits according to hint.
* Arguments:   - word32 a: input element
*              - word32 hint: hint bit
*
* Returns corrected high bits.
*/
word32 Dilithium::UseHint(word32 a, word32 hint) {
  word32 a0, a1;

  a1 = Decompose(a, &a0);
  if(hint == 0)
    return a1;
  else if(a0 > mQ)
    return (a1 + 1) & 0xF;
  else
    return (a1 - 1) & 0xF;

  //Note: This was commented out in the reference implementation
  /* If decompose does not divide out ALPHA:
  if(hint == 0)
    return a1;
  else if(a0 > Q)
    return (a1 + ALPHA) % (Q - 1);
  else
    return (a1 - ALPHA) % (Q - 1);
  */
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
* Bit-pack secret key sk = (rho, key, tr, s1, s2, t0).
*
* Arguments:   - byte *sk: pointer to output byte array
*              - const byte *rho: pointer to byte array containing rho
*              - const byte *key: pointer to byte array containing key
*              - const byte *tr: byte array containing tr
*              - const polyvec *s1: pointer to vector s1
*              - const polyvec *s2: pointer to vector s2
*              - const polyvec *t0: pointer to vector t0
*/
void Dilithium::PackSk(byte *sk, const byte *rho, const byte *key, const byte *tr, const polyvec *s1,
             const polyvec *s2, const polyvec *t0)
{
  word32 i;

  for(i = 0; i < mSeedBytes; ++i)
    sk[i] = rho[i];
  sk += mSeedBytes;

  for(i = 0; i < mSeedBytes; ++i)
    sk[i] = key[i];
  sk += mSeedBytes;

  for(i = 0; i < mCrhBytes; ++i)
    sk[i] = tr[i];
  sk += mCrhBytes;

  for(i = 0; i < mL; ++i)
    PolyEtaPack(sk + i*mPolyEtaPackedBytes, &s1->at(i));
  sk += mL*mPolyEtaPackedBytes;

  for(i = 0; i < mK; ++i)
    PolyEtaPack(sk + i*mPolyEtaPackedBytes, &s2->at(i));
  sk += mK*mPolyEtaPackedBytes;

  for(i = 0; i < mK; ++i)
    Polyt0Pack(sk + i*mPolyt0PackedBytes, &t0->at(i));
}

/* Unpack secret key sk = (rho, key, tr, s1, s2, t0).
* Arguments:   - const byte *rho: pointer to output byte array for rho
*              - const byte *key: pointer to output byte array for key
*              - const byte *tr: pointer to output byte array for tr
*              - const polyvec *s1: pointer to output vector s1
*              - const polyvec *s2: pointer to output vector s2
*              - const polyvec *r0: pointer to output vector t0
*              - byte *sk: byte array containing bit-packed sk
*/
void Dilithium::UnpackSk(byte *rho, byte *key, byte *tr, 
                          polyvec *s1, polyvec *s2, polyvec *t0, const byte *sk)
{
  word32 i;

  for(i = 0; i < mSeedBytes; ++i)
    rho[i] = sk[i];
  sk += mSeedBytes;

  for(i = 0; i < mSeedBytes; ++i)
    key[i] = sk[i];
  sk += mSeedBytes;

  for(i = 0; i < mCrhBytes; ++i)
    tr[i] = sk[i];
  sk += mCrhBytes;
  
  for(i=0; i < mL; ++i)
    PolyEtaUnpack(&s1->at(i), sk + i*mPolyEtaPackedBytes);
  sk += mL*mPolyEtaPackedBytes;

  for(i=0; i < mK; ++i)
    PolyEtaUnpack(&s2->at(i), sk + i*mPolyEtaPackedBytes);
  sk += mK*mPolyEtaPackedBytes;
  
  for(i=0; i < mK; ++i)
    Polyt0Unpack(&t0->at(i), sk + i*mPolyt0PackedBytes);
}

/* Bit-pack signature sig = (z, h, c).
* Arguments:   - byte *sig: pointer to output byte array
*              - const polyvec *z: pointer to vector z
*              - const polyvec *h: pointer to hint vector h
*              - const poly *c: pointer to challenge polynomial
*/
void Dilithium::PackSig(byte *sig, const polyvec *z, const polyvec *h, const poly *c)
{
  word32 i, j, k;
  word64 signs, mask;

  for(i = 0; i < mL; ++i)
    PolyzPack(sig + i*mPolyzPackedBytes, &z->at(i));
  sig += mL*mPolyzPackedBytes;

  /* Encode h */
  k = 0;
  for(i = 0; i < mK; ++i) {
    for(j = 0; j < mN; ++j)
      if(h->at(i).at(j) != 0)
        sig[k++] = j;

    sig[mOmega + i] = k;
  }
  while(k < mOmega) sig[k++] = 0;
  sig += mOmega + mK;

  /* Encode c */
  signs = 0;
  mask = 1;
  for(i = 0; i < mN/8; ++i) {
    sig[i] = 0;
    for(j = 0; j < 8; ++j) {
      if(c->at(8*i+j) != 0) {
        sig[i] |= (1U << j);
        if(c->at(8*i+j) == (mQ - 1)) signs |= mask;
        mask <<= 1;
      }
    }
  }
  sig += mN/8;
  for(i = 0; i < 8; ++i)
    sig[i] = signs >> 8*i;
}


/* Unpack signature sig = (z, h, c).
* Arguments:   - polyvec *z: pointer to output vector z
*              - polyvec *h: pointer to output hint vector h
*              - poly *c: pointer to output challenge polynomial
*              - const byte *sig: pointer to byte array containing
*                bit-packed signature
*
* Returns 1 in case of malformed signature; otherwise 0.
*/
int Dilithium::UnpackSig(polyvec *z, polyvec *h, poly *c, const byte *sig)
{
  word32 i, j, k;
  word64 signs;

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

  sig += mOmega + mK;

  /* Decode c */
  for(i = 0; i < mN; ++i)
    c->at(i) = 0;

  signs = 0;
  for(i = 0; i < 8; ++i)
    signs |= (uint64_t)sig[mN/8+i] << 8*i;

  /* Extra sign bits are zero for strong unforgeability */
  if(signs >> 60)
    return 1;

  for(i = 0; i < mN/8; ++i) {
    for(j = 0; j < 8; ++j) {
      if((sig[i] >> j) & 0x01) {
        c->at(8*i+j) = 1;
        c->at(8*i+j) ^= -(signs & 1) & (1 ^ (mQ-1));
        signs >>= 1;
      }
    }
  }

  return 0;
}


/*
* For finite field element a with 0 <= a <= Q*2^32,
* compute r \equiv a*2^{-32} (mod Q) such that 0 <= r < 2*Q.
* Arguments:   - word64: finite field element a
*
* Returns r.
*/
word32 Dilithium::MontgomeryReduce(word64 a) {
  word64 t;

  t = a * mQInv;
  t &= ((word64)1 << 32) - 1;
  t *= mQ;
  t = a + t;
  t >>= 32;
  return t;
}

/* 
* For finite field element a, compute r \equiv a (mod Q)
* such that 0 <= r < 2*Q.
* Arguments:   - word32: finite field element a
*
* Returns r.
*/
word32 Dilithium::Reduce32(word32 a) {
  word32 t;

  t = a & 0x7FFFFF;
  a >>= 23;
  t += (a << 13) - a;
  return t;
}

/*
* Subtract Q if input coefficient is bigger than Q.
*
* Arguments:   - word32: finite field element a
*
* Returns r.
*/
word32 Dilithium::CSubQ(word32 a) {
  a -= mQ;
  a += ((sword32)a >> 31) & mQ;
  return a;
}

/* 
*  For finite field element a, compute standard
* representative r = a mod Q.
* Arguments:   - word32: finite field element a
*
* Returns r.
*/
word32 Dilithium::Freeze(word32 a) {
  a = Reduce32(a);
  a = CSubQ(a);
  return a;
}


/*
 * Initiate Shake-128 Stream.
 * Arguments: - *keccakState 
 * 
 * 
 */
void Dilithium::Shake128StreamInit(keccakState *state, const byte *seed, word16 nonce)
{
  byte t[2];
  t[0] = nonce;
  t[1] = nonce >> 8;

  Shake128Init(state);
  Shake128Absorb(state, seed, mSeedBytes);
  Shake128Absorb(state, t, 2);
  Shake128Finalize(state);
}

void Dilithium::Shake256StreamInit(keccakState *state, const byte *seed, word16 nonce)
{
  byte t[2];
  t[0] = nonce;
  t[1] = nonce >> 8;

  Shake256Init(state);
  Shake256Absorb(state, seed, mCrhBytes);
  Shake256Absorb(state, t, 2);
  Shake256Finalize(state);
}


#define NROUNDS 24
#define ROL(a, offset) ((a << offset) ^ (a >> (64-offset)))

/*
* Load 8 bytes into word64 in little-endian order
*
* Arguments:   - const byte *x: pointer to input byte array
*
* Returns the loaded 64-bit unsigned integer
**************************************************/
static word64 Load64(const byte x[8]) {
  word32 i;
  word64 r = 0;

  for(i=0;i<8;i++)
    r |= (word64)x[i] << 8*i;

  return r;
}

//TODO: JATKA TÄSTÄ MA

/*
* Store a 64-bit integer to array of 8 bytes in little-endian order
*
* Arguments:   - byte *x: pointer to the output byte array (allocated)
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
* Arguments:   - keccakState *state: pointer to Keccak state
*/
void Dilithium::KeccakInit(keccakState *state)
{
  word32 i;
  for(i=0;i<25;i++)
    state->s[i] = 0;
  state->pos = 0;
}

/*
* Absorb step of Keccak; incremental.
*
* Arguments:   - word64 *s:      pointer to Keccak state
*              - word32 r:   rate in bytes (e.g., 168 for SHAKE128)
*              - word32 pos: position in current block to be absorbed
*              - const byte *m: pointer to input to be absorbed into s
*              - size_t mlen:      length of input in bytes
*
* Returns new position pos in current block
*/
word32 Dilithium::KeccakAbsorb(word64 s[25], word32 r, word32 pos, const byte *m, size_t mlen)
{
  word32 i;
  byte t[8] = {0};

  if(pos & 7) {
    i = pos & 7;
    while(i < 8 && mlen > 0) {
      t[i++] = *m++;
      mlen--;
      pos++;
    }
    s[(pos-i)/8] ^= Load64(t);
  }

  if(pos && mlen >= r-pos) {
    for(i=0;i<(r-pos)/8;i++)
      s[pos/8+i] ^= Load64(m+8*i);
    m += r-pos;
    mlen -= r-pos;
    pos = 0;
    KeccakF1600(s);
  }

  while(mlen >= r) {
    for(i=0;i<r/8;i++)
      s[i] ^= Load64(m+8*i);
    m += r;
    mlen -= r;
    KeccakF1600(s);
  }

  for(i=0;i<mlen/8;i++)
    s[pos/8+i] ^= Load64(m+8*i);
  m += 8*i;
  mlen -= 8*i;
  pos += 8*i;

  if(mlen) {
    for(i=0;i<8;i++)
      t[i] = 0;
    for(i=0;i<mlen;i++)
      t[i] = m[i];
    s[pos/8] ^= Load64(t);
    pos += mlen;
  }

  return pos;
}

/*
* Finalize absorb step.
*
* Arguments:   - word64 *s:      pointer to Keccak state
*              - word32 r:   rate in bytes (e.g., 168 for SHAKE128)
*              - word32 pos: position in current block to be absorbed
*              - byte p:        domain separation byte
*/
void Dilithium::KeccakFinalize(word64 s[25], word32 r, word32 pos, byte p)
{
  word32 i,j;

  i = pos >> 3;
  j = pos & 7;
  s[i] ^= (uint64_t)p << 8*j;
  s[r/8-1] ^= 1ULL << 63;
}

/*
* Squeeze step of Keccak. Squeezes full blocks of r bytes each.
* Modifies the state. Can be called multiple times to keep
* squeezing, i.e., is incremental. Assumes zero bytes of current
* block have already been squeezed.
*
* Arguments:   - byte *out:   pointer to output blocks
*              - size_t nBlocks: number of blocks to be squeezed (written to out)
*              - word64 *s:    pointer to input/output Keccak state
*              - word32 r: rate in bytes (e.g., 168 for SHAKE128)
*/
void Dilithium::KeccakSqueezeBlocks(byte *out, size_t nBlocks, word64 s[25], word32 r)
{
  word32 i;

  while(nBlocks > 0) {
    KeccakF1600(s);
    for(i=0;i<r/8;i++)
      Store64(out + 8*i, s[i]);
    out += r;
    nBlocks--;
  }
}

/*
* Squeeze step of Keccak. Squeezes arbitratrily many bytes.
* Modifies the state. Can be called multiple times to keep
* squeezing, i.e., is incremental.
*
* Arguments:   - byte *out:     pointer to output
*              - size_t outLen:    number of bytes to be squeezed (written to out)
*              - word64 *s:      pointer to input/output Keccak state
*              - word32 r:   rate in bytes (e.g., 168 for SHAKE128)
*              - word32 pos: number of bytes in current block already squeezed
*
* Returns new position pos in current block
*/
unsigned int Dilithium::KeccakSqueeze(byte *out, size_t outLen, word64 s[25], word32 r, word32 pos)
{
  word32 i;
  byte t[8];

  if(pos & 7) {
    Store64(t,s[pos/8]);
    i = pos & 7;
    while(i < 8 && outLen > 0) {
      *out++ = t[i++];
      outLen--;
      pos++;
    }
  }

  if(pos && outLen >= r-pos) {
    for(i=0;i<(r-pos)/8;i++)
      Store64(out+8*i,s[pos/8+i]);
    out += r-pos;
    outLen -= r-pos;
    pos = 0;
  }

  while(outLen >= r) {
    KeccakF1600(s);
    for(i=0;i<r/8;i++)
      Store64(out+8*i,s[i]);
    out += r;
    outLen -= r;
  }

  if(!outLen)
    return pos;
  else if(!pos)
    KeccakF1600(s);

  for(i=0;i<outLen/8;i++)
    Store64(out+8*i,s[pos/8+i]);
  out += 8*i;
  outLen -= 8*i;
  pos += 8*i;

  Store64(t,s[pos/8]);
  for(i=0;i<outLen;i++)
    out[i] = t[i];
  pos += outLen;
  return pos;
}

/*
* Initilizes Keccak state for use as SHAKE128 XOF
*
* Arguments:   - keccakState *state: pointer to (uninitialized)
*                                     Keccak state
*/
void Dilithium::Shake128Init(keccakState *state)
{
  KeccakInit(state);
}

/*
* Absorb step of the SHAKE128 XOF; incremental.
*
* Arguments:   - keccakState *state: pointer to (initialized) output
*                                     Keccak state
*              - const byte *in:   pointer to input to be absorbed into s
*              - size_t inLen:        length of input in bytes
*/
void Dilithium::Shake128Absorb(keccakState *state, const byte *in, size_t inLen)
{
  state->pos = KeccakAbsorb(state->s, SHAKE128_RATE, state->pos, in, inLen);
}

/*
* Finalize absorb step of the SHAKE128 XOF.
*
* Arguments:   - keccakState *state: pointer to Keccak state
*/
void Dilithium::Shake128Finalize(keccakState *state)
{
  KeccakFinalize(state->s, SHAKE128_RATE, state->pos, 0x1F);
  state->pos = 0;
}

/*
* Squeeze step of SHAKE128 XOF. Squeezes full blocks of
* SHAKE128_RATE bytes each. Can be called multiple times
* to keep squeezing. Assumes zero bytes of current block
* have already been squeezed (state->pos = 0).
*
* Arguments:   - byte *out:    pointer to output blocks
*              - size_t nbBocks:  number of blocks to be squeezed (written to output)
*              - keccakState *s: pointer to input/output Keccak state
*/
void Dilithium::Shake128SqueezeBlocks(byte *out, size_t nBlocks, keccakState *state)
{
  KeccakSqueezeBlocks(out, nBlocks, state->s, SHAKE128_RATE);
}

/*
* Squeeze step of SHAKE128 XOF. Squeezes arbitraily many
* bytes. Can be called multiple times to keep squeezing.
*
* Arguments:   - byte *out:    pointer to output blocks
*              - size_t outlen :  number of bytes to be squeezed
*                                 (written to output)
*              - keccakState *s: pointer to input/output Keccak state
**************************************************/
void Dilithium::Shake128Squeeze(byte *out, size_t outLen, keccakState *state)
{
  state->pos = KeccakSqueeze(out, outLen, state->s, SHAKE128_RATE, state->pos);
}

/*
* Initilizes Keccak state for use as SHAKE256 XOF
*
* Arguments:   - keccakState *state: pointer to (uninitialized) Keccak state
*/
void Dilithium::Shake256Init(keccakState *state)
{
  KeccakInit(state);
}

/*
* Absorb step of the SHAKE256 XOF; incremental.
*
* Arguments:   - keccakState *state: pointer to (initialized) output Keccak state
*              - const byte *in:   pointer to input to be absorbed into s
*              - size_t inLen:        length of input in bytes
*/
void Dilithium::Shake256Absorb(keccakState *state, const byte *in, size_t inLen)
{
  state->pos = KeccakAbsorb(state->s, SHAKE256_RATE, state->pos, in, inLen);
}

/*
* Finalize absorb step of the SHAKE256 XOF.
*
* Arguments:   - keccakState *state: pointer to Keccak state
*/
void Dilithium::Shake256Finalize(keccakState *state)
{
  KeccakFinalize(state->s, SHAKE256_RATE, state->pos, 0x1F);
  state->pos = 0;
}

/*
* Squeeze step of SHAKE256 XOF. Squeezes full blocks of
* SHAKE256_RATE bytes each. Can be called multiple times
* to keep squeezing. Assumes zero bytes of current block
* have already been squeezed (state->pos = 0).
*
* Arguments:   - byte *out:    pointer to output blocks
*              - size_t nBlocks:  number of blocks to be squeezed
*                                 (written to output)
*              - keccakState *s: pointer to input/output Keccak state
*/
void Dilithium::Shake256SqueezeBlocks(byte *out, size_t nBlocks, keccakState *state)
{
  KeccakSqueezeBlocks(out, nBlocks, state->s, SHAKE256_RATE);
}

/*
* Squeeze step of SHAKE256 XOF. Squeezes arbitraily many
* bytes. Can be called multiple times to keep squeezing.
*
* Arguments:   - byte *out:    pointer to output blocks
*              - size_t outLen :  number of bytes to be squeezed
*                                 (written to output)
*              - keccakState *s: pointer to input/output Keccak state
*/
void Dilithium::Shake256Squeeze(byte *out, size_t outLen, keccakState *state)
{
  state->pos = KeccakSqueeze(out, outLen, state->s, SHAKE256_RATE, state->pos);
}


NAMESPACE_END