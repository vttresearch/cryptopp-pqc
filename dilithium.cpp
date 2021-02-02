/* The C++ file for CRYSTALS-Dilithium. Adapted by Julius Hekkala
from the public-domain reference
implementation of Dilithium by the CRYSTALS team 
(https://github.com/pq-crystals/dilithium)
 */

#include "pch.h"
#include "dilithium.h"
#include "shake.h"

NAMESPACE_BEGIN(CryptoPP)


/*
* Implementation of H. Samples polynomial with 60 nonzero
* coefficients in {-1,1} using the output stream of
* SHAKE256(mu|w1).
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const byte *mu: pointer to byte array containing mu
*              - const polyvec *w1: pointer to vector w1
*/
void Dilithium::Challenge(poly *c, const byte *mu, const polyvec *w1)
{
  word32 i, b, pos;
  word64 signs;
  word32 bufLen = mCrhBytes + mK*mPolyw1PackedBytes;
  //byte buf[bufLen];
  std::vector<byte> bufVector (bufLen);
  std::vector<byte> *buf = &bufVector; 
  keccakState state;

  for(i = 0; i < mCrhBytes; ++i)
    buf->at(i) = mu[i];
  for(i = 0; i < mK; ++i)
    Polyw1Pack(buf->data() + mCrhBytes + i*mPolyw1PackedBytes, &w1->at(i));

  Shake256Init(&state);
  Shake256Absorb(&state, buf->data(), bufLen);
  Shake256Finalize(&state);
  Shake256SqueezeBlocks(buf->data(), 1, &state);



  signs = 0;
  for(i = 0; i < 8; ++i)
    signs |= (word64)buf->at(i) << 8*i;

  pos = 8;

  for(i = 0; i < mN; ++i)
    c->at(i) = 0;

  for(i = 196; i < 256; ++i) {
    do {
      if(pos >= SHAKE256_RATE) {
        Shake256SqueezeBlocks(buf->data(), 1, &state);
        pos = 0;
      }

      b = buf->at(pos++);
    } while(b > i);

    c->at(i) = c->at(b);
    c->at(b) = 1;
    c->at(b) ^= -((word32)signs & 1) & (1 ^ (mQ-1));
    signs >>= 1;
  }
}

/*
* Generates public and private key.
*
* Arguments:   - byte *pk: pointer to output public key (allocated
*                             array of PUBLICKEYBYTES bytes)
*              - byte *sk: pointer to output private key (allocated
*                             array of SECRETKEYBYTES bytes)
*
* Returns 0 (success)
*/
int Dilithium::Keypair(byte *pk, byte *sk) {
  unsigned int i;
  //byte seedbuf[3*mSeedBytes];
  std::vector<byte> seedBufVector(3*mSeedBytes);
  std::vector<byte> *seedBuf = &seedBufVector;
  //byte tr[mCrhBytes];
  std::vector<byte> trVector(mCrhBytes);
  std::vector<byte> *tr = &trVector;
  const byte *rho, *rhoprime, *key;
  word16 nonce = 0;

  polyvec matPolyvec (mL);
  std::vector<polyvec> mat (mK, matPolyvec);
  polyvec s1 (mL), s1hat (mL);
  polyvec s2 (mK), t1 (mK), t0 (mK);

  /* Get randomness for rho, rhoprime and key */
  RandomBytes(seedBuf->data(), 3*mSeedBytes);
  rho = seedBuf->data();
  rhoprime = seedBuf->data() + mSeedBytes;
  key = seedBuf->data() + 2*mSeedBytes;

  /* Expand matrix */
  ExpandMat(&mat, rho);

  /* Sample short vectors s1 and s2 */
  for(i = 0; i < mL; ++i)
    PolyUniformEta(&s1.at(i), rhoprime, nonce++);
  for(i = 0; i < mK; ++i)
    PolyUniformEta(&s2.at(i), rhoprime, nonce++);

  /* Matrix-vector multiplication */
  s1hat = s1;
  PolyvecNtt(&s1hat, mL);
  for(i = 0; i < mK; ++i) {
    PolyvecPointwiseAccMontgomery(&t1.at(i), &mat.at(i), &s1hat, mL);
    PolyReduce(&t1.at(i));
    PolyInvNttToMont(&t1.at(i));
  }

  /* Add error vector s2 */
  PolyvecAdd(&t1, &t1, &s2, mK);

  /* Extract t1 and write public key */
  PolyvecFreeze(&t1, mK);
  PolyvecPower2Round(&t1, &t0, &t1, mK);
  PackPk(pk, rho, &t1);

  /* Compute CRH(rho, t1) and write secret key */
  SHAKE256 shake = SHAKE256(mCrhBytes);
  shake.Update(pk, mPublicKeyBytes);
  shake.Final(tr->data());
  PackSk(sk, rho, key, tr->data(), &s1, &s2, &t0);

  return 0;
}

/*
* Computes signature.
*
* Arguments:   - byte *sig:   pointer to output signature (of length CRYPTO_BYTES)
*              - size_t *siglen: pointer to output length of signed message
*              - byte *m:     pointer to message to be signed
*              - size_t mlen:    length of message
*              - byte *sk:    pointer to bit-packed secret key
*
* Returns 0 (success)
*/
int Dilithium::Signature(byte *sig, size_t *sigLen, const byte *m, size_t mLen, const byte *sk)
{
  unsigned int i, n;
  //byte seedbuf[2*mSeedBytes + 3*mCrhBytes];
  std::vector<byte> seedBufVector(2*mSeedBytes + 3*mCrhBytes);
  std::vector<byte> *seedBuf = &seedBufVector;
  byte *rho, *tr, *key, *mu, *rhoprime;
  word16 nonce = 0;
  poly c, chat;
  polyvec matPolyvec (mL);
  std::vector<polyvec> mat (mK, matPolyvec);
  polyvec s1 (mL), y (mL), z (mL);
  polyvec t0 (mK), s2 (mK), w1 (mK), w0 (mK), h (mK);
  keccakState state;

  rho = seedBuf->data();
  tr = rho + mSeedBytes;
  key = tr + mCrhBytes;
  mu = key + mSeedBytes;
  rhoprime = mu + mCrhBytes;
  UnpackSk(rho, key, tr, &s1, &s2, &t0, sk);
  
  /* Compute CRH(tr, msg) */
  Shake256Init(&state);
  Shake256Absorb(&state, tr, mCrhBytes);
  Shake256Absorb(&state, m, mLen);
  Shake256Finalize(&state);
  Shake256Squeeze(mu, mCrhBytes, &state);

#ifdef DILITHIUM_RANDOMIZED_SIGNING
  RandomBytes(rhoprime, mCrhBytes);
#else
  SHAKE256 shake = SHAKE256(mCrhBytes);
  shake.Update(key, mSeedBytes + mCrhBytes);
  shake.Final(rhoprime);
#endif

  /* Expand matrix and transform vectors */
  ExpandMat(&mat, rho);
  PolyvecNtt(&s1, mL);
  PolyvecNtt(&s2, mK);
  PolyvecNtt(&t0, mK);
  
  bool rejected = true;
  while (rejected) {
  
      /* Sample intermediate vector y */
    for(i = 0; i < mL; ++i)
      PolyUniformGamma1m1(&y.at(i), rhoprime, nonce++);

    /* Matrix-vector multiplication */
    z = y;
    PolyvecNtt(&z, mL);
    for(i = 0; i < mK; ++i) {
      PolyvecPointwiseAccMontgomery(&w1.at(i), &mat.at(i), &z, mL);
      PolyReduce(&w1.at(i));
      PolyInvNttToMont(&w1.at(i));
    }

    /* Decompose w and call the random oracle */
    PolyvecCsubQ(&w1, mK);
    PolyvecDecompose(&w1, &w0, &w1, mK);
    Challenge(&c, mu, &w1);
    chat = c;
    PolyNtt(&chat);

    /* Compute z, reject if it reveals secret */
    for(i = 0; i < mL; ++i) {
      PolyPointwiseMontgomery(&z.at(i), &chat, &s1.at(i));
      PolyInvNttToMont(&z.at(i));
    }
    PolyvecAdd(&z, &z, &y, mL);
    PolyvecFreeze(&z, mL);
    if(PolyvecChkNorm(&z, mGamma1 - mBeta, mL)) {
      continue;
    }
      
    /* Check that subtracting cs2 does not change high bits of w and low bits
    * do not reveal secret information */
    for(i = 0; i < mK; ++i) {
      PolyPointwiseMontgomery(&h.at(i), &chat, &s2.at(i));
      PolyInvNttToMont(&h.at(i));
    }
    PolyvecSub(&w0, &w0, &h, mK);
    PolyvecFreeze(&w0, mK);
    if(PolyvecChkNorm(&w0, mGamma2 - mBeta, mK)) {
      continue;
    }
      

    /* Compute hints for w1 */
    for(i = 0; i < mK; ++i) {
      PolyPointwiseMontgomery(&h.at(i), &chat, &t0.at(i));
      PolyInvNttToMont(&h.at(i));
    }
    PolyvecCsubQ(&h, mK);
    if(PolyvecChkNorm(&h, mGamma2, mK)) {
      continue;
    }
      

    PolyvecAdd(&w0, &w0, &h, mK);
    PolyvecCsubQ(&w0, mK);
    n = PolyvecMakeHint(&h, &w0, &w1, mK);
    if(n > mOmega) {
      continue;
    }
      
    //Everything is ok, algorithm can proceed 
    rejected = false;
    break;

  }
  
  /* Write signature */
  PackSig(sig, &z, &h, &c);
  *sigLen = mBytes;
  return 0;
}

/*
* Compute signed message.
*
* Arguments:   - byte *sm: pointer to output signed message (allocated
*                             array with CRYPTO_BYTES + mlen bytes),
*                             can be equal to m
*              - size_t *smLen: pointer to output length of signed
*                               message
*              - const byte *m: pointer to message to be signed
*              - size_t mLen: length of message
*              - const byte *sk: pointer to bit-packed secret key
*
* Returns 0 (success)
*/
int Dilithium::Sign(byte *sm, size_t *smLen, const byte *m, size_t mLen, const byte *sk)
{
  size_t i;

  for(i = 0; i < mLen; ++i)
    sm[mBytes + mLen - 1 - i] = m[mLen - 1 - i];
  Signature(sm, smLen, sm + mBytes, mLen, sk);
  *smLen += mLen;
  return 0;
}

/*
* Verifies signature.
*
* Arguments:   - byte *m: pointer to input signature
*              - size_t sigLen: length of signature
*              - const byte *m: pointer to message
*              - size_t mLen: length of message
*              - const byte *pk: pointer to bit-packed public key
*
* Returns 0 if signature could be verified correctly and -1 otherwise
*/
int Dilithium::Verify(const byte *sig, size_t sigLen, const byte *m, size_t mLen, const byte *pk)
{
  unsigned int i;
  //byte rho[mSeedBytes];
  std::vector<byte> rhoVector(mSeedBytes);
  std::vector<byte> *rho = &rhoVector;
  std::vector<byte> muVector(mCrhBytes);
  std::vector<byte> *mu = &muVector;
  poly c, cp;
  polyvec matPolyvec(mL);
  std::vector<polyvec> mat(mK, matPolyvec);
  polyvec z(mL);
  polyvec t1(mK), h(mK), w1(mK);
  keccakState state;

  if(sigLen != mBytes)
    return -1;

  UnpackPk(rho->data(), &t1, pk);
  if(UnpackSig(&z, &h, &c, sig))
    return -1;
  if(PolyvecChkNorm(&z, mGamma1 - mBeta, mL))
    return -1;

  /* Compute CRH(CRH(rho, t1), msg) */
  SHAKE256 shake = SHAKE256(mCrhBytes);
  shake.Update(pk, mPublicKeyBytes);
  shake.Final(mu->data());
  
  //TODO: VOiko vaihtaa CryptoPP:n shakella
  Shake256Init(&state);
  Shake256Absorb(&state, mu->data(), mCrhBytes);
  Shake256Absorb(&state, m, mLen);
  Shake256Finalize(&state);
  Shake256Squeeze(mu->data(), mCrhBytes, &state);

  /* Matrix-vector multiplication; compute Az - c2^dt1 */
  ExpandMat(&mat, rho->data());

  PolyvecNtt(&z, mL);
  for(i = 0; i < mK ; ++i)
    PolyvecPointwiseAccMontgomery(&w1.at(i), &mat.at(i), &z, mL);

  cp = c;
  PolyNtt(&cp);
  PolyvecShiftL(&t1, mK);
  PolyvecNtt(&t1, mK);
  for(i = 0; i < mK; ++i)
    PolyPointwiseMontgomery(&t1.at(i), &cp, &t1.at(i));

  PolyvecSub(&w1, &w1, &t1, mK);
  PolyvecReduce(&w1, mK);
  PolyvecInvNttToMont(&w1, mK);

  /* Reconstruct w1 */
  PolyvecCsubQ(&w1, mK);
  PolyvecUseHint(&w1, &w1, &h, mK);

  /* Call random oracle and verify challenge */
  Challenge(&cp, mu->data(), &w1);
  for(i = 0; i < mN; ++i)
    if(c.at(i) != cp.at(i))
      return -1;

  return 0;
}

/*
* Verify signed message.
*
* Arguments:   - byte *m: pointer to output message (allocated
*                            array with smlen bytes), can be equal to sm
*              - size_t *mLen: pointer to output length of message
*              - const byte *sm: pointer to signed message
*              - size_t smLen: length of signed message
*              - const byte *pk: pointer to bit-packed public key
*
* Returns 0 if signed message could be verified correctly and -1 otherwise
*/
int Dilithium::Open(byte *m, size_t *mLen, const byte *sm, size_t smLen, const byte *pk)
{
  size_t i;

  bool badSignature = false;

  if(smLen < mBytes)
    badSignature = true;

  if (!badSignature) {
    *mLen = smLen - mBytes;
    if(Verify(sm, mBytes, sm + mBytes, *mLen, pk))
      badSignature =  true;
    else {
      /* All good, copy msg, return 0 */
      for(i = 0; i < *mLen; ++i)
        m[i] = sm[mBytes + i];
      return 0;
    }

  }

  /* Signature verification failed */
  *mLen = -1;
  for(i = 0; i < smLen; ++i)
    m[i] = 0;

  return -1;
}


//Constructor
Dilithium::Dilithium(byte mode, byte eta, word16 beta, word16 omega, word16 polyEtaPacked,
                        word16 publicKeyBytes, word16 secretKeyBytes, word16 bytes) :
    mK(2 + mode),
    mL(1 + mode),
    mEta(eta),
    mBeta(beta),
    mOmega(omega),
    mPolyEtaPackedBytes(polyEtaPacked),
    mPublicKeyBytes(publicKeyBytes),
    mSecretKeyBytes(secretKeyBytes),
    mBytes(bytes)
    {

    }






NAMESPACE_END