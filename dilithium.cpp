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
* Implementation of H. Samples polynomial with TAU nonzero
* coefficients in {-1,1} using the output stream of
* SHAKE256(mu).
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const byte *mu: pointer to byte array containing mu
*/
void Dilithium::Challenge(poly *c, const byte *mu)
{
  word32 i, b, pos;
  word64 signs;
  std::vector<byte> bufVector (SHAKE256_RATE);
  std::vector<byte> *buf = &bufVector; 
  keccakState state;

  Shake256Init(&state);
  Shake256Absorb(&state, mu, mSeedBytes);
  Shake256Finalize(&state);
  Shake256SqueezeBlocks(buf->data(), 1, &state);

  signs = 0;
  for(i = 0; i < 8; ++i)
    signs |= (word64)buf->at(i) << 8*i;

  pos = 8;

  for(i = 0; i < mN; ++i)
    c->at(i) = 0;

  for(i = mN - mTau; i < mN; ++i) {
    do {
        if(pos >= SHAKE256_RATE) {
        Shake256SqueezeBlocks(buf->data(), 1, &state);
        pos = 0;
      }
      b = buf->at(pos++); 
    } while(b > i);

    c->at(i) = c->at(b);
    c->at(b) = 1 - 2*(signs & 1);
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
  std::vector<byte> seedBufVector(2*mSeedBytes + mCrhBytes);
  std::vector<byte> *seedBuf = &seedBufVector;
  //byte tr[mCrhBytes];
  std::vector<byte> trVector(mSeedBytes);
  std::vector<byte> *tr = &trVector;
  const byte *rho, *rhoprime, *key;

  polyvec matPolyvec (mL);
  std::vector<polyvec> mat (mK, matPolyvec);
  polyvec s1 (mL), s1hat (mL);
  polyvec s2 (mK), t1 (mK), t0 (mK);

  /* Get randomness for rho, rhoprime and key */
  RandomBytes(seedBuf->data(), mSeedBytes);
  SHAKE256 shake = SHAKE256(2*mSeedBytes + mCrhBytes);
  shake.Update(seedBuf->data(), mSeedBytes);
  shake.Final(seedBuf->data());
  rho = seedBuf->data();
  rhoprime = rho + mSeedBytes;
  key = rhoprime + mCrhBytes;
  /* Expand matrix */
  ExpandMat(&mat, rho);

  /* Sample short vectors s1 and s2 */
  PolyvecUniformEta(&s1, rhoprime, 0, mL);
  PolyvecUniformEta(&s2, rhoprime, mL, mK);
  /* Matrix-vector multiplication */
  s1hat = s1;
  PolyvecNtt(&s1hat, mL);
  PolyvecMatrixPointwiseMontgomery(&t1, &mat, &s1hat);
  PolyvecReduce(&t1, mK);
  PolyvecInvNttToMont(&t1, mK);
  

  /* Add error vector s2 */
  PolyvecAdd(&t1, &t1, &s2, mK);

  /* Extract t1 and write public key */
  PolyvecCAddQ(&t1, mK);
  PolyvecPower2Round(&t1, &t0, &t1, mK);
  PackPk(pk, rho, &t1);

  /* Compute CRH(rho, t1) and write secret key */
  SHAKE256 shakeTr = SHAKE256(mSeedBytes);
  shakeTr.Update(pk, mPublicKeyBytes);
  shakeTr.Final(tr->data());
  PackSk(sk, rho, tr->data(), key, &t0, &s1, &s2);
  return 0;
}

/*
* Computes signature.
*
* Arguments:   - byte *sig:   pointer to output signature (of length CRYPTO_BYTES)
*              - size_t *siglen: pointer to output length of signature
*              - byte *m:     pointer to message to be signed
*              - size_t mlen:    length of message
*              - byte *sk:    pointer to bit-packed secret key
*
* Returns 0 (success)
*/
int Dilithium::Signature(byte *sig, size_t *sigLen, const byte *m, size_t mLen, const byte *sk)
{
  word32 i, n;
  //byte seedbuf[2*mSeedBytes + 3*mCrhBytes];
  std::vector<byte> seedBufVector(3*mSeedBytes + 2*mCrhBytes);
  std::vector<byte> *seedBuf = &seedBufVector;
  byte *rho, *tr, *key, *mu, *rhoprime;
  word16 nonce = 0;
  poly challengePoly;
  polyvec matPolyvec (mL);
  std::vector<polyvec> mat (mK, matPolyvec);
  polyvec s1 (mL), y (mL), z (mL);
  polyvec t0 (mK), s2 (mK), w1 (mK), w0 (mK), h (mK);
  keccakState state;

  rho = seedBuf->data();
  tr = rho + mSeedBytes;
  key = tr + mSeedBytes;
  mu = key + mSeedBytes;
  rhoprime = mu + mCrhBytes;
  UnpackSk(rho, tr, key, &t0, &s1, &s2, sk);
  /* Compute CRH(tr, msg) */
  Shake256Init(&state);
  Shake256Absorb(&state, tr, mSeedBytes);
  Shake256Absorb(&state, m, mLen);
  Shake256Finalize(&state);
  Shake256Squeeze(mu, mCrhBytes, &state);
//#ifdef DILITHIUM_RANDOMIZED_SIGNING
//  RandomBytes(rhoprime, mCrhBytes);
//#else§
  SHAKE256 shake = SHAKE256(mCrhBytes);
  shake.Update(key, mSeedBytes + mCrhBytes);
  shake.Final(rhoprime);
//#endif

  /* Expand matrix and transform vectors */
  ExpandMat(&mat, rho);
  PolyvecNtt(&s1, mL);
  PolyvecNtt(&s2, mK);
  PolyvecNtt(&t0, mK);
  
  bool rejected = true;
  while (rejected) {
  
      /* Sample intermediate vector y */
    PolyvecUniformGamma1(&y, rhoprime, nonce++);

    /* Matrix-vector multiplication */
    z = y;
    PolyvecNtt(&z, mL);
    PolyvecMatrixPointwiseMontgomery(&w1, &mat, &z);
    PolyvecReduce(&w1, mK);
    PolyvecInvNttToMont(&w1, mK);
    

    /* Decompose w and call the random oracle */
    PolyvecCAddQ(&w1, mK);
    PolyvecDecompose(&w1, &w0, &w1, mK);



    Polyvecw1Pack(sig, &w1, mK);

    Shake256Init(&state);
    Shake256Absorb(&state, mu, mCrhBytes);
    Shake256Absorb(&state, sig, mK*mPolyw1PackedBytes);
    Shake256Finalize(&state);
    Shake256Squeeze(sig, mSeedBytes, &state);
    //TODO: mahdollinen virhekohta
    Challenge(&challengePoly, sig);

    PolyNtt(&challengePoly);

    /* Compute z, reject if it reveals secret */
    
    PolyvecPointwisePolyMontgomery(&z, &challengePoly, &s1, mL);
    PolyvecInvNttToMont(&z, mL);
    
    PolyvecAdd(&z, &z, &y, mL);
    PolyvecReduce(&z, mL);
    int ret = PolyvecChkNorm(&z, mGamma1 - mBeta, mL);
    if(ret) {
      continue;
    }
      
    /* Check that subtracting cs2 does not change high bits of w and low bits
    * do not reveal secret information */
    PolyvecPointwisePolyMontgomery(&h, &challengePoly, &s2, mK);
    PolyvecInvNttToMont(&h, mK);
    PolyvecSub(&w0, &w0, &h, mK);
    PolyvecReduce(&w0, mK);
    if(PolyvecChkNorm(&w0, mGamma2 - mBeta, mK)) {
      continue;
    }
      

    /* Compute hints for w1 */
    PolyvecPointwisePolyMontgomery(&h, &challengePoly, &t0, mK);
    PolyvecInvNttToMont(&h, mK);
    PolyvecReduce(&h, mK);
    if(PolyvecChkNorm(&h, mGamma2, mK)) {
      continue;
    }
      

    PolyvecAdd(&w0, &w0, &h, mK);
    n = PolyvecMakeHint(&h, &w0, &w1, mK);
    if(n > mOmega) {
      continue;
    }
    //Everything is ok, algorithm can proceed 
    rejected = false;
    break;

  }
  /* Write signature */
  PackSig(sig, sig, &z, &h);
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
  word32 i;
  //byte rho[mSeedBytes];
  std::vector<byte> bufVector(mK*mPolyw1PackedBytes);
  std::vector<byte> *buf = &bufVector;
  std::vector<byte> rhoVector(mSeedBytes);
  std::vector<byte> *rho = &rhoVector;
  std::vector<byte> muVector(mCrhBytes);
  std::vector<byte> *mu = &muVector;
  std::vector<byte> cVector(mSeedBytes);
  std::vector<byte> *c = &cVector;
  std::vector<byte> c2Vector(mSeedBytes);
  std::vector<byte> *c2 = &c2Vector;
  poly challengePoly;
  polyvec matPolyvec(mL);
  std::vector<polyvec> mat(mK, matPolyvec);
  polyvec z(mL);
  polyvec t1(mK), h(mK), w1(mK);
  keccakState state;
  if(sigLen != mBytes)
    return -1;

  UnpackPk(rho->data(), &t1, pk);
  if(UnpackSig(c->data(), &z, &h, sig))
    return -1;
  if(PolyvecChkNorm(&z, mGamma1 - mBeta, mL))
    return -1;

  /* Compute CRH(H(rho, t1), msg) */
  SHAKE256 shake = SHAKE256(mSeedBytes);
  shake.Update(pk, mPublicKeyBytes);
  shake.Final(mu->data());
  Shake256Init(&state);
  Shake256Absorb(&state, mu->data(), mSeedBytes);
  Shake256Absorb(&state, m, mLen);
  Shake256Finalize(&state);
  Shake256Squeeze(mu->data(), mCrhBytes, &state);

  /* Matrix-vector multiplication; compute Az - c2^dt1 */
  Challenge(&challengePoly, c->data());
  ExpandMat(&mat, rho->data());


  PolyvecNtt(&z, mL);
  PolyvecMatrixPointwiseMontgomery(&w1, &mat, &z);

  PolyNtt(&challengePoly);
  PolyvecShiftL(&t1, mK);
  PolyvecNtt(&t1, mK);
  PolyvecPointwisePolyMontgomery(&t1, &challengePoly, &t1, mK);

  PolyvecSub(&w1, &w1, &t1, mK);
  PolyvecReduce(&w1, mK);
  PolyvecInvNttToMont(&w1, mK);

  /* Reconstruct w1 */
  PolyvecCAddQ(&w1, mK);
  PolyvecUseHint(&w1, &w1, &h, mK);
  Polyvecw1Pack(buf->data(), &w1, mK);


  /* Call random oracle and verify challenge */
  Shake256Init(&state);
  Shake256Absorb(&state, mu->data(), mCrhBytes);
  Shake256Absorb(&state, buf->data(), mK*mPolyw1PackedBytes);
  Shake256Finalize(&state);
  Shake256Squeeze(c2 -> data(), mSeedBytes, &state);
  for(i = 0; i < mSeedBytes; ++i)
    if(c->at(i) != c2->at(i))
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
Dilithium::Dilithium(byte k, byte l, byte eta, byte tau, word16 beta, word32 gamma1, word32 gamma2, word16 omega, word16 polyEtaPacked, 
        word16 polyzPackedBytes, byte polyw1PackedBytes):
    mK(k),
    mL(l),
    mTau(tau),
    mEta(eta),
    mBeta(beta),
    mGamma1(gamma1),
    mGamma2(gamma2),
    mOmega(omega),
    mPolyEtaPackedBytes(polyEtaPacked),
    mPolyzPackedBytes(polyzPackedBytes),
    mPolyw1PackedBytes(polyw1PackedBytes),
    mPublicKeyBytes(mSeedBytes + k * mPolyt1PackedBytes),
    mSecretKeyBytes(3 * mSeedBytes + l * polyEtaPacked + k * polyEtaPacked + k * mPolyt0PackedBytes),
    mBytes(mSeedBytes + l * polyzPackedBytes + omega + k)
    {

    }






NAMESPACE_END