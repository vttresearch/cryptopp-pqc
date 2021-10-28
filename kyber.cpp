/* The C++ file for CRYSTALS-Kyber. Adapted from the reference
implementation of kyber by the CRYSTALS team 
(https://github.com/pq-crystals/kyber) */

#include "pch.h"
#include "kyber.h"
#include "misc.h"

NAMESPACE_BEGIN(CryptoPP)

// Generates public and private key
// for CCA-secure Kyber key encapsulation mechanism
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
sword32 Kyber<T_K, T_C, T_C2, T_ETA1>::KemKeypair(byte  *pk, byte *sk)
{
  size_t i; 
  IndCpaKeypair(pk, sk);
  for(i=0;i<INDCPA_PUBLICKEYBYTES;i++) {
    sk[i+INDCPA_SECRETKEYBYTES] = pk[i];
  }
  
  Sha3_256(sk+SECRETKEYBYTES-2*SYMBYTES, pk, PUBLICKEYBYTES);
  /* Value z for pseudo-random output on reject */
 
  RandomBytes(sk+SECRETKEYBYTES-SYMBYTES, SYMBYTES);
  return 0;
}

// Generates cipher text and shared
// secret for given public key
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
sword32 Kyber<T_K, T_C, T_C2, T_ETA1>::KemEnc(byte *ct, byte *ss, const byte *pk)
{
  byte buf[2*SYMBYTES];
  /* Will contain key, coins */
  byte kr[2*SYMBYTES];
  
  RandomBytes(buf, SYMBYTES);
  /* Don't release system RNG output */
  Sha3_256(buf, buf, SYMBYTES);

  /* Multitarget countermeasure for coins + contributory KEM */
  Sha3_256(buf + SYMBYTES, pk, PUBLICKEYBYTES);
  
  Sha3_512(kr, buf, 2 * SYMBYTES);
  /* coins are in kr+KYBER_SYMBYTES */
  IndCpaEnc(ct, buf, pk, kr+ SYMBYTES); 

  /* overwrite coins in kr with H(c) */
  Sha3_256(kr + SYMBYTES, ct, CIPHERTEXTBYTES);
  /* hash concatenation of pre-k and H(c) to k */
  Shake256(ss, 32, kr, 2* SYMBYTES);
  return 0;
}


// Generates shared secret for given
// cipher text and private key
// On failure, ss will contain a pseudo-random value.
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
sword32 Kyber<T_K, T_C, T_C2, T_ETA1>::KemDec(byte *ss, const byte *ct, const byte *sk)
{
  size_t i;
  sword32 fail;
  byte buf[2*SYMBYTES];
  /* Will contain key, coins */
  byte kr[2*SYMBYTES];
  byte cmp[CIPHERTEXTBYTES];
  const byte *pk = sk+INDCPA_SECRETKEYBYTES;

  IndCpaDec(buf, ct, sk);

  /* Multitarget countermeasure for coins + contributory KEM */
  for(i=0;i<SYMBYTES;i++) {
    buf[SYMBYTES+i] = sk[SECRETKEYBYTES-2*SYMBYTES+i];
  }
  Sha3_512(kr, buf, 2*SYMBYTES);

  /* coins are in kr+KYBER_SYMBYTES */
  IndCpaEnc(cmp, buf, pk, kr+SYMBYTES);
  fail = !(VerifyBufsEqual(ct, cmp, CIPHERTEXTBYTES));

  /* overwrite coins in kr with H(c) */
  Sha3_256(kr + SYMBYTES, ct, CIPHERTEXTBYTES);

  /* Overwrite pre-k with z on re-encryption failure */
  Cmov(kr, sk+SECRETKEYBYTES-SYMBYTES, SYMBYTES, fail);

  /* hash concatenation of pre-k and H(c) to k */
  Shake256(ss, 32, kr, 2*SYMBYTES);
  return 0;
}


// Generates public and private key for the CPA-secure
// public-key encryption scheme underlying Kyber
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::IndCpaKeypair(byte *pk, byte *sk) {
  word32 i;
  byte buf[2*SYMBYTES];
  const byte *publicseed = buf;
  const byte *noiseseed = buf+SYMBYTES;
  byte nonce = 0;


  polyvec a[K];

  polyvec e, pkpv, skpv;

  RandomBytes(buf, SYMBYTES);
  Sha3_512(buf, buf, SYMBYTES);

  GenMatrix(a, publicseed, 0);

  for(i=0;i<K;i++) {
    PolyGetNoise(&skpv.vec[i], noiseseed, nonce++, ETA1);
  }
      
  for(i=0;i<K;i++) {
    PolyGetNoise(&e.vec[i], noiseseed, nonce++, ETA1);
  }

  PolyvecNtt(&skpv);
  PolyvecNtt(&e);

  // matrix-vector multiplication
  for(i=0;i<K;i++) {
      PolyvecBasemulAccMontgomery(&pkpv.vec[i], &a[i], &skpv);
      PolyToMont(&pkpv.vec[i]);
  }

  PolyvecAdd(&pkpv, &pkpv, &e);
  PolyvecReduce(&pkpv);

  PackSk(sk, &skpv);
  PackPk(pk, &pkpv, publicseed);


};

// Encryption function of the CPA-secure
// public-key encryption scheme underlying Kyber.
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::IndCpaEnc(byte *c, const byte *m, const byte *pk, const byte *coins)
{
  word32 i;
  byte seed[SYMBYTES];
  byte nonce = 0;
  polyvec sp, pkpv, ep, at[K], bp;
  poly v, k, epp;

  UnpackPk(&pkpv, seed, pk);
  PolyFromMsg(&k, m);
  GenMatrix(at, seed, 1);

  for(i=0;i<K;i++) {
    PolyGetNoise(sp.vec+i, coins, nonce++, ETA1);
  }
  for(i=0;i<K;i++) {
    PolyGetNoise(ep.vec+i, coins, nonce++, ETA2);
  }
    
  PolyGetNoise(&epp, coins, nonce++, ETA2);
  PolyvecNtt(&sp);
  // matrix-vector multiplication
  for(i=0;i<K;i++) {
    PolyvecBasemulAccMontgomery(&bp.vec[i], &at[i], &sp);
  }
    

  PolyvecBasemulAccMontgomery(&v, &pkpv, &sp);

  PolyvecInvNttToMont(&bp);
  PolyInvNttToMont(&v);

  PolyvecAdd(&bp, &bp, &ep);
  PolyAdd(&v, &v, &epp);
  PolyAdd(&v, &v, &k);
  PolyvecReduce(&bp);
  PolyReduce(&v);

  PackCiphertext(c, &bp, &v);
}

// Decryption function of the CPA-secure
// public-key encryption scheme underlying Kyber.
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::IndCpaDec(byte *m, const byte *c, const byte *sk)
{
  polyvec bp, skpv;
  poly v, mp;

  UnpackCiphertext(&bp, &v, c);
  UnpackSk(&skpv, sk);

  PolyvecNtt(&bp);
  PolyvecBasemulAccMontgomery(&mp, &skpv, &bp);
  PolyInvNttToMont(&mp);

  PolySub(&mp, &v, &mp);
  PolyReduce(&mp);

  PolyToMsg(m, &mp);
}




/* Run rejection sampling on uniform random bytes to generate
* uniform random integers mod q (rej_uniform function of 
* example implementation)*/
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
word32 Kyber<T_K, T_C, T_C2, T_ETA1>::RejUniform(sword16 *r, word32 len, const byte *buf, word32 buflen) {
  word32 ctr, pos;
  word16 val0, val1;

  ctr = pos = 0;
  while(ctr < len && pos + 3 <= buflen) {
    val0 = ((buf[pos+0] >> 0) | ((word16)buf[pos+1] << 8)) & 0xFFF;
    val1 = ((buf[pos+1] >> 4) | ((word16)buf[pos+2] << 4)) & 0xFFF;
    pos += 3;

    if(val0 < Q) {
      r[ctr++] = val0;
    }
    if (ctr < len && val1 < Q) {
      r[ctr++] = val1;
    }
  }
  return ctr;
};


//Deterministically generate matrix A (or the transpose of A)
//from a seed. Entries of the matrix are polynomials that look
//uniformly random. Performs rejection sampling on output of
//a XOF
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::GenMatrix(polyvec *a, const byte *seed, sword32 transposed)
{
  word32 ctr, i, j, k;
  word32 bufLen, off;
  const sword32 xofBlockBytes = 168;
  const size_t genMatrixNBlocks = ((12*N/8*(1 << 12)/Q + xofBlockBytes)/xofBlockBytes);
  byte buf[genMatrixNBlocks*xofBlockBytes+2];

  Shake128Init();
  for(i=0;i<K;i++) {
    for(j=0;j<K;j++) {
      if(transposed) {
        KyberShake128Absorb(seed, i, j);
      }

      else {

        KyberShake128Absorb(seed, j, i);
      }
      
      Shake128SqueezeBlocks(buf, genMatrixNBlocks);
      bufLen = genMatrixNBlocks*xofBlockBytes;
      ctr = RejUniform(a[i].vec[j].coeffs, N, buf, bufLen);
     
      while(ctr < N) {
        off = bufLen % 3;
        for(k = 0; k < off; k++)
          buf[k] = buf[bufLen - off + k];
        Shake128SqueezeBlocks(buf + off, 1);
        bufLen = off + xofBlockBytes;
        ctr += RejUniform(a[i].vec[j].coeffs + ctr, N - ctr, buf, bufLen);
      }
    }
  }
};


// Serialize the public key as concatenation of the
// serialized vector of polynomials pk
// and the public seed used to generate the matrix A.
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PackPk(byte *r, polyvec *pk, const byte *seed)
{
  size_t i;
  PolyvecToBytes(r, pk);
  for(i=0;i<SYMBYTES;i++)
    r[i+POLYVECBYTES] = seed[i];
}

// De-serialize public key from a byte array;
// approximate inverse of pack_pk
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::UnpackPk(polyvec *pk, byte *seed, const byte *packedpk)
{
  size_t i;
  PolyvecFromBytes(pk, packedpk);
  for(i=0;i<SYMBYTES;i++)
    seed[i] = packedpk[i+POLYVECBYTES];
}


//Serialize the secret key
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PackSk(byte *r, polyvec *sk)
{
  PolyvecToBytes(r, sk);
}


//De-serialize the secret key
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::UnpackSk(polyvec *sk, const byte *packedsk)
{
  PolyvecFromBytes(sk, packedsk);
}

// Serialize the ciphertext as concatenation of the
// compressed and serialized vector of polynomials b
// and the compressed and serialized polynomial v
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PackCiphertext(byte *r, polyvec *b, poly *v)
{
  PolyvecCompress(r, b);
  PolyCompress(r+POLYVECCOMPRESSEDBYTES, v);
}

//De-serialize and decompress ciphertext from a byte array;
//approximate inverse of PackCiphertext
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::UnpackCiphertext(polyvec *b, poly *v, const byte *c)
{
  PolyvecDecompress(b, c);
  PolyDecompress(v, c+POLYVECCOMPRESSEDBYTES);
}

template class Kyber<2, 320, 128, 3>;
template class Kyber<3, 320, 128, 2>;
template class Kyber<4, 352, 160, 2>;
 


NAMESPACE_END