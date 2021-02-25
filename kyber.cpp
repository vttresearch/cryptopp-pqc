/* The C++ file for CRYSTALS-Kyber. Adapted from the reference
implementation of kyber by the CRYSTALS team 
(https://github.com/pq-crystals/kyber) */

#include "pch.h"
#include "kyber.h"

#include "sha3.h"
#include "shake.h"
#include "misc.h"

NAMESPACE_BEGIN(CryptoPP)

// Generates public and private key
// for CCA-secure Kyber key encapsulation mechanism
template<int T_K, unsigned int T_Compr>
int Kyber<T_K, T_Compr>::KemKeypair(unsigned char *pk, unsigned char *sk)
{
  size_t i; 
  IndCpaKeypair(pk, sk);
  for(i=0;i<INDCPA_PUBLICKEYBYTES;i++) {
    sk[i+INDCPA_SECRETKEYBYTES] = pk[i];
  }
  SHA3_256 sha3_256 = SHA3_256();
  sha3_256.Update(pk, PUBLICKEYBYTES);
  sha3_256.Final(sk+SECRETKEYBYTES-2*SYMBYTES);
  /* Value z for pseudo-random output on reject */
 
  randombytes(sk+SECRETKEYBYTES-SYMBYTES, SYMBYTES);
  return 0;
}

// Generates cipher text and shared
// secret for given public key
template<int T_K, unsigned int T_Compr>
int Kyber<T_K,T_Compr>::KemEnc(unsigned char *ct, unsigned char *ss, const unsigned char *pk)
{
  uint8_t buf[2*SYMBYTES];
  /* Will contain key, coins */
  uint8_t kr[2*SYMBYTES];
  
  randombytes(buf, SYMBYTES);
  /* Don't release system RNG output */
  SHA3_256 sha3 = SHA3_256();
  sha3.Update(buf, SYMBYTES);
  sha3.Final(buf);

  sha3.Restart();

  /* Multitarget countermeasure for coins + contributory KEM */
  sha3.Update(pk, PUBLICKEYBYTES);
  sha3.Final(buf+ SYMBYTES);

  sha3.Restart();

  SHA3_512 sha3_512 = SHA3_512();
  sha3_512.Update(buf, 2* SYMBYTES);
  sha3_512.Final(kr);

  /* coins are in kr+KYBER_SYMBYTES */
  IndCpaEnc(ct, buf, pk, kr+ SYMBYTES);

  /* overwrite coins in kr with H(c) */
  sha3.Update(ct, CIPHERTEXTBYTES);
  sha3.Final(kr+SYMBYTES);
  /* hash concatenation of pre-k and H(c) to k */
  SHAKE256 shake = SHAKE256(32);
  shake.Update(kr, 2*SYMBYTES);
  shake.Final(ss);
  return 0;
}


// Generates shared secret for given
// cipher text and private key
// On failure, ss will contain a pseudo-random value.
template<int T_K, unsigned int T_Compr>
int Kyber<T_K,T_Compr>::KemDec(unsigned char *ss, const unsigned char *ct, const unsigned char *sk)
{
  size_t i;
  int fail;
  uint8_t buf[2*SYMBYTES];
  /* Will contain key, coins */
  uint8_t kr[2*SYMBYTES];
  uint8_t cmp[CIPHERTEXTBYTES];
  const uint8_t *pk = sk+INDCPA_SECRETKEYBYTES;

  IndCpaDec(buf, ct, sk);

  /* Multitarget countermeasure for coins + contributory KEM */
  for(i=0;i<SYMBYTES;i++) {
    buf[SYMBYTES+i] = sk[SECRETKEYBYTES-2*SYMBYTES+i];
  }
  SHA3_512 sha3_512 = SHA3_512();
  sha3_512.Update(buf, 2*SYMBYTES);
  sha3_512.Final(kr);

  /* coins are in kr+KYBER_SYMBYTES */
  IndCpaEnc(cmp, buf, pk, kr+SYMBYTES);
  fail = !(VerifyBufsEqual(ct, cmp, CIPHERTEXTBYTES));

  /* overwrite coins in kr with H(c) */
  SHA3_256 sha3_256 = SHA3_256();
  sha3_256.Update(ct, CIPHERTEXTBYTES);
  sha3_256.Final(kr+SYMBYTES);

  /* Overwrite pre-k with z on re-encryption failure */
  Cmov(kr, sk+SECRETKEYBYTES-SYMBYTES, SYMBYTES, fail);

  /* hash concatenation of pre-k and H(c) to k */
  SHAKE256 shake = SHAKE256(32);
  shake.Update(kr, 2*SYMBYTES);
  shake.Final(ss);
  return 0;
}


// Generates public and private key for the CPA-secure
// public-key encryption scheme underlying Kyber
template<int T_K, unsigned int T_Compr>
void Kyber<T_K,T_Compr>::IndCpaKeypair(uint8_t *pk, uint8_t *sk) {
  unsigned int i;
  uint8_t buf[2*SYMBYTES];
  const uint8_t *publicseed = buf;
  const uint8_t *noiseseed = buf+SYMBYTES;
  uint8_t nonce = 0;


  polyvec a[K];

  //   poly a[kyber_k][kyber_k];
  polyvec e, pkpv, skpv;

  randombytes(buf, SYMBYTES);

  SHA3_512 sha3_512 = SHA3_512();

  sha3_512.Update(buf, SYMBYTES);
  sha3_512.Final(buf);


  GenMatrix(a, publicseed, 0);

  for(i=0;i<K;i++) {
    PolyGetNoise(&skpv.vec[i], noiseseed, nonce++);
  }
      
  for(i=0;i<K;i++) {
    PolyGetNoise(&e.vec[i], noiseseed, nonce++);
  }

  PolyvecNtt(&skpv);
  PolyvecNtt(&e);

  // matrix-vector multiplication
  for(i=0;i<K;i++) {
      PolyvecPointwiseAccMontgomery(&pkpv.vec[i], &a[i], &skpv);
      PolyToMont(&pkpv.vec[i]);
  }

  PolyvecAdd(&pkpv, &pkpv, &e);
  PolyvecReduce(&pkpv);

  PackSk(sk, &skpv);
  PackPk(pk, &pkpv, publicseed);


};

// Encryption function of the CPA-secure
// public-key encryption scheme underlying Kyber.
template<int T_K, unsigned int T_Compr>
void Kyber<T_K,T_Compr>::IndCpaEnc(uint8_t *c, const uint8_t *m, const uint8_t *pk, const uint8_t *coins)
{
  unsigned int i;
  uint8_t seed[SYMBYTES];
  uint8_t nonce = 0;
  polyvec sp, pkpv, ep, at[K], bp;
  poly v, k, epp;

  UnpackPk(&pkpv, seed, pk);
  PolyFromMsg(&k, m);
  GenMatrix(at, seed, 1);

  for(i=0;i<K;i++) {
    PolyGetNoise(sp.vec+i, coins, nonce++);
  }
  for(i=0;i<K;i++) {
    PolyGetNoise(ep.vec+i, coins, nonce++);
  }
    
  PolyGetNoise(&epp, coins, nonce++);

  PolyvecNtt(&sp);

  // matrix-vector multiplication
  for(i=0;i<K;i++) {
    PolyvecPointwiseAccMontgomery(&bp.vec[i], &at[i], &sp);
  }
    

  PolyvecPointwiseAccMontgomery(&v, &pkpv, &sp);

  PolyvecInvNttToMont(&bp);
  PolyInvNttToMont(&v);

  PolyvecAdd(&bp, &bp, &ep);
  PolyAdd(&v, &v, &epp);
  PolyAdd(&v, &v, &k);
  PolyvecReduce(&bp);
  PolyReduce(&v);

  PackCiphertext(c, &bp, &v);
}

// ecryption function of the CPA-secure
// public-key encryption scheme underlying Kyber.

template<int T_K, unsigned int T_Compr>
void Kyber<T_K,T_Compr>::IndCpaDec(uint8_t *m, const uint8_t *c, const uint8_t *sk)
{
  polyvec bp, skpv;
  poly v, mp;

  UnpackCiphertext(&bp, &v, c);
  UnpackSk(&skpv, sk);

  PolyvecNtt(&bp);
  PolyvecPointwiseAccMontgomery(&mp, &skpv, &bp);
  PolyInvNttToMont(&mp);

  PolySub(&mp, &v, &mp);
  PolyReduce(&mp);

  PolyToMsg(m, &mp);
}




/* Run rejection sampling on uniform random bytes to generate
* uniform random integers mod q (rej_uniform function of 
* example implementation)*/
template<int T_K, unsigned int T_Compr>
unsigned int Kyber<T_K,T_Compr>::RejUniform(int16_t *r, unsigned int len, const uint8_t *buf, unsigned int buflen) {
  unsigned int ctr, pos;
  uint16_t val;

  ctr = pos = 0;
  while(ctr < len && pos + 2 <= buflen) {
    val = buf[pos] | ((uint16_t)buf[pos+1] << 8);
    pos += 2;

    if(val < 19*Q) {
      val -= (val >> 12)*Q; // Barrett reduction
      r[ctr++] = (int16_t)val;
    }
  }

  return ctr;

};


//Deterministically generate matrix A (or the transpose of A)
//from a seed. Entries of the matrix are polynomials that look
//uniformly random. Performs rejection sampling on output of
//a XOF
template<int T_K, unsigned int T_Compr>
void Kyber<T_K,T_Compr>::GenMatrix(polyvec *a, const uint8_t *seed, int transposed)
{
  unsigned int ctr, i, j;
  const int xof_blockbytes = 168;
  const size_t gen_matrix_nblocks = ((2*N*(1U << 16)/(19*Q) + xof_blockbytes)/xof_blockbytes);
  uint8_t buf[gen_matrix_nblocks*xof_blockbytes];
  keccak_state state;


  for(i=0;i<K;i++) {
    for(j=0;j<K;j++) {
      if(transposed) {
        XOFAbsorb(&state, seed, i, j);
      }

      else {

        XOFAbsorb(&state, seed, j, i);
      }
      
      XOFSqueezeBlocks(buf, gen_matrix_nblocks, &state);

      ctr = RejUniform(a[i].vec[j].coeffs, N, buf, sizeof(buf));
      
      while(ctr < N) {
        XOFSqueezeBlocks(buf, 1, &state);
        ctr += RejUniform(a[i].vec[j].coeffs + ctr, N - ctr, buf,
                           xof_blockbytes);
      }
    }
  }
};


// Serialize the public key as concatenation of the
// serialized vector of polynomials pk
// and the public seed used to generate the matrix A.
template<int T_K, unsigned int T_Compr>
void Kyber<T_K,T_Compr>::PackPk(uint8_t *r, polyvec *pk, const uint8_t *seed)
{
  size_t i;
  PolyvecToBytes(r, pk);
  for(i=0;i<SYMBYTES;i++)
    r[i+POLYVECBYTES] = seed[i];
}

// De-serialize public key from a byte array;
// approximate inverse of pack_pk

template<int T_K, unsigned int T_Compr>
void Kyber<T_K,T_Compr>::UnpackPk(polyvec *pk, uint8_t *seed, const uint8_t *packedpk)
{
  size_t i;
  PolyvecFromBytes(pk, packedpk);
  for(i=0;i<SYMBYTES;i++)
    seed[i] = packedpk[i+POLYVECBYTES];
}


//Serialize the secret key

template<int T_K, unsigned int T_Compr>
void Kyber<T_K,T_Compr>::PackSk(uint8_t *r, polyvec *sk)
{
  PolyvecToBytes(r, sk);
}


//De-serialize the secret key
template<int T_K, unsigned int T_Compr>
void Kyber<T_K,T_Compr>::UnpackSk(polyvec *sk, const uint8_t *packedsk)
{
  PolyvecFromBytes(sk, packedsk);
}

// Serialize the ciphertext as concatenation of the
// compressed and serialized vector of polynomials b
// and the compressed and serialized polynomial v
template<int T_K, unsigned int T_Compr>
void Kyber<T_K,T_Compr>::PackCiphertext(uint8_t *r, polyvec *b, poly *v)
{
  PolyvecCompress(r, b);
  PolyCompress(r+POLYVECCOMPRESSEDBYTES, v);
}

//De-serialize and decompress ciphertext from a byte array;
//approximate inverse of PackCiphertext
template<int T_K, unsigned int T_Compr>
void Kyber<T_K,T_Compr>::UnpackCiphertext(polyvec *b, poly *v, const uint8_t *c)
{
  PolyvecDecompress(b, c);
  PolyDecompress(v, c+POLYVECCOMPRESSEDBYTES);
}

template class Kyber<2, 320>;
template class Kyber<3, 320>;
template class Kyber<4, 352>;
 


NAMESPACE_END