/* The C++ file for CRYSTALS-Kyber. Adapted from the reference
implementation of kyber by the CRYSTALS team 
(https://github.com/pq-crystals/kyber) */

#include <iterator>
#include <iostream>

#include "kyber.h"

#include "sha3.h"
#include "shake.h"
#include "misc.h"

NAMESPACE_BEGIN(CryptoPP)

// Generates public and private key
// for CCA-secure Kyber key encapsulation mechanism
int Kyber::KemKeypair(unsigned char *pk, unsigned char *sk)
{
  size_t i;
  IndCpaKeypair(pk, sk);
  for(i=0;i<KYBER_INDCPA_PUBLICKEYBYTES;i++) {
    sk[i+KYBER_INDCPA_SECRETKEYBYTES] = pk[i];
  }
  SHA3_256 sha3_256 = SHA3_256();
  sha3_256.Update(pk, KYBER_PUBLICKEYBYTES);
  sha3_256.Final(sk+KYBER_SECRETKEYBYTES-2*KYBER_SYMBYTES);
  /* Value z for pseudo-random output on reject */
 
  randombytes(sk+KYBER_SECRETKEYBYTES-KYBER_SYMBYTES, KYBER_SYMBYTES);
  return 0;
}

// Generates cipher text and shared
// secret for given public key

int Kyber::KemEnc(unsigned char *ct,
                   unsigned char *ss,
                   const unsigned char *pk)
{
  uint8_t buf[2*KYBER_SYMBYTES];
  /* Will contain key, coins */
  uint8_t kr[2*KYBER_SYMBYTES];
  
  randombytes(buf, KYBER_SYMBYTES);
  /* Don't release system RNG output */
  SHA3_256 sha3 = SHA3_256();
  sha3.Update(buf, KYBER_SYMBYTES);
  sha3.Final(buf);

  sha3.Restart();

  /* Multitarget countermeasure for coins + contributory KEM */
  sha3.Update(pk, KYBER_PUBLICKEYBYTES);
  sha3.Final(buf+KYBER_SYMBYTES);

  sha3.Restart();

  SHA3_512 sha3_512 = SHA3_512();
  sha3_512.Update(buf, 2*KYBER_SYMBYTES);
  sha3_512.Final(kr);

  /* coins are in kr+KYBER_SYMBYTES */
  IndCpaEnc(ct, buf, pk, kr+KYBER_SYMBYTES);

  /* overwrite coins in kr with H(c) */
  sha3.Update(ct, KYBER_CIPHERTEXTBYTES);
  sha3.Final(kr+KYBER_SYMBYTES);
  /* hash concatenation of pre-k and H(c) to k */
  SHAKE256 shake = SHAKE256(32);
  shake.Update(kr, 2*KYBER_SYMBYTES);
  shake.Final(ss);
  return 0;
}


// Generates shared secret for given
// cipher text and private key
// On failure, ss will contain a pseudo-random value.
int Kyber::KemDec(unsigned char *ss,
                   const unsigned char *ct,
                   const unsigned char *sk)
{
  size_t i;
  int fail;
  uint8_t buf[2*KYBER_SYMBYTES];
  /* Will contain key, coins */
  uint8_t kr[2*KYBER_SYMBYTES];
  uint8_t cmp[KYBER_CIPHERTEXTBYTES];
  const uint8_t *pk = sk+KYBER_INDCPA_SECRETKEYBYTES;

  IndCpaDec(buf, ct, sk);

  /* Multitarget countermeasure for coins + contributory KEM */
  for(i=0;i<KYBER_SYMBYTES;i++) {
    buf[KYBER_SYMBYTES+i] = sk[KYBER_SECRETKEYBYTES-2*KYBER_SYMBYTES+i];
  }
  SHA3_512 sha3_512 = SHA3_512();
  sha3_512.Update(buf, 2*KYBER_SYMBYTES);
  sha3_512.Final(kr);

  /* coins are in kr+KYBER_SYMBYTES */
  IndCpaEnc(cmp, buf, pk, kr+KYBER_SYMBYTES);
  fail = !(VerifyBufsEqual(ct, cmp, KYBER_CIPHERTEXTBYTES));

  /* overwrite coins in kr with H(c) */
  SHA3_256 sha3_256 = SHA3_256();
  sha3_256.Update(ct, KYBER_CIPHERTEXTBYTES);
  sha3_256.Final(kr+KYBER_SYMBYTES);

  /* Overwrite pre-k with z on re-encryption failure */
  Cmov(kr, sk+KYBER_SECRETKEYBYTES-KYBER_SYMBYTES, KYBER_SYMBYTES, fail);

  /* hash concatenation of pre-k and H(c) to k */
  SHAKE256 shake = SHAKE256(32);
  shake.Update(kr, 2*KYBER_SYMBYTES);
  shake.Final(ss);
  return 0;
}


// Generates public and private key for the CPA-secure
//spublic-key encryption scheme underlying Kyber
void Kyber::IndCpaKeypair(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES], uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES]) {
  unsigned int i;
  uint8_t buf[2*KYBER_SYMBYTES];
  const uint8_t *publicseed = buf;
  const uint8_t *noiseseed = buf+KYBER_SYMBYTES;
  uint8_t nonce = 0;


  polyvec a[KYBER_K];

  //   poly a[kyber_k][kyber_k];
  polyvec e, pkpv, skpv;

  randombytes(buf, KYBER_SYMBYTES);

  SHA3_512 sha3_512 = SHA3_512();

  sha3_512.Update(buf, KYBER_SYMBYTES);
  sha3_512.Final(buf);


  GenMatrix(a, publicseed, 0);

  for(i=0;i<KYBER_K;i++) {
    PolyGetNoise(&skpv.vec[i], noiseseed, nonce++);
  }
      
  for(i=0;i<KYBER_K;i++) {
    PolyGetNoise(&e.vec[i], noiseseed, nonce++);
  }

  PolyvecNtt(&skpv);
  PolyvecNtt(&e);

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++) {
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
void Kyber::IndCpaEnc(uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                const uint8_t coins[KYBER_SYMBYTES])
{
  unsigned int i;
  uint8_t seed[KYBER_SYMBYTES];
  uint8_t nonce = 0;
  polyvec sp, pkpv, ep, at[KYBER_K], bp;
  poly v, k, epp;

  UnpackPk(&pkpv, seed, pk);
  PolyFromMsg(&k, m);
  GenMatrix(at, seed, 1);

  for(i=0;i<KYBER_K;i++) {
    PolyGetNoise(sp.vec+i, coins, nonce++);
  }
  for(i=0;i<KYBER_K;i++) {
    PolyGetNoise(ep.vec+i, coins, nonce++);
  }
    
  PolyGetNoise(&epp, coins, nonce++);

  PolyvecNtt(&sp);

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++) {
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

void Kyber::IndCpaDec(uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES])
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
unsigned int Kyber::RejUniform(int16_t *r, unsigned int len, const uint8_t *buf, unsigned int buflen) {
  unsigned int ctr, pos;
  uint16_t val;

  ctr = pos = 0;
  while(ctr < len && pos + 2 <= buflen) {
    val = buf[pos] | ((uint16_t)buf[pos+1] << 8);
    pos += 2;

    if(val < 19*KYBER_Q) {
      val -= (val >> 12)*KYBER_Q; // Barrett reduction
      r[ctr++] = (int16_t)val;
    }
  }

  return ctr;

};


//Deterministically generate matrix A (or the transpose of A)
//from a seed. Entries of the matrix are polynomials that look
//uniformly random. Performs rejection sampling on output of
//a XOF
void Kyber::GenMatrix(polyvec *a, const uint8_t seed[KYBER_SYMBYTES], int transposed)
{
  unsigned int ctr, i, j;
  int xof_blockbytes = 168;
  size_t gen_matrix_nblocks = ((2*KYBER_N*(1U << 16)/(19*KYBER_Q) + xof_blockbytes)/xof_blockbytes);
  uint8_t buf[gen_matrix_nblocks*xof_blockbytes];
  keccak_state state;


  for(i=0;i<KYBER_K;i++) {
    for(j=0;j<KYBER_K;j++) {
      if(transposed) {
        XOFAbsorb(&state, seed, i, j);
      }

      else {

        XOFAbsorb(&state, seed, j, i);
      }
      
      XOFSqueezeBlocks(buf, gen_matrix_nblocks, &state);

      ctr = RejUniform(a[i].vec[j].coeffs, KYBER_N, buf, sizeof(buf));
      
      while(ctr < KYBER_N) {
        XOFSqueezeBlocks(buf, 1, &state);
        ctr += RejUniform(a[i].vec[j].coeffs + ctr, KYBER_N - ctr, buf,
                           xof_blockbytes);
      }
    }
  }
};


// Serialize the public key as concatenation of the
// serialized vector of polynomials pk
// and the public seed used to generate the matrix A.
void Kyber::PackPk(uint8_t r[KYBER_INDCPA_PUBLICKEYBYTES],
                    polyvec *pk,
                    const uint8_t seed[KYBER_SYMBYTES])
{
  size_t i;
  PolyvecToBytes(r, pk);
  for(i=0;i<KYBER_SYMBYTES;i++)
    r[i+KYBER_POLYVECBYTES] = seed[i];
}

// De-serialize public key from a byte array;
// approximate inverse of pack_pk

void Kyber::UnpackPk(polyvec *pk,
                      uint8_t seed[KYBER_SYMBYTES],
                      const uint8_t packedpk[KYBER_INDCPA_PUBLICKEYBYTES])
{
  size_t i;
  PolyvecFromBytes(pk, packedpk);
  for(i=0;i<KYBER_SYMBYTES;i++)
    seed[i] = packedpk[i+KYBER_POLYVECBYTES];
}


//Serialize the secret key

void Kyber::PackSk(uint8_t r[KYBER_INDCPA_SECRETKEYBYTES], polyvec *sk)
{
  PolyvecToBytes(r, sk);
}


//De-serialize the secret key
void Kyber::UnpackSk(polyvec *sk,
                      const uint8_t packedsk[KYBER_INDCPA_SECRETKEYBYTES])
{
  PolyvecFromBytes(sk, packedsk);
}

// Serialize the ciphertext as concatenation of the
// compressed and serialized vector of polynomials b
// and the compressed and serialized polynomial v
void Kyber::PackCiphertext(uint8_t r[KYBER_INDCPA_BYTES],
                            polyvec *b,
                            poly *v)
{
  PolyvecCompress(r, b);
  PolyCompress(r+KYBER_POLYVECCOMPRESSEDBYTES, v);
}

//De-serialize and decompress ciphertext from a byte array;
//approximate inverse of PackCiphertext
void Kyber::UnpackCiphertext(polyvec *b,
                              poly *v,
                              const uint8_t c[KYBER_INDCPA_BYTES])
{
  PolyvecDecompress(b, c);
  PolyDecompress(v, c+KYBER_POLYVECCOMPRESSEDBYTES);
}





NAMESPACE_END