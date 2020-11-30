/* The C++ file for CRYSTALS-Kyber. Adapted from the reference
implementation of kyber by the CRYSTALS team 
(https://github.com/pq-crystals/kyber) */

#include <iterator>

#include "kyber.h"
#include "kyber_utils.h"

#include "sha3.h"
#include "shake.h"

NAMESPACE_BEGIN(CryptoPP)

/*************************************************
* Name:        crypto_kem_keypair
*
* Description: Generates public and private key
*              for CCA-secure Kyber key encapsulation mechanism
*
* Arguments:   - unsigned char *pk: pointer to output public key
*                (an already allocated array of CRYPTO_PUBLICKEYBYTES bytes)
*              - unsigned char *sk: pointer to output private key
*                (an already allocated array of CRYPTO_SECRETKEYBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int Kyber::CryptoKemKeypair(unsigned char *pk, unsigned char *sk)
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
  RandomBytes rnd = RandomBytes();
  rnd.randombytes(sk+KYBER_SECRETKEYBYTES-KYBER_SYMBYTES, KYBER_SYMBYTES);
  return 0;
}


void Kyber::IndCpaKeypair(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES], uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES]) {
  unsigned int i;
  uint8_t buf[2*KYBER_SYMBYTES];
  const uint8_t *publicseed = buf;
  const uint8_t *noiseseed = buf+KYBER_SYMBYTES;
  uint8_t nonce = 0;


  polyvec a[KYBER_K];

  //   poly a[kyber_k][kyber_k];
  polyvec e, pkpv, skpv;

  RandomBytes rnd = RandomBytes();

  rnd.randombytes(buf, KYBER_SYMBYTES);

  SHA3_512 sha3_512 = SHA3_512();

  sha3_512.Update(buf, KYBER_SYMBYTES);
  sha3_512.Final(buf);

  // a = AllocMatrix();

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



/* Run rejection sampling on uniform random bytes to generate
* uniform random integers mod q (rej_uniform function of 
* example implementation)*/
unsigned int Kyber::RejUniform(int16_t *r, unsigned int len, const uint8_t *buf, unsigned int buflen) {
  unsigned int ctr,pos;
  uint16_t val;

  ctr = pos = 0;
  while(ctr < KYBER_N && pos + 2 <= sizeof(buf)) {
    val = buf[pos] | ((uint16_t)buf[pos+1] << 8);
    pos += 2;

    if(val < 19*KYBER_Q) {
      val -= (val >> 12)*KYBER_Q; // Barrett reduction
      r[ctr++] = (int16_t)val;
    }
  }

  return ctr;

};


void Kyber::GenMatrix(polyvec *a, const uint8_t seed[], int transposed)
{
  unsigned int ctr, i, j;
  int xof_blockbytes = 168;
  int gen_matrix_nblocks = ((2*KYBER_N*(1U << 16)/(19*KYBER_Q) + xof_blockbytes)/xof_blockbytes);
  uint8_t buf[gen_matrix_nblocks*xof_blockbytes];
  keccak_state state;
  SHAKE128 xof = SHAKE128();

  uint8_t extSeed[KYBER_SYMBYTES + 2];
  std::copy(seed, seed + KYBER_SYMBYTES, extSeed);



  for(i=0;i<KYBER_K;i++) {
    for(j=0;j<KYBER_K;j++) {
      if(transposed) {
        extSeed[KYBER_SYMBYTES] = i;
        extSeed[KYBER_SYMBYTES + 1] = j;
        Shake128Absorb(&state, extSeed, sizeof(extSeed));
      }

       // xof_absorb(&state, seed, i, j);
      else {
        extSeed[KYBER_SYMBYTES] = j;
        extSeed[KYBER_SYMBYTES + 1] = i;
        Shake128Absorb(&state, extSeed, sizeof(extSeed));
      }
      
      //TODO: KeccakF1600 oma toteutus
      Shake128SqueezeBlocks(buf, gen_matrix_nblocks, &state);


      ctr = RejUniform(a[i].vec[j].coeffs, KYBER_N, buf, sizeof(buf));
      

      while(ctr < KYBER_N) {
        Shake128SqueezeBlocks(buf, 1, &state);
        ctr += RejUniform(a[i].vec[j].coeffs + ctr, KYBER_N - ctr, buf,
                           xof_blockbytes);
      }
    }
  }
};


/*************************************************
* Name:        pack_pk
*
* Description: Serialize the public key as concatenation of the
*              serialized vector of polynomials pk
*              and the public seed used to generate the matrix A.
*
* Arguments:   uint8_t *r:          pointer to the output serialized public key
*              polyvec *pk:         pointer to the input public-key polyvec
*              const uint8_t *seed: pointer to the input public seed
**************************************************/
void Kyber::PackPk(uint8_t r[KYBER_INDCPA_PUBLICKEYBYTES],
                    polyvec *pk,
                    const uint8_t seed[KYBER_SYMBYTES])
{
  size_t i;
  PolyvecToBytes(r, pk);
  for(i=0;i<KYBER_SYMBYTES;i++)
    r[i+KYBER_POLYVECBYTES] = seed[i];
}

/*************************************************
* Name:        unpack_pk
*
* Description: De-serialize public key from a byte array;
*              approximate inverse of pack_pk
*
* Arguments:   - polyvec *pk:             pointer to output public-key
*                                         polynomial vector
*              - uint8_t *seed:           pointer to output seed to generate
*                                         matrix A
*              - const uint8_t *packedpk: pointer to input serialized public key
**************************************************/
void Kyber::UnpackPk(polyvec *pk,
                      uint8_t seed[KYBER_SYMBYTES],
                      const uint8_t packedpk[KYBER_INDCPA_PUBLICKEYBYTES])
{
  size_t i;
  PolyvecFromBytes(pk, packedpk);
  for(i=0;i<KYBER_SYMBYTES;i++)
    seed[i] = packedpk[i+KYBER_POLYVECBYTES];
}

/*************************************************
* Name:        pack_sk
*
* Description: Serialize the secret key
*
* Arguments:   - uint8_t *r:  pointer to output serialized secret key
*              - polyvec *sk: pointer to input vector of polynomials (secret key)
**************************************************/
void Kyber::PackSk(uint8_t r[KYBER_INDCPA_SECRETKEYBYTES], polyvec *sk)
{
  PolyvecToBytes(r, sk);
}

/*************************************************
* Name:        unpack_sk
*
* Description: De-serialize the secret key;
*              inverse of pack_sk
*
* Arguments:   - polyvec *sk:             pointer to output vector of
*                                         polynomials (secret key)
*              - const uint8_t *packedsk: pointer to input serialized secret key
**************************************************/
void Kyber::UnpackSk(polyvec *sk,
                      const uint8_t packedsk[KYBER_INDCPA_SECRETKEYBYTES])
{
  PolyvecFromBytes(sk, packedsk);
}

/*************************************************
* Name:        pack_ciphertext
*
* Description: Serialize the ciphertext as concatenation of the
*              compressed and serialized vector of polynomials b
*              and the compressed and serialized polynomial v
*
* Arguments:   uint8_t *r: pointer to the output serialized ciphertext
*              poly *pk:   pointer to the input vector of polynomials b
*              poly *v:    pointer to the input polynomial v
**************************************************/
void Kyber::PackCiphertext(uint8_t r[KYBER_INDCPA_BYTES],
                            polyvec *b,
                            poly *v)
{
  PolyvecCompress(r, b);
  PolyCompress(r+KYBER_POLYVECCOMPRESSEDBYTES, v);
}

/*************************************************
* Name:        unpack_ciphertext
*
* Description: De-serialize and decompress ciphertext from a byte array;
*              approximate inverse of pack_ciphertext
*
* Arguments:   - polyvec *b:       pointer to the output vector of polynomials b
*              - poly *v:          pointer to the output polynomial v
*              - const uint8_t *c: pointer to the input serialized ciphertext
**************************************************/
void Kyber::UnpackCiphertext(polyvec *b,
                              poly *v,
                              const uint8_t c[KYBER_INDCPA_BYTES])
{
  PolyvecDecompress(b, c);
  PolyDecompress(v, c+KYBER_POLYVECCOMPRESSEDBYTES);
}





NAMESPACE_END