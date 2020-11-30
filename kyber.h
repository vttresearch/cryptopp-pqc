/* The header file for CRYSTALS-Kyber. Adapted from the reference
implementation of kyber by the CRYSTALS team 
(https://github.com/pq-crystals/kyber) */

#ifndef CRYPTOPP_KYBER_H
#define CRYPTOPP_KYBER_H

#include "cryptlib.h"
#include "kyber_params.h"

NAMESPACE_BEGIN(CryptoPP)




class Kyber : public Algorithm
{
public:
    int CryptoKemKeypair(unsigned char *pk, unsigned char *sk);

    

protected:




    
    Kyber(word16 kyber_version) {};


    
    //Return the name of the algorithm
    std::string AlgorithmName() const {return "CRYSTALS-Kyber-768";}

    /*
    * Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
    * coeffs[0] + X*coeffs[1] + X^2*xoeffs[2] + ... + X^{n-1}*coeffs[n-1]
    */
    typedef struct{
        int16_t coeffs[256];
    } poly;

    typedef struct {
        poly vec[KYBER_K];
    } polyvec;
    



    int CheckKeySize(unsigned int key[], unsigned int keySize);
    void IndCpaKeypair(uint8_t pk[], uint8_t sk[]);
    
    
    unsigned int RejUniform(int16_t *r, unsigned int len, const uint8_t *buf, unsigned int buflen);

    void GenMatrix(polyvec *a, const uint8_t seed[], int transposed);
    typedef struct {
        uint64_t s[25];
    } keccak_state;

    void Shake128Absorb(keccak_state *state, const uint8_t *in, int inlen);
    void KeccakAbsorb(uint64_t s[25], unsigned int r, const uint8_t *m, int mlen, uint8_t p);
    void Shake128SqueezeBlocks(uint8_t *out, size_t nblocks, keccak_state *state);



//Handles the polynomials, functions can be found in kyber_poly.cpp
    void PolyCompress(uint8_t r[KYBER_POLYCOMPRESSEDBYTES], poly *a);
    void PolyDecompress(poly *r, const uint8_t a[KYBER_POLYCOMPRESSEDBYTES]);
    void PolyToBytes(uint8_t r[KYBER_POLYBYTES], poly *a);
    void PolyFromBytes(poly *r, const uint8_t a[KYBER_POLYBYTES]);
    void PolyFromMsg(poly *r, const uint8_t msg[KYBER_INDCPA_MSGBYTES]);
    void PolyToMsg(uint8_t msg[KYBER_INDCPA_MSGBYTES], poly *a);
    void PolyGetNoise(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce);
    void PolyNtt(poly *r);
    void PolyInvNttToMont(poly *r);
    void PolyBasemulMontgomery(poly *r, const poly *a, const poly *b);
    void PolyToMont(poly *r);
    void PolyReduce(poly *r);
    void PolyCsubq(poly *r);
    void PolyAdd(poly *r, const poly *a, const poly *b);
    void PolySub(poly *r, const poly *a, const poly *b);

    void PolyvecCompress(uint8_t r[KYBER_POLYVECCOMPRESSEDBYTES], polyvec *a);
    void PolyvecDecompress(polyvec *r,
                        const uint8_t a[KYBER_POLYVECCOMPRESSEDBYTES]);
    void PolyvecToBytes(uint8_t r[KYBER_POLYVECBYTES], polyvec *a);
    void PolyvecFromBytes(polyvec *r, const uint8_t a[KYBER_POLYVECBYTES]);
    void PolyvecNtt(polyvec *r);
    void PolyvecInvNttToMont(polyvec *r);
    void PolyvecPointwiseAccMontgomery(poly *r,
                                      const polyvec *a,
                                      const polyvec *b);
    void PolyvecReduce(polyvec *r);
    void PolyvecCsubq(polyvec *r);
    void PolyvecAdd(polyvec *r, const polyvec *a, const polyvec *b);

    //Ntt
    static const int16_t zetas[128];
    static const int16_t zetas_inv[128];

    int16_t Fqmul(int16_t a, int16_t b);
    void Ntt(int16_t r[256]);
    void InvNtt(int16_t r[256]);
    void Basemul(int16_t r[2], const int16_t a[2], const int16_t b[2], int16_t zeta);


    //Reduce
    int16_t MontgomeryReduce(int32_t a);
    int16_t BarrettReduce(int16_t a);
    int16_t Csubq(int16_t a);

    void Shake256Prf(uint8_t *out, size_t outlen, const uint8_t key[KYBER_SYMBYTES], uint8_t nonce);
    void Cbd(poly *r, const uint8_t buf[KYBER_ETA*KYBER_N/4]);

    void PackPk(uint8_t r[KYBER_INDCPA_PUBLICKEYBYTES],
                    polyvec *pk,
                    const uint8_t seed[KYBER_SYMBYTES]);
    void UnpackPk(polyvec *pk,
                      uint8_t seed[KYBER_SYMBYTES],
                      const uint8_t packedpk[KYBER_INDCPA_PUBLICKEYBYTES]);
    void PackSk(uint8_t r[KYBER_INDCPA_SECRETKEYBYTES], polyvec *sk);
    void UnpackSk(polyvec *sk, const uint8_t packedsk[KYBER_INDCPA_SECRETKEYBYTES]);
    void PackCiphertext(uint8_t r[KYBER_INDCPA_BYTES],
                            polyvec *b,
                            poly *v);
    void UnpackCiphertext(polyvec *b,
                              poly *v,
                              const uint8_t c[KYBER_INDCPA_BYTES]);

    



};


   
    
  



//Kyber-512
class Kyber512 : public Kyber
{
public:
    Kyber512() : Kyber(2) {};
};

//Kyber-768
class Kyber768 : public Kyber
{
public:
    Kyber768() : Kyber(3) {};

};

//Kyber-1024
class Kyber1024 : public Kyber
{
public:
    Kyber1024() : Kyber(4) {};

};


NAMESPACE_END

#endif
