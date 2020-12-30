/* The header file for CRYSTALS-Kyber. Adapted from the reference
implementation of kyber by the CRYSTALS team 
(https://github.com/pq-crystals/kyber) */

#ifndef CRYPTOPP_KYBER_H
#define CRYPTOPP_KYBER_H

#include "cryptlib.h"

NAMESPACE_BEGIN(CryptoPP)


template<int T_K, unsigned int T_Compr>
class Kyber : public Algorithm
{
public:
    //Kyber parameters that the user needs to be able to use 
    CRYPTOPP_CONSTANT(PUBLICKEYBYTES = T_K * 384 + 32);
    CRYPTOPP_CONSTANT(SECRETKEYBYTES = T_K * 384 + T_K * 384 + 32 + 2 * 32);
    CRYPTOPP_CONSTANT(CIPHERTEXTBYTES = 96 + (T_K - 2) * 32 + T_K * T_Compr);
    CRYPTOPP_CONSTANT(SHAREDSECRETBYTES = 32);

    //The user should use these
    int KemKeypair(unsigned char *pk, unsigned char *sk);
    int KemEnc(unsigned char *ct, unsigned char *ss, const unsigned char *pk);
    int KemDec(unsigned char *ss, const unsigned char *ct, const unsigned char *sk);


    //Public for testing purposes
    void randombytes(uint8_t *out, size_t outlen);

    void IndCpaKeypair(uint8_t *pk, uint8_t *sk);
    void IndCpaEnc(uint8_t *c, const uint8_t *m, const uint8_t *pk, const uint8_t *coins);
    void IndCpaDec(uint8_t *m, const uint8_t *c, const uint8_t *sk);

    
protected:
    /*
    * Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
    * coeffs[0] + X*coeffs[1] + X^2*xoeffs[2] + ... + X^{n-1}*coeffs[n-1]
    */
    typedef struct{
        int16_t coeffs[256];
    } poly;

    typedef struct {
        poly vec[T_K];
    } polyvec;

    typedef struct {
        uint64_t s[25];
    } keccak_state;
    
    unsigned int RejUniform(int16_t *r, unsigned int len, const uint8_t *buf, unsigned int buflen);

    void GenMatrix(polyvec *a, const uint8_t *seed, int transposed);

    //Implementations in kyber_utils.cpp
    void XOFAbsorb(keccak_state *state, const uint8_t *seed, uint8_t x, uint8_t y);
    void Shake128Absorb(keccak_state *state, const uint8_t *in, size_t inlen);
    void KeccakAbsorb(uint64_t s[25], unsigned int r, const uint8_t *m, int mlen, uint8_t p);
    void XOFSqueezeBlocks(uint8_t *out, size_t nblocks, keccak_state *state);
    void Shake128SqueezeBlocks(uint8_t *out, size_t nblocks, keccak_state *state);
    void KeccakSqueezeBlocks(uint8_t *out, size_t nblocks, uint64_t s[25], unsigned int r);

    //Handles the polynomials, functions can be found in kyber_poly.cpp
    void PolyCompress(uint8_t *r, poly *a);
    void PolyDecompress(poly *r, const uint8_t *a);
    void PolyToBytes(uint8_t *r, poly *a);
    void PolyFromBytes(poly *r, const uint8_t *a);
    void PolyFromMsg(poly *r, const uint8_t *msg);
    void PolyToMsg(uint8_t *msg, poly *a);
    void PolyGetNoise(poly *r, const uint8_t *seed, uint8_t nonce);
    void PolyNtt(poly *r);
    void PolyInvNttToMont(poly *r);
    void PolyBasemulMontgomery(poly *r, const poly *a, const poly *b);
    void PolyToMont(poly *r);
    void PolyReduce(poly *r);
    void PolyCsubq(poly *r);
    void PolyAdd(poly *r, const poly *a, const poly *b);
    void PolySub(poly *r, const poly *a, const poly *b);

    void PolyvecCompress(uint8_t *r, polyvec *a);
    void PolyvecDecompress(polyvec *r, const uint8_t *a);
    void PolyvecToBytes(uint8_t *r, polyvec *a);
    void PolyvecFromBytes(polyvec *r, const uint8_t *a);
    void PolyvecNtt(polyvec *r);
    void PolyvecInvNttToMont(polyvec *r);
    void PolyvecPointwiseAccMontgomery(poly *r, const polyvec *a, const polyvec *b);
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

    void Shake256Prf(uint8_t *out, size_t outlen, const uint8_t *key, uint8_t nonce);
    void Cbd(poly *r, const uint8_t *buf);

    void PackPk(uint8_t *r, polyvec *pk, const uint8_t *seed);
    void UnpackPk(polyvec *pk, uint8_t *seed, const uint8_t *packedpk);
    void PackSk(uint8_t *r, polyvec *sk);
    void UnpackSk(polyvec *sk, const uint8_t *packedsk);
    void PackCiphertext(uint8_t *r, polyvec *b, poly *v);
    void UnpackCiphertext(polyvec *b, poly *v, const uint8_t *c);


    void Cmov(uint8_t *r, const uint8_t *x, size_t len, uint8_t b);

    //Kyber parameters that do not need to be public
    CRYPTOPP_CONSTANT(K = T_K);
    CRYPTOPP_CONSTANT(N = 256);
    CRYPTOPP_CONSTANT(Q= 3329);
    CRYPTOPP_CONSTANT(ETA = 2);
    CRYPTOPP_CONSTANT(SYMBYTES = 32);
    CRYPTOPP_CONSTANT(POLYBYTES = 384);
    CRYPTOPP_CONSTANT(POLYVECBYTES = T_K * 384);
    CRYPTOPP_CONSTANT(POLYCOMPRESSEDBYTES = 96 + (T_K - 2) * 32);
    CRYPTOPP_CONSTANT(POLYVECCOMPRESSEDBYTES = T_K * T_Compr);
    CRYPTOPP_CONSTANT(INDCPA_MSGBYTES = 32);
    CRYPTOPP_CONSTANT(INDCPA_PUBLICKEYBYTES = T_K * 384 + 32);
    CRYPTOPP_CONSTANT(INDCPA_SECRETKEYBYTES = T_K * 384);
    CRYPTOPP_CONSTANT(INDCPA_BYTES = 96 + (T_K - 2) * 32 + T_K * T_Compr);
    CRYPTOPP_CONSTANT(MONT = 2285);
    CRYPTOPP_CONSTANT(QINV = 62209);

};

  



//Kyber-512
class Kyber512 : public Kyber<2,320>
{
public:
    //Return the name of the algorithm
    std::string AlgorithmName() const {return "CRYSTALS-Kyber-512";}
};

//Kyber-768
class Kyber768 : public Kyber<3,320>
{
public:
    //Return the name of the algorithm
    std::string AlgorithmName() const {return "CRYSTALS-Kyber-768";}    
};

//Kyber-1024
class Kyber1024 : public Kyber<4,352>
{
public:
    //Return the name of the algorithm
    std::string AlgorithmName() const {return "CRYSTALS-Kyber-1024";}
};


NAMESPACE_END

#endif
