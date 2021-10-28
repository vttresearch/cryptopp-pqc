/* The header file for CRYSTALS-Kyber. Adapted from the reference
implementation of kyber by the CRYSTALS team 
(https://github.com/pq-crystals/kyber) */

#ifndef CRYPTOPP_KYBER_H
#define CRYPTOPP_KYBER_H

#include "cryptlib.h"
#include "pqc_fips202.h"

NAMESPACE_BEGIN(CryptoPP)


template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
class Kyber : public Algorithm, public PQCFips202
{
public:
    

    //Kyber parameters that the user needs to be able to use 
    CRYPTOPP_CONSTANT(PUBLICKEYBYTES = T_K * 384 + 32);
    CRYPTOPP_CONSTANT(SECRETKEYBYTES = T_K * 384 + T_K * 384 + 32 + 2 * 32);
    CRYPTOPP_CONSTANT(CIPHERTEXTBYTES = T_C2 + T_K * T_C);
    CRYPTOPP_CONSTANT(SHAREDSECRETBYTES = 32);

    //The user should use these
    sword32 KemKeypair(byte *pk, byte *sk);
    sword32 KemEnc(byte *ct, byte *ss, const byte *pk);
    sword32 KemDec(byte *ss, const byte *ct, const byte *sk);


    //Public for testing purposes
    void RandomBytes(byte *out, size_t outLen);

    void IndCpaKeypair(byte *pk, byte *sk);
    void IndCpaEnc(byte *c, const byte *m, const byte *pk, const byte *coins);
    void IndCpaDec(byte *m, const byte *c, const byte *sk);

    
protected:
    /*
    * Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
    * coeffs[0] + X*coeffs[1] + X^2*coeffs[2] + ... + X^{n-1}*coeffs[n-1]
    */
    typedef struct{
        sword16 coeffs[256];
    } poly;

    typedef struct {
        poly vec[T_K];
    } polyvec;
    
    word32 RejUniform(sword16 *r, word32 len, const byte *buf, word32 buflen);

    void GenMatrix(polyvec *a, const byte *seed, sword32 transposed);

    //Handles the polynomials, functions can be found in kyber_poly.cpp
    void PolyCompress(byte *r, poly *a);
    void PolyDecompress(poly *r, const byte *a);
    void PolyToBytes(byte *r, const poly *a);
    void PolyFromBytes(poly *r, const byte *a);
    void PolyFromMsg(poly *r, const byte *msg);
    void PolyToMsg(byte *msg, poly *a);
    void PolyGetNoise(poly *r, const byte *seed, byte nonce, const sword32 eta);
    void PolyNtt(poly *r);
    void PolyInvNttToMont(poly *r);
    void PolyBasemulMontgomery(poly *r, const poly *a, const poly *b);
    void PolyToMont(poly *r);
    void PolyReduce(poly *r);
    void PolyAdd(poly *r, const poly *a, const poly *b);
    void PolySub(poly *r, const poly *a, const poly *b);

    void PolyvecCompress(byte *r, polyvec *a);
    void PolyvecDecompress(polyvec *r, const byte *a);
    void PolyvecToBytes(byte *r, const polyvec *a);
    void PolyvecFromBytes(polyvec *r, const byte *a);
    void PolyvecNtt(polyvec *r);
    void PolyvecInvNttToMont(polyvec *r);
    void PolyvecBasemulAccMontgomery(poly *r, const polyvec *a, const polyvec *b);
    void PolyvecReduce(polyvec *r);
    void PolyvecAdd(polyvec *r, const polyvec *a, const polyvec *b);

    //Ntt
    static const sword16 zetas[128];
    static const sword16 zetas_inv[128];

    sword16 Fqmul(sword16 a, sword16 b);
    void Ntt(sword16 r[256]);
    void InvNtt(sword16 r[256]);
    void Basemul(sword16 r[2], const sword16 a[2], const sword16 b[2], sword16 zeta);


    //Reduce
    sword16 MontgomeryReduce(sword32 a);
    sword16 BarrettReduce(sword16 a);

    void KyberShake256PRF(byte *output, size_t outputLen, const byte *key, byte nonce);
    void KyberShake128Absorb(const byte *seed, byte x, byte y);

    word32 Load32LittleEndian(const byte x[4]);
    word32 Load24LittleEndian(const byte x[3]);

    void Cbd(poly *r, const byte *buf, const byte eta);

    void PackPk(byte *r, polyvec *pk, const byte *seed);
    void UnpackPk(polyvec *pk, byte *seed, const byte *packedpk);
    void PackSk(byte *r, polyvec *sk);
    void UnpackSk(polyvec *sk, const byte *packedsk);
    void PackCiphertext(byte *r, polyvec *b, poly *v);
    void UnpackCiphertext(polyvec *b, poly *v, const byte *c);


    void Cmov(byte *r, const byte *x, size_t len, byte b);

    //Kyber parameters that do not need to be public
    CRYPTOPP_CONSTANT(K = T_K);
    CRYPTOPP_CONSTANT(N = 256);
    CRYPTOPP_CONSTANT(Q = 3329);
    CRYPTOPP_CONSTANT(ETA1 = T_ETA1);
    CRYPTOPP_CONSTANT(ETA2 = 2);
    CRYPTOPP_CONSTANT(SYMBYTES = 32);
    CRYPTOPP_CONSTANT(POLYBYTES = 384);
    CRYPTOPP_CONSTANT(POLYVECBYTES = T_K * 384);
    CRYPTOPP_CONSTANT(POLYCOMPRESSEDBYTES = T_C2);
    CRYPTOPP_CONSTANT(POLYVECCOMPRESSEDBYTES = T_K * T_C);
    CRYPTOPP_CONSTANT(INDCPA_MSGBYTES = 32);
    CRYPTOPP_CONSTANT(INDCPA_PUBLICKEYBYTES = T_K * 384 + 32);
    CRYPTOPP_CONSTANT(INDCPA_SECRETKEYBYTES = T_K * 384);
    CRYPTOPP_CONSTANT(INDCPA_BYTES = T_C2 + T_K * T_C);
    CRYPTOPP_CONSTANT(MONT = -1044);
    CRYPTOPP_CONSTANT(QINV = -3327);
};

  



//Kyber-512
class Kyber512 : public Kyber<2,320,128,3>
{
public:
    //Return the name of the algorithm
    std::string AlgorithmName() const {return "CRYSTALS-Kyber-512";}
};

//Kyber-768
class Kyber768 : public Kyber<3,320,128,2>
{
public:
    //Return the name of the algorithm
    std::string AlgorithmName() const {return "CRYSTALS-Kyber-768";}    
};

//Kyber-1024
class Kyber1024 : public Kyber<4,352,160,2>
{
public:
    //Return the name of the algorithm
    std::string AlgorithmName() const {return "CRYSTALS-Kyber-1024";}
};


NAMESPACE_END

#endif
