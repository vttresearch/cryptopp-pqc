/* The header file for the post quantum KEM algorithm SABER. Adapted by Julius Hekkala 
from the public-domain reference
implementation of SABER by the SABER team 
(https://github.com/KULeuven-COSIC/SABER) 
 */
#ifndef CRYPTOPP_SABER_H
#define CRYPTOPP_SABER_H

#include "cryptlib.h"

NAMESPACE_BEGIN(CryptoPP)


/**
 * Template class to support different modes of SABER
 */
template<byte T_L, byte T_Et>
class SABER {
public:
    //Generate a keypair consisting of public and secret key
    int KemKeypair(byte *pk, byte *sk);

    //Key encapsulation
    int KemEnc(byte *ct, byte *ss, const byte *pk);

    //Key decapsulation
    int KemDec(byte *ss, const byte *ct, const byte *sk);


private:
    //SABER parameters
    CRYPTOPP_CONSTANT(L=T_L);
    CRYPTOPP_CONSTANT(ET=T_Et);
    CRYPTOPP_CONSTANT(MU=10 - (T_L - 2) * 2);
    CRYPTOPP_CONSTANT(EQ=13);
    CRYPTOPP_CONSTANT(EP=10);
    CRYPTOPP_CONSTANT(N=256);
    CRYPTOPP_CONSTANT(SEEDBYTES=32);
    CRYPTOPP_CONSTANT(NOISE_SEEDBYTES=32);
    CRYPTOPP_CONSTANT(KEYBYTES=32);
    CRYPTOPP_CONSTANT(SHAREDSECRETBYTES=KEYBYTES);
    CRYPTOPP_CONSTANT(HASHBYTES=32);
    CRYPTOPP_CONSTANT(POLYCOINBYTES=MU * N / 8);
    CRYPTOPP_CONSTANT(POLYBYTES= EQ * N / 8);
    CRYPTOPP_CONSTANT(POLYVECBYTES=L*POLYBYTES);
    CRYPTOPP_CONSTANT(POLYCOMPRESSEDBYTES=EP * N /8);
    CRYPTOPP_CONSTANT(POLYVECCOMPRESSEDBYTES=L * POLYCOMPRESSEDBYTES);
    CRYPTOPP_CONSTANT(SCALEBYTES_KEM = ET * N / 8);
    CRYPTOPP_CONSTANT(INDCPA_PUBLICKEYBYTES=POLYVECCOMPRESSEDBYTES + SEEDBYTES);
    CRYPTOPP_CONSTANT(INDCPA_SECRETKEYBYTES= POLYVECBYTES);
    CRYPTOPP_CONSTANT(PUBLICKEYBYTES=INDCPA_PUBLICKEYBYTES);
    CRYPTOPP_CONSTANT(SECRETKEYBYTES=INDCPA_SECRETKEYBYTES + INDCPA_PUBLICKEYBYTES + HASHBYTES + KEYBYTES);
    CRYPTOPP_CONSTANT(CCA_DEC = POLYVECCOMPRESSEDBYTES + SCALEBYTES_KEM);
    CRYPTOPP_CONSTANT(CIPHERTEXTBYTES = CCA_DEC);
    CRYPTOPP_CONSTANT(h1 = 1 << (EQ - EP - 1));
    CRYPTOPP_CONSTANT(h2 = (1 << (EP - 2)) - (1 << (EP - ET - 1)) + (1 << (EQ - EP - 1)));

    CRYPTOPP_CONSTANT(N_RES = N << 1);
    CRYPTOPP_CONSTANT(N_SB = N >> 2);
    CRYPTOPP_CONSTANT(N_SB_RES = 2 * N_SB - 1);

    
    void IndCpaKemKeypair(byte *pk, byte *sk);
    void IndCpaKemEnc(const byte *m, const byte *seedSp, const byte *pk, byte *cipherText);
    void IndCpaKemDec(const byte *sk, const byte *cipherText, byte *m);



    void RandomBytes(byte *out, word64 outLen);
    void Cmov(byte *r, const byte *x, size_t len, byte b);

    void MatrixVectorMul(const word16 A[L][L][N], const word16 s[L][N], word16 res[L][N], word16 transpose);
    void InnerProd(const word16 b[L][N], const word16 s[L][N], word16 *res);
    void GenMatrix(word16 a[L][L][N], const byte seed[SEEDBYTES]);
    void GenSecret(word16 s[L][N], const byte seed[NOISE_SEEDBYTES]);
    void PolyMulAcc(const word16 *a, const word16 *b, word16 *res);
    void KaratsubaSimple(const word16 *a1, const word16 *b1, word16 *resultFinal);
    void ToomCook4Way (const word16 *a1, const word16 *b1, word16 *result);

    void Cbd(word16 *s, const byte *buf);


    //Packing utils
    void POLT2BS(byte *bytes, const word16 *data);
    void BS2POLT(const byte *bytes, word16 *data);
    void POLq2BS(byte *bytes, const word16 *data);
    void BS2POLq(const byte *bytes, word16 *data);
    void POLp2BS(byte *bytes, const word16 *data);
    void BS2POLp(const byte *bytes, word16 *data);
    void POLVECq2BS(byte *bytes, const word16 *data);
    void BS2POLVECq(const byte bytes[POLYVECBYTES], word16 data[L][N]);
    void POLVECq2BS(byte bytes[POLYVECBYTES], const word16 data[L][N]);
    void POLVECp2BS(byte bytes[POLYVECCOMPRESSEDBYTES], const word16 data[L][N]);
    void BS2POLVECp(const byte bytes[POLYVECCOMPRESSEDBYTES], word16 data[L][N]);
    void BS2POLmsg(const byte bytes[KEYBYTES], word16 data[N]);
    void POLmsg2BS(byte bytes[KEYBYTES], const word16 data[N]);





};


//Different modes of SABER

class LightSaber : public SABER<2,3> {
public:
    //Parameters for the user
    CRYPTOPP_CONSTANT(PUBLICKEYBYTES = 672);
    CRYPTOPP_CONSTANT(SECRETKEYBYTES = 1568);
    CRYPTOPP_CONSTANT(CIPHERTEXTBYTES = 736);
    CRYPTOPP_CONSTANT(SHAREDSECRETBYTES = 32);
};

class Saber : public SABER<3,4> {
public:
    //Parameters for the user
    CRYPTOPP_CONSTANT(PUBLICKEYBYTES = 992);
    CRYPTOPP_CONSTANT(SECRETKEYBYTES = 2304);
    CRYPTOPP_CONSTANT(CIPHERTEXTBYTES = 1088);
    CRYPTOPP_CONSTANT(SHAREDSECRETBYTES = 32);
};

class FireSaber : public SABER<4,6> {
public:
    //Parameters for the user
    CRYPTOPP_CONSTANT(PUBLICKEYBYTES = 1312);
    CRYPTOPP_CONSTANT(SECRETKEYBYTES = 3040);
    CRYPTOPP_CONSTANT(CIPHERTEXTBYTES = 1472);
    CRYPTOPP_CONSTANT(SHAREDSECRETBYTES = 32);
};



























NAMESPACE_END

#endif