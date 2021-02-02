/* The header file for CRYSTALS-Dilithium. Adapted by Julius Hekkala 
from the public-domain reference
implementation of Dilithium by the CRYSTALS team 
(https://github.com/pq-crystals/dilithium) 
 */
#ifndef CRYPTOPP_DILITHIUM_H
#define CRYPTOPP_DILITHIUM_H

#include "cryptlib.h"
#include <array>


NAMESPACE_BEGIN(CryptoPP)

#define SHAKE128_RATE (168)
#define SHAKE256_RATE (136)
#define SHA3_128_RATE (72)
#define SHA3_256_RATE (168)

class Dilithium
{
public:

    //Constructor. This is used by classes extending Dilithium (practically different modes of Dilithium)
    //Constructor is defined in dilithium.cpp
    Dilithium(byte mode, byte eta, word16 beta, word16 omega, word16 polyEtaPacked,
        word16 publicKeyBytes, word16 secretKeyBytes, word16 bytes);

    
    int Keypair(byte *pk, byte *sk);
    int Signature(byte *sig, size_t *sigLen, const byte *m, size_t mLen, const byte *sk);
    int Sign(byte *sm, size_t *smLen, const byte *m, size_t mLen, const byte *sk);
    int Verify(const byte *sig, size_t sigLen, const byte *m, size_t mlen, const byte *pk);
    int Open(byte *m, size_t *mLen, const byte *sm, size_t smLen, const byte *pk);

    void RandomBytes(byte *out, size_t outlen);


private:
    //Dilithium parameters
    byte mSeedBytes = 32;
    byte mCrhBytes = 48;
    word16 mN = 256;
    word32 mQ = 8380417;
    word16 mRootOfUnity = 1753;
    byte mD = 14;
    word32 mGamma1 = (mQ - 1) / 16;
    word32 mGamma2 = mGamma1 / 2;
    word32 mAlpha = 2 * mGamma2;
    byte mK;
    byte mL;
    byte mEta;
    word16 mBeta;
    byte mOmega;
    word16 mPolyt1PackedBytes = 288;
    word16 mPolyt0PackedBytes = 448;
    word16 mPolyzPackedBytes = 640;
    byte mPolyw1PackedBytes = 128;
    word16 mPolyEtaPackedBytes;
    word16 mPublicKeyBytes;
    word16 mSecretKeyBytes;
    word16 mBytes;
    word32 mMont = 4193792;
    word64 mQInv = 4236238847;


    //Polynomial stuff
    typedef std::array<word32, 256> poly;
    typedef std::vector<poly> polyvec;

    //Keccak state struct 
    typedef struct {
        word64 s[25];
        word32 pos;
    } keccakState;

    void Challenge(poly *c, const byte *mu, const polyvec *w1);


    //Utility
    word32 Power2Round(word32 a, word32 *a0);
    word32 Decompose(word32 a, word32 *a0);
    word32 MakeHint(word32 a0, word32 a1);
    word32 UseHint(word32 a, word32 hint);

    void PackPk(byte *pk, const byte *rho, const polyvec *t1);
    void UnpackPk(byte *rho, polyvec *t1, const byte *pk);
    void PackSk(byte *sk, const byte *rho, const byte *key, const byte *tr, const polyvec *s1,
             const polyvec *s2, const polyvec *t0);
    void UnpackSk(byte *rho, byte *key, byte *tr, polyvec *s1, polyvec *s2, polyvec *t0, const byte *sk);
    void PackSig(byte *sig, const polyvec *z, const polyvec *h, const poly *c);
    int UnpackSig(polyvec *z, polyvec *h, poly *c, const byte *sig);

    word32 MontgomeryReduce(word64 a);
    word32 Reduce32(word32 a);
    word32 CSubQ(word32 a);
    word32 Freeze(word32 a);

    //Custom SHAKE 
    void Shake128StreamInit(keccakState *state, const byte *seed, word16 nonce);
    void Shake256StreamInit(keccakState *state, const byte *seed, word16 nonce);
    void KeccakInit(keccakState *state);
    word32 KeccakAbsorb(word64 s[25], word32 r, word32 pos, const byte *m, size_t mlen);
    void KeccakFinalize(word64 s[25], word32 r, word32 pos, byte p);
    void KeccakSqueezeBlocks(byte *out, size_t nBlocks, word64 s[25], word32 r);
    word32 KeccakSqueeze(byte *out, size_t outlen, word64 s[25], word32 r, word32 pos);
    void Shake128Init(keccakState *state);
    void Shake128Absorb(keccakState *state, const byte *in, size_t inlen);
    void Shake128Finalize(keccakState *state);
    void Shake128SqueezeBlocks(byte *out, size_t nBlocks, keccakState *state);
    void Shake128Squeeze(byte *out, size_t outLen, keccakState *state);
    void Shake256Init(keccakState *state);
    void Shake256Absorb(keccakState *state, const byte *in, size_t inLen);
    void Shake256Finalize(keccakState *state);
    void Shake256SqueezeBlocks(byte *out, size_t nBlocks, keccakState *state);
    void Shake256Squeeze(byte *out, size_t outLen, keccakState *state);
    



    //Polynomial implementations can be found in dilithium_poly.cpp
    //Poly
    void PolyReduce(poly *a);
    void PolyCsubQ(poly *a);
    void PolyFreeze(poly *a);
    void PolyAdd(poly *c, const poly *a, const poly *b);
    void PolySub(poly *c, const poly *a, const poly *b);
    void PolyShiftL(poly *a);
    void PolyNtt(poly *a);
    void PolyInvNttToMont(poly *a);
    void PolyPointwiseMontgomery(poly *c, const poly *a, const poly *b);
    void PolyPower2Round(poly *a1, poly *a0, const poly *a);
    void PolyDecompose(poly *a1, poly *a0, const poly *a);
    word32 PolyMakeHint(poly *h, const poly *a0, const poly *a1);
    void PolyUseHint(poly *b, const poly *a, const poly *h);
    int PolyChkNorm(const poly *a, word32 B);
    word32 RejUniform(word32 *a, word32 len, const byte *buf, word32 buflen);
    void PolyUniform(poly *a, const byte *seed, word16 nonce);
    word32 RejEta(word32 *a, word32 len, const byte *buf, word32 buflen);
    void PolyUniformEta(poly *a, const byte *seed, word16 nonce);
    word32 RejGamma1m1(word32 *a, word32 len, const byte *buf, word32 buflen);
    void PolyUniformGamma1m1(poly *a, const byte *seed, word16 nonce);
    void PolyEtaPack(byte *r, const poly *a);
    void PolyEtaUnpack(poly *r, const byte *a);
    void Polyt1Pack(byte *r, const poly *a);
    void Polyt1Unpack(poly *r, const byte *a);
    void Polyt0Pack(byte *r, const poly *a);
    void Polyt0Unpack(poly *r, const byte *a);
    void PolyzPack(byte *r, const poly *a);
    void PolyzUnpack(poly *r, const byte *a);
    void Polyw1Pack(byte *r, const poly *a);
    //Polyvec
    void ExpandMat(std::vector<polyvec> *mat, const byte *rho);
    void PolyvecFreeze(polyvec *v, byte vectorLen);
    void PolyvecAdd(polyvec *w, const polyvec *u, const polyvec *v, byte vectorLen);
    void PolyvecNtt(polyvec *v, byte vectorLen);
    void PolyvecPointwiseAccMontgomery(poly *w, const polyvec *u, const polyvec *v, byte vectorLen);
    int PolyvecChkNorm(const polyvec *v, word32 bound, byte vectorLen);
    void PolyvecReduce(polyvec *v, byte vectorLen);
    void PolyvecCsubQ(polyvec *v, byte vectorLen);
    void PolyvecSub(polyvec *w, const polyvec *u, const polyvec *v, byte vectorLen);
    void PolyvecShiftL(polyvec *v, byte vectorLen);
    void PolyvecInvNttToMont(polyvec *v, byte vectorLen);
    void PolyvecPower2Round(polyvec *v1, polyvec *v0, const polyvec *v, byte vectorLen);
    void PolyvecDecompose(polyvec *v1, polyvec *v0, const polyvec *v, byte vectorLen);
    word32 PolyvecMakeHint(polyvec *h, const polyvec *v0, const polyvec *v1, byte vectorLen);
    void PolyvecUseHint(polyvec *w, const polyvec *u, const polyvec *h, byte vectorLen);

    //Ntt
    void Ntt(word32 *p);
    void InvNttToMont(word32 *p);



    

 


};


//Different modes of Dilithium

class Dilithium1 : public Dilithium {
public:
    //Parameters for mode 1 dilithium
    Dilithium1() : Dilithium(1, 7, 375, 64, 128, 896, 2096, 1387) {};

    //Parameters for the help of the user
    CRYPTOPP_CONSTANT(BYTES = 1387);
    CRYPTOPP_CONSTANT(PUBLICKEYBYTES = 896);
    CRYPTOPP_CONSTANT(SECRETKEYBYTES = 2096);

};

class Dilithium2: public Dilithium {
public:
    //Parameters for mode 2 dilithum
    Dilithium2() : Dilithium(2, 6, 325, 80, 128, 1184, 2800, 2044) {};

    //Parameters for the help of the user
    CRYPTOPP_CONSTANT(BYTES = 2044);
    CRYPTOPP_CONSTANT(PUBLICKEYBYTES = 1184);
    CRYPTOPP_CONSTANT(SECRETKEYBYTES = 2800);
};

class Dilithium3: public Dilithium {
public:
    //Paramters for mode 3 Dilithium
    Dilithium3() : Dilithium(3, 5, 275, 96, 128, 1472, 3500, 2701) {};

    //Parameters for the help of the user
    CRYPTOPP_CONSTANT(BYTES = 2701);
    CRYPTOPP_CONSTANT(PUBLICKEYBYTES = 1472);
    CRYPTOPP_CONSTANT(SECRETKEYBYTES = 3500);
};

class Dilithium4: public Dilithium {
public:
    //Parametrs for mode 4 Dilithium
    Dilithium4() : Dilithium(4, 3, 175, 120, 96, 1760, 3856, 3366) {};

    //Parameters for the help of the user
    CRYPTOPP_CONSTANT(BYTES = 3366);
    CRYPTOPP_CONSTANT(PUBLICKEYBYTES = 1760);
    CRYPTOPP_CONSTANT(SECRETKEYBYTES = 3856);
};



NAMESPACE_END

#endif
