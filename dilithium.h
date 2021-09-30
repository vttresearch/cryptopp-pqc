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
    Dilithium(byte k, byte l, byte eta, byte tau, word16 beta, word32 gamma1, word32 gamma2, word16 omega, word16 polyEtaPacked, 
        word16 polyzPackedBytes, byte polyw1PackedBytes);

    
    int Keypair(byte *pk, byte *sk);
    int Signature(byte *sig, size_t *sigLen, const byte *m, size_t mLen, const byte *sk);
    int Sign(byte *sm, size_t *smLen, const byte *m, size_t mLen, const byte *sk);
    int Verify(const byte *sig, size_t sigLen, const byte *m, size_t mlen, const byte *pk);
    int Open(byte *m, size_t *mLen, const byte *sm, size_t smLen, const byte *pk);

    void RandomBytes(byte *out, size_t outlen);
    //Set to randomized mode by calling this
    int setRandomizedSignature() {
        mRandomized = true;
        return 0;
    }

private:
    //Dilithium parameters
    byte mSeedBytes = 32;
    byte mCrhBytes = 64;
    word16 mN = 256;
    sword32 mQ = 8380417;
    word16 mRootOfUnity = 1753;
    byte mD = 13;
    sword32 mGamma1;
    sword32 mGamma2;
    byte mK;
    byte mL;
    byte mTau;
    byte mEta;
    word16 mBeta;
    byte mOmega;
    word16 mPolyt1PackedBytes = 320;
    word16 mPolyt0PackedBytes = 416;
    word16 mPolyzPackedBytes;
    byte mPolyw1PackedBytes;
    word16 mPolyEtaPackedBytes;
    word16 mPublicKeyBytes;
    word16 mSecretKeyBytes;
    word16 mBytes;
    sword32 mMont = -4186625;
    sword32 mQInv = 58728449;

    //Randomized signing. Deterministic by default
    bool mRandomized = false;

    

    //Polynomial stuff
    typedef std::array<sword32, 256> poly;
    typedef std::vector<poly> polyvec;

    //Keccak state struct 
    typedef struct {
        word64 s[25];
        word32 pos;
    } keccakState;

    void Challenge(poly *c, const byte *mu);


    //Utility
    sword32 Power2Round(sword32 *a0, sword32 a);
    sword32 Decompose(sword32 *a0, sword32 a);
    word32 MakeHint(sword32 a0, sword32 a1);
    sword32 UseHint(sword32 a, word32 hint);

    void PackPk(byte *pk, const byte *rho, const polyvec *t1);
    void UnpackPk(byte *rho, polyvec *t1, const byte *pk);
    void PackSk(byte *sk, const byte *rho, const byte *tr, const byte *key, const polyvec *t0,
             const polyvec *s1, const polyvec *s2);
    void UnpackSk(byte *rho, byte *tr, byte *key, polyvec *t0,  
                          polyvec *s1, polyvec *s2, const byte *sk);
    void PackSig(byte *sig, const byte *c, const polyvec *z, const polyvec *h);
    int UnpackSig(byte *c, polyvec *z, polyvec *h, const byte *sig);

    sword32 MontgomeryReduce(sword64 a);
    sword32 Reduce32(sword32 a);
    //word32 CSubQ(word32 a);
    sword32 Freeze(sword32 a);
    sword32 CAddQ(sword32 a);

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
    void PolyCAddQ(poly *a);
    void PolyShiftL(poly *a);
    void PolyNtt(poly *a);
    void PolyInvNttToMont(poly *a);
    void PolyPointwiseMontgomery(poly *c, const poly *a, const poly *b);
    void PolyPower2Round(poly *a1, poly *a0, const poly *a);
    void PolyDecompose(poly *a1, poly *a0, const poly *a);
    word32 PolyMakeHint(poly *h, const poly *a0, const poly *a1);
    void PolyUseHint(poly *b, const poly *a, const poly *h);
    int PolyChkNorm(const poly *a, sword32 B);
    word32 RejUniform(sword32 *a, word32 len, const byte *buf, word32 buflen);
    void PolyUniform(poly *a, const byte *seed, word16 nonce);
    word32 RejEta(sword32 *a, word32 len, const byte *buf, word32 buflen);
    void PolyUniformEta(poly *a, const byte *seed, word16 nonce);
    word32 RejGamma1m1(word32 *a, word32 len, const byte *buf, word32 buflen);
    void PolyUniformGamma1(poly *a, const byte *seed, word16 nonce);
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
    void Polyvecw1Pack(byte *r, const polyvec *h, const byte vectorLen);
    void ExpandMat(std::vector<polyvec> *mat, const byte *rho);
    void PolyvecFreeze(polyvec *v, byte vectorLen);
    void PolyvecUniformEta(polyvec *v, const byte *seed, word16 nonce, const byte vectorLen);
    void PolyvecUniformGamma1(polyvec *v, const byte *seed, word16 nonce);
    void PolyvecMatrixPointwiseMontgomery(polyvec *t, const polyvec *mat, const polyvec *v);
    void PolyvecAdd(polyvec *w, const polyvec *u, const polyvec *v, byte vectorLen);
    void PolyvecCAddQ(polyvec *v, byte vectorLen);
    void PolyvecNtt(polyvec *v, byte vectorLen);
    void PolyvecInvNttToMont(polyvec *v, const byte vectorLen);
    void PolyvecMatrixPointwiseMontgomery(polyvec *t, const std::vector<polyvec> *mat, const polyvec *v);
    void PolyvecPointwisePolyMontgomery(polyvec *r, const poly *a, const polyvec *v, const byte vectorLen); 
    void PolyvecPointwiseAccMontgomery(poly *w, const polyvec *u, const polyvec *v, byte vectorLen);
    int PolyvecChkNorm(const polyvec *v, sword32 bound, byte vectorLen);
    void PolyvecReduce(polyvec *v, byte vectorLen);
    void PolyvecCsubQ(polyvec *v, byte vectorLen);
    void PolyvecSub(polyvec *w, const polyvec *u, const polyvec *v, byte vectorLen);
    void PolyvecShiftL(polyvec *v, byte vectorLen);
    //void PolyvecInvNttToMont(polyvec *v, byte vectorLen);
    void PolyvecPower2Round(polyvec *v1, polyvec *v0, const polyvec *v, byte vectorLen);
    void PolyvecDecompose(polyvec *v1, polyvec *v0, const polyvec *v, byte vectorLen);
    word32 PolyvecMakeHint(polyvec *h, const polyvec *v0, const polyvec *v1, byte vectorLen);
    void PolyvecUseHint(polyvec *w, const polyvec *u, const polyvec *h, byte vectorLen);

    //Ntt
    void Ntt(sword32 *a);
    void InvNttToMont(sword32 *p);



    

 


};


//Different modes of Dilithium

class Dilithium2: public Dilithium {
public:
    //Parameters for mode 2 dilithum
    Dilithium2() : Dilithium(4, 4, 2, 39, 78, 1 << 17, 8380416 / 88, 80, 96, 576, 192) {};

    //Parameters for the help of the user
    CRYPTOPP_CONSTANT(BYTES = 2420);
    CRYPTOPP_CONSTANT(PUBLICKEYBYTES = 1312);
    CRYPTOPP_CONSTANT(SECRETKEYBYTES = 2528);
};

class Dilithium3: public Dilithium {
public:
    //Parameters for mode 3 Dilithium
    Dilithium3() : Dilithium(6, 5, 4, 49, 196, 1 << 19, 8380416/32, 55, 128, 640, 128) {};

    //Parameters for the help of the user
    CRYPTOPP_CONSTANT(BYTES = 3293);
    CRYPTOPP_CONSTANT(PUBLICKEYBYTES = 1952);
    CRYPTOPP_CONSTANT(SECRETKEYBYTES = 4000);
};

class Dilithium5: public Dilithium {
public:
    //Parameters for mode 5 Dilithium
    Dilithium5() : Dilithium(8, 7, 2, 60, 120, 1 << 19, 8380416/32, 75, 96, 640, 128) {};

    //Parameters for the help of the user
    CRYPTOPP_CONSTANT(BYTES = 4595);
    CRYPTOPP_CONSTANT(PUBLICKEYBYTES = 2592);
    CRYPTOPP_CONSTANT(SECRETKEYBYTES = 4864);
};



NAMESPACE_END

#endif
