/* The header file for FrodoKEM. Adapted
from the implementation of FrodoKEM by the FrodoKEM team and Microsoft  
(https://github.com/Microsoft/PQCrypto-LWEKE)

Original license 
MIT License

    Copyright (c) Microsoft Corporation. All rights reserved.

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE
 */


#ifndef CRYPTOPP_FRODO_KEM_H
#define CRYPTOPP_FRODO_KEM_H

#include "cryptlib.h"
#include "pqc_fips202.h"

// Definition of operating system

#define OS_WIN       1
#define OS_NIX       2

#if defined(__WIN64)            // Microsoft Windows
    #define OS_TARGET OS_WIN
#else              // Unix-like operative systems
    #define OS_TARGET OS_NIX 
#endif


// Definition of compiler

#define COMPILER_VC      1
#define COMPILER_GCC     2
#define COMPILER_CLANG   3

#if defined(_MSC_VER)           // Microsoft Visual C compiler
    #define COMPILER COMPILER_VC
#elif defined(__GNUC__)         // GNU GCC compiler
    #define COMPILER COMPILER_GCC   
#elif defined(__clang__)        // Clang compiler
    #define COMPILER COMPILER_CLANG
#else
    #error -- "Unsupported COMPILER"
#endif


// Definition of the targeted architecture and basic data types
    
#define TARGET_AMD64        1
#define TARGET_x86          2
#define TARGET_ARM          3
#define TARGET_PPC          4
#define TARGET_S390X        5


#if (defined(__arm64__) || defined(__x64__) || defined(__x86_64__)) 
    #define TARGET TARGET_AMD64
#elif defined(__x86__)
    #define TARGET TARGET_x86
#endif


#if defined(__WIN64)
    #define ALIGN_HEADER(N) __declspec(align(N))
    #define ALIGN_FOOTER(N) 
#else
    #define ALIGN_HEADER(N)
    #define ALIGN_FOOTER(N) __attribute__((aligned(N)))
#endif



// Configuration for endianness
#if (TARGET == TARGET_PPC || TARGET == TARGET_S390X) || (defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__))
    #define LE_TO_UINT16(n) (((((unsigned short)(n) & 0xFF)) << 8) | (((unsigned short)(n) & 0xFF00) >> 8))
    #define UINT16_TO_LE(n) (((((unsigned short)(n) & 0xFF)) << 8) | (((unsigned short)(n) & 0xFF00) >> 8))
#elif (TARGET == TARGET_x86 || TARGET == TARGET_AMD64) || (defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__))
    #define LE_TO_UINT16(n) (n)
    #define UINT16_TO_LE(n) (n)
#else
    #define LE_TO_UINT16(n) (((uint8_t *) &(n))[0] | (((uint8_t *) &(n))[1] << 8))
    static inline uint16_t UINT16_TO_LE(const uint16_t x) {
        uint16_t y;
        uint8_t *z = (uint8_t *) &y;
        z[0] = x & 0xFF;
        z[1] = x >> 8;
        return y;
    }
#endif

NAMESPACE_BEGIN(CryptoPP)

template<word16 T_SK_BYTES, word16 T_PK_BYTES, byte T_BYTES, word16 T_CT_BYTES, word16 T_N, byte T_LOGQ, byte T_EB>
class FrodoKEM : public Algorithm, public PQCFips202
{
    public: 
        //User needs to be able to use these
        CRYPTOPP_CONSTANT(SECRETKEYBYTES=T_SK_BYTES);
        CRYPTOPP_CONSTANT(PUBLICKEYBYTES=T_PK_BYTES);
        CRYPTOPP_CONSTANT(BYTES=T_BYTES);
        CRYPTOPP_CONSTANT(CIPHERTEXTBYTES=T_CT_BYTES);


         //The user should use these
        int KemKeypair(byte *pk, byte *sk);
        int KemEnc(byte *ct, byte *ss, const byte *pk);
        int KemDec(byte *ss, const byte *ct, const byte *sk);

        //Public for testing purposes
        void RandomBytes(byte *out, size_t outLen);
    
    protected:
        CRYPTOPP_CONSTANT(N = T_N);
        CRYPTOPP_CONSTANT(NBAR = 8);
        CRYPTOPP_CONSTANT(LOGQ = T_LOGQ);
        CRYPTOPP_CONSTANT(Q =  (1 << T_LOGQ));
        CRYPTOPP_CONSTANT(EXTRACTED_BITS = T_EB);
        CRYPTOPP_CONSTANT(STRIPE_STEP = 8);
        CRYPTOPP_CONSTANT(PARALLEL = 4);
        CRYPTOPP_CONSTANT(BYTES_SEED_A = 16);
        CRYPTOPP_CONSTANT(BYTES_MU = (T_EB*NBAR*NBAR)/8);
        CRYPTOPP_CONSTANT(BYTES_PKHASH = T_BYTES);


        int FrodoMulAddAsPlusE(word16 *out, const word16 *s, const word16 *e, const byte *seed_A);
        int FrodoMulAddsAPlusE(word16 *out, const word16 *s, word16 *e, const byte *seedA); 
        void FrodoMulBs(word16 *out, const word16 *b, const word16 *s);
        void FrodoMulAddSbPlusE(word16 *out, const word16 *b, const word16 *s, const word16 *e);
        void FrodoAdd(word16 *out, const word16 *a, const word16 *b);
        void FrodoSub(word16 *out, const word16 *a, const word16 *b);
        void FrodoKeyEncode(word16 *out, const word16 *in);
        void FrodoKeyDecode(word16 *out, const word16 *in);
        void SampleN(word16 *s, const size_t n, word16 version);
        void Pack(byte *out, const size_t outLen, const word16 *in, const size_t inLen, const byte lsb);
        void Unpack(word16 *out, const size_t outLen, const byte *in, const size_t inLen, const byte lsb);
        sbyte CtVerify(const word16 *a, const word16 *b, size_t len);
        void CtSelect(byte *r, const byte *a, const byte *b, size_t len, sbyte selector);
        void ClearBytes(byte *mem, size_t n);


};


class Frodo640 : public FrodoKEM<19888, 9616, 16, 9720, 640, 15, 2> 
{
    public:
        //Return the name of the algorithm
        std::string AlgorithmName() const {return "FrodoKEM-640";}
};

class Frodo976 : public FrodoKEM<31296, 15632, 24, 15744, 976, 16, 3> 
{
    public:
        //Return the name of the algorithm
        std::string AlgorithmName() const {return "FrodoKEM-976";}
};

class Frodo1344 : public FrodoKEM<43088, 21520, 32, 21632, 1344, 16, 4> 
{
    public:
        //Return the name of the algorithm
        std::string AlgorithmName() const {return "FrodoKEM-1344";}
};

NAMESPACE_END

#endif