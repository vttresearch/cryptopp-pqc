/* The C++ file for FrodoKEM. Adapted
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

#include "frodo_kem.h"
#include "osrng.h"

NAMESPACE_BEGIN(CryptoPP)

template<word16 T_SK_BYTES, word16 T_PK_BYTES, byte T_BYTES, word16 T_CT_BYTES, word16 T_N, byte T_LOGQ, byte T_EB>
int FrodoKEM<T_SK_BYTES, T_PK_BYTES, T_BYTES, T_CT_BYTES, T_N, T_LOGQ, T_EB>::KemKeypair(byte* pk, byte* sk)
{ // FrodoKEM's key generation
  // Outputs: public key pk (               BYTES_SEED_A + (PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8 bytes)
  //          secret key sk (CRYPTO_BYTES + BYTES_SEED_A + (PARAMS_LOGQ*PARAMS_N*PARAMS_NBAR)/8 + 2*PARAMS_N*PARAMS_NBAR + BYTES_PKHASH bytes)
    byte *pkSeedA = &pk[0];
    byte *pk_b = &pk[BYTES_SEED_A];
    byte *sk_s = &sk[0];
    byte *sk_pk = &sk[BYTES];
    byte *sk_S = &sk[BYTES + PUBLICKEYBYTES];
    byte *sk_pkh = &sk[BYTES + PUBLICKEYBYTES + 2*N*NBAR];
    word16 B[N*NBAR] = {0};
    word16 S[2*N*NBAR] = {0};               // contains secret data
    word16 *E = (word16 *)&S[N*NBAR];     // contains secret data
    byte randomness[2*BYTES + BYTES_SEED_A];      // contains secret data via randomness_s and randomness_seedSE
    byte *randomness_s = &randomness[0];                 // contains secret data
    byte *randomness_seedSE = &randomness[BYTES]; // contains secret data
    byte *randomness_z = &randomness[2*BYTES];
    byte shake_input_seedSE[1 + BYTES];           // contains secret data

    // Generate the secret value s, the seed for S and E, and the seed for the seed for A. Add seed_A to the public key
    RandomBytes(randomness, BYTES + BYTES + BYTES_SEED_A);
    

    //If Frodo640
    if (N == 640) {
        Shake128(pkSeedA, BYTES_SEED_A, randomness_z, BYTES_SEED_A);
    } else {
        //If Frodo976 or Frodo1344
        Shake256(pkSeedA, BYTES_SEED_A, randomness_z, BYTES_SEED_A);
    }
    

    // Generate S and E, and compute B = A*S + E. Generate A on-the-fly
    shake_input_seedSE[0] = 0x5F;
    std::memcpy(&shake_input_seedSE[1], randomness_seedSE, BYTES);

    //If Frodo640
    if (N == 640) {
        Shake128((byte*)S, 2*N*NBAR*sizeof(word16), shake_input_seedSE, 1 + BYTES);
    } else {
        //If Frodo976 or Frodo1344
        Shake256((byte*)S, 2*N*NBAR*sizeof(word16), shake_input_seedSE, 1 + BYTES);
    }
    for (size_t i = 0; i < 2 * N * NBAR; i++) {
        S[i] = LE_TO_UINT16(S[i]);
    }

    SampleN(S, N*NBAR, N);
    SampleN(E, N*NBAR, N);
    FrodoMulAddAsPlusE(B, S, E, pk);

    // Encode the second part of the public key
    Pack(pk_b, PUBLICKEYBYTES - BYTES_SEED_A, B, N*NBAR, LOGQ);

    // Add s, pk and S to the secret key
    std::memcpy(sk_s, randomness_s, BYTES);
    std::memcpy(sk_pk, pk, PUBLICKEYBYTES);
    for (size_t i = 0; i < N * NBAR; i++) {
        S[i] = UINT16_TO_LE(S[i]);
    }
    std::memcpy(sk_S, S, 2*N*NBAR);

    // Add H(pk) to the secret key
    //If Frodo640
    if (N == 640) {
        Shake128(sk_pkh, BYTES_PKHASH, pk, PUBLICKEYBYTES);
    } else {
        //If Frodo976 or Frodo1344
        Shake256(sk_pkh, BYTES_PKHASH, pk, PUBLICKEYBYTES);
    }

    // Cleanup:
    ClearBytes((byte *)S, N*NBAR*sizeof(word16));
    ClearBytes((byte *)E, N*NBAR*sizeof(word16));
    ClearBytes(randomness, 2*BYTES);
    ClearBytes(shake_input_seedSE, 1 + BYTES);
    return 0;
}

template<word16 T_SK_BYTES, word16 T_PK_BYTES, byte T_BYTES, word16 T_CT_BYTES, word16 T_N, byte T_LOGQ, byte T_EB>
int FrodoKEM<T_SK_BYTES, T_PK_BYTES, T_BYTES, T_CT_BYTES, T_N, T_LOGQ, T_EB>::KemEnc(byte *ct, byte *ss, const byte *pk)
{ // FrodoKEM's key encapsulation
    const byte *pkSeedA = &pk[0];
    const byte *pk_b = &pk[BYTES_SEED_A];
    byte *ct_c1 = &ct[0];
    byte *ct_c2 = &ct[(LOGQ*N*NBAR)/8];
    word16 B[N*NBAR] = {0};
    word16 V[NBAR*NBAR]= {0};                 // contains secret data
    word16 C[NBAR*NBAR] = {0};
    ALIGN_HEADER(32) word16 Bp[N*NBAR] ALIGN_FOOTER(32) = {0};
    ALIGN_HEADER(32) word16 Sp[(2*N+NBAR)*NBAR] ALIGN_FOOTER(32) = {0};  // contains secret data
    word16 *Ep = (word16 *)&Sp[N*NBAR];     // contains secret data
    word16 *Epp = (word16 *)&Sp[2*N*NBAR];  // contains secret data
    byte G2in[BYTES_PKHASH + BYTES_MU];                    // contains secret data via mu
    byte *pkh = &G2in[0];
    byte *mu = &G2in[BYTES_PKHASH];                        // contains secret data
    byte G2out[2*BYTES];                            // contains secret data
    byte *seedSE = &G2out[0];                              // contains secret data
    byte *k = &G2out[BYTES];                        // contains secret data
    byte Fin[CIPHERTEXTBYTES + BYTES];       // contains secret data via Fin_k
    byte *Fin_ct = &Fin[0];
    byte *Fin_k = &Fin[CIPHERTEXTBYTES];            // contains secret data
    byte shakeInputSeedSE[1 + BYTES];             // contains secret data

    // pkh <- G_1(pk), generate random mu, compute (seedSE || k) = G_2(pkh || mu)

    //If Frodo640
    if (N == 640) {
        Shake128(pkh, BYTES_PKHASH, pk, PUBLICKEYBYTES);
    } else {
        //If Frodo976 or Frodo1344
        Shake256(pkh, BYTES_PKHASH, pk, PUBLICKEYBYTES);
    }
    RandomBytes(mu, BYTES_MU);

    //If Frodo640
    if (N == 640) {
        Shake128(G2out, BYTES + BYTES, G2in, BYTES_PKHASH + BYTES_MU);
    } else {
        //If Frodo976 or Frodo1344
        Shake256(G2out, BYTES + BYTES, G2in, BYTES_PKHASH + BYTES_MU);
    }

    // Generate Sp and Ep, and compute Bp = Sp*A + Ep. Generate A on-the-fly
    shakeInputSeedSE[0] = 0x96;
    std::memcpy(&shakeInputSeedSE[1], seedSE, BYTES);
    //If Frodo640
    if (N == 640) {
        Shake128((byte*)Sp, (2*N+NBAR)*NBAR*sizeof(word16), shakeInputSeedSE, 1 + BYTES);
    } else {
        //If Frodo976 or Frodo1344
        Shake256((byte*)Sp, (2*N+NBAR)*NBAR*sizeof(word16), shakeInputSeedSE, 1 + BYTES);
    }

    for (size_t i = 0; i < (2 * N + NBAR) * NBAR; i++) {
        Sp[i] = LE_TO_UINT16(Sp[i]);
    }
    SampleN(Sp, N*NBAR, N);
    SampleN(Ep, N*NBAR, N);
    FrodoMulAddsAPlusE(Bp, Sp, Ep, pkSeedA);
    Pack(ct_c1, (LOGQ*N*NBAR)/8, Bp, N*NBAR, LOGQ);

    // Generate Epp, and compute V = Sp*B + Epp
    SampleN(Epp, NBAR*NBAR, N);
    Unpack(B, N*NBAR, pk_b, PUBLICKEYBYTES - BYTES_SEED_A, LOGQ);
    FrodoMulAddSbPlusE(V, B, Sp, Epp);

    // Encode mu, and compute C = V + enc(mu) (mod q)
    FrodoKeyEncode(C, (word16*)mu);
    FrodoAdd(C, V, C);
    Pack(ct_c2, (LOGQ*NBAR*NBAR)/8, C, NBAR*NBAR, LOGQ);

    // Compute ss = F(ct||KK)
    std::memcpy(Fin_ct, ct, CIPHERTEXTBYTES);
    std::memcpy(Fin_k, k, BYTES);
    //If Frodo640
    if (N == 640) {
        Shake128(ss, BYTES, Fin, CIPHERTEXTBYTES + BYTES);
    } else {
        //If Frodo976 or Frodo1344
        Shake256(ss, BYTES, Fin, CIPHERTEXTBYTES + BYTES);
    }

    // Cleanup:
    ClearBytes((byte *)V, NBAR*NBAR*sizeof(word16));
    ClearBytes((byte *)Sp, N*NBAR*sizeof(word16));
    ClearBytes((byte *)Ep, N*NBAR*sizeof(word16));
    ClearBytes((byte *)Epp, NBAR*NBAR*sizeof(word16));
    ClearBytes(mu, BYTES_MU);
    ClearBytes(G2out, 2*BYTES);
    ClearBytes(Fin_k, BYTES);
    ClearBytes(shakeInputSeedSE, 1 + BYTES);
    return 0;
}


template<word16 T_SK_BYTES, word16 T_PK_BYTES, byte T_BYTES, word16 T_CT_BYTES, word16 T_N, byte T_LOGQ, byte T_EB>
int FrodoKEM<T_SK_BYTES, T_PK_BYTES, T_BYTES, T_CT_BYTES, T_N, T_LOGQ, T_EB>::KemDec(byte *ss, const byte *ct, const byte *sk)
{ // FrodoKEM's key decapsulation
    word16 B[N*NBAR] = {0};
    word16 Bp[N*NBAR] = {0};
    word16 W[NBAR*NBAR] = {0};                // contains secret data
    word16 C[NBAR*NBAR] = {0};
    word16 CC[NBAR*NBAR] = {0};
    ALIGN_HEADER(32) word16 BBp[N*NBAR] ALIGN_FOOTER(32) = {0};
    ALIGN_HEADER(32) word16 Sp[(2*N+NBAR)*NBAR] ALIGN_FOOTER(32) = {0};  // contains secret data
    word16 *Ep = (word16 *)&Sp[N*NBAR];     // contains secret data
    word16 *Epp = (word16 *)&Sp[2*N*NBAR];  // contains secret data
    const byte *ct_c1 = &ct[0];
    const byte *ct_c2 = &ct[(LOGQ*N*NBAR)/8];
    const byte *sk_s = &sk[0];
    const byte *sk_pk = &sk[BYTES];
    const word16 *sk_S = (word16 *) &sk[BYTES + PUBLICKEYBYTES];
    word16 S[N * NBAR];                      // contains secret data
    const byte *sk_pkh = &sk[BYTES + PUBLICKEYBYTES + 2*N*NBAR];
    const byte *pk_seedA = &sk_pk[0];
    const byte *pk_b = &sk_pk[BYTES_SEED_A];
    byte G2in[BYTES_PKHASH + BYTES_MU];                   // contains secret data via muprime
    byte *pkh = &G2in[0];
    byte *muprime = &G2in[BYTES_PKHASH];                  // contains secret data
    byte G2out[2*BYTES];                           // contains secret data
    byte *seedSEprime = &G2out[0];                        // contains secret data
    byte *kprime = &G2out[BYTES];                  // contains secret data
    byte Fin[CIPHERTEXTBYTES + BYTES];      // contains secret data via Fin_k
    byte *Fin_ct = &Fin[0];
    byte *Fin_k = &Fin[CIPHERTEXTBYTES];           // contains secret data
    byte shakeInputSeedSEprime[1 + BYTES];       // contains secret data

    for (size_t i = 0; i < N * NBAR; i++) {
        S[i] = LE_TO_UINT16(sk_S[i]);
    }

    // Compute W = C - Bp*S (mod q), and decode the randomness mu
    Unpack(Bp, N*NBAR, ct_c1, (LOGQ*N*NBAR)/8, LOGQ);
    Unpack(C, NBAR*NBAR, ct_c2, (LOGQ*NBAR*NBAR)/8, LOGQ);
    FrodoMulBs(W, Bp, S);
    FrodoSub(W, C, W);
    FrodoKeyDecode((word16*)muprime, W);

    // Generate (seedSE' || k') = G_2(pkh || mu')
    std::memcpy(pkh, sk_pkh, BYTES_PKHASH);
    //If Frodo640
    if (N == 640) {
        Shake128(G2out, BYTES + BYTES, G2in, BYTES_PKHASH + BYTES_MU);
    } else {
        //If Frodo976 or Frodo1344
        Shake256(G2out, BYTES + BYTES, G2in, BYTES_PKHASH + BYTES_MU);
    }

    // Generate Sp and Ep, and compute BBp = Sp*A + Ep. Generate A on-the-fly
    shakeInputSeedSEprime[0] = 0x96;
    std::memcpy(&shakeInputSeedSEprime[1], seedSEprime, BYTES);
    //If Frodo640
    if (N == 640) {
        Shake128((byte*)Sp, (2*N+NBAR)*NBAR*sizeof(word16), shakeInputSeedSEprime, 1 + BYTES);
    } else {
        //If Frodo976 or Frodo1344
        Shake256((byte*)Sp, (2*N+NBAR)*NBAR*sizeof(word16), shakeInputSeedSEprime, 1 + BYTES);
    }

    for (size_t i = 0; i < (2*N+NBAR)*NBAR; i++) {
        Sp[i] = LE_TO_UINT16(Sp[i]);
    }
    SampleN(Sp, N*NBAR, N);
    SampleN(Ep, N*NBAR, N);
    FrodoMulAddsAPlusE(BBp, Sp, Ep, pk_seedA);

    // Generate Epp, and compute W = Sp*B + Epp
    SampleN(Epp, NBAR*NBAR, N);
    Unpack(B, N*NBAR, pk_b, PUBLICKEYBYTES - BYTES_SEED_A, LOGQ);
    FrodoMulAddSbPlusE(W, B, Sp, Epp);

    // Encode mu, and compute CC = W + enc(mu') (mod q)
    FrodoKeyEncode(CC, (word16*)muprime);
    FrodoAdd(CC, W, CC);

    // Prepare input to F
    std::memcpy(Fin_ct, ct, CIPHERTEXTBYTES);

    // Reducing BBp modulo q
    for (int i = 0; i < N*NBAR; i++) BBp[i] = BBp[i] & ((1 << LOGQ)-1);

    // If (Bp == BBp & C == CC) then ss = F(ct || k'), else ss = F(ct || s)
    // Needs to avoid branching on secret data as per:
    //     Qian Guo, Thomas Johansson, Alexander Nilsson. A key-recovery timing attack on post-quantum 
    //     primitives using the Fujisaki-Okamoto transformation and its application on FrodoKEM. In CRYPTO 2020.
    sbyte selector = CtVerify(Bp, BBp, N*NBAR) | CtVerify(C, CC, NBAR*NBAR);
    // If (selector == 0) then load k' to do ss = F(ct || k'), else if (selector == -1) load s to do ss = F(ct || s)
    CtSelect((byte*)Fin_k, (byte*)kprime, (byte*)sk_s, BYTES, selector);

    //If Frodo640
    if (N == 640) {
        Shake128(ss, BYTES, Fin, CIPHERTEXTBYTES + BYTES);
    } else {
        //If Frodo976 or Frodo1344
        Shake256(ss, BYTES, Fin, CIPHERTEXTBYTES + BYTES);
    }

    // Cleanup:
    ClearBytes((byte *)W, NBAR*NBAR*sizeof(word16));
    ClearBytes((byte *)Sp, N*NBAR*sizeof(word16));
    ClearBytes((byte *)S, N*NBAR*sizeof(word16));
    ClearBytes((byte *)Ep, N*NBAR*sizeof(word16));
    ClearBytes((byte *)Epp, NBAR*NBAR*sizeof(word16));
    ClearBytes(muprime, BYTES_MU);
    ClearBytes(G2out, 2*BYTES);
    ClearBytes(Fin_k, BYTES);
    ClearBytes(shakeInputSeedSEprime, 1 + BYTES);
    return 0;
}


template<word16 T_SK_BYTES, word16 T_PK_BYTES, byte T_BYTES, word16 T_CT_BYTES, word16 T_N, byte T_LOGQ, byte T_EB>
int FrodoKEM<T_SK_BYTES, T_PK_BYTES, T_BYTES, T_CT_BYTES, T_N, T_LOGQ, T_EB>::FrodoMulAddAsPlusE(word16 *out, const word16 *s, const word16 *e, const byte *seedA) 
{ // Generate-and-multiply: generate matrix A (N x N) row-wise, multiply by s on the right.
  // Inputs: s, e (N x N_BAR)
  // Output: out = A*s + e (N x N_BAR)
    int i, j, k;
    sword16 A[N * N] = {0};       
    // Matrix A generation using SHAKE128, done per 16*N-bit row   
    byte seedASeparated[2 + BYTES_SEED_A];
    word16* seedAOrigin = (word16*)&seedASeparated;
    std::memcpy(&seedASeparated[2], seedA, BYTES_SEED_A);
    for (i = 0; i < N; i++) {
        seedAOrigin[0] = UINT16_TO_LE((word16) i);
        Shake128((byte*)(A + i*N), (unsigned long long)(2*N), seedASeparated, 2 + BYTES_SEED_A);
    }
    
    for (i = 0; i < N * N; i++) {
        A[i] = LE_TO_UINT16(A[i]);
    }
    std::memcpy(out, e, NBAR * N * sizeof(word16));  

    for (i = 0; i < N; i++) {                            // Matrix multiplication-addition A*s + e
        for (k = 0; k < NBAR; k++) {
            word16 sum = 0;
            for (j = 0; j < N; j++) {                                
                sum += A[i*N + j] * s[k*N + j];  
            }
            out[i*NBAR + k] += sum;                      // Adding e. No need to reduce modulo 2^15, extra bits are taken care of during packing later on.
        }
    }
    
    return 1;
}

template<word16 T_SK_BYTES, word16 T_PK_BYTES, byte T_BYTES, word16 T_CT_BYTES, word16 T_N, byte T_LOGQ, byte T_EB>
int FrodoKEM<T_SK_BYTES, T_PK_BYTES, T_BYTES, T_CT_BYTES, T_N, T_LOGQ, T_EB>::FrodoMulAddsAPlusE(word16 *out, const word16 *s, word16 *e, const byte *seedA) 
{ // Generate-and-multiply: generate matrix A (N x N) column-wise, multiply by s' on the left.
  // Inputs: s', e' (N_BAR x N)
  // Output: out = s'*A + e' (N_BAR x N)
    int i, j, k;
    int16_t A[N * N] = {0};        
  
    // Matrix A generation using SHAKE128, done per 16*N-bit row
    byte seedASeparated[2 + BYTES_SEED_A];
    word16* seedAOrigin = (word16*)&seedASeparated;
    std::memcpy(&seedASeparated[2], seedA, BYTES_SEED_A);
    for (i = 0; i < N; i++) {
        seedAOrigin[0] = UINT16_TO_LE((word16) i);
        Shake128((byte*)(A + i*N), (unsigned long long)(2*N), seedASeparated, 2 + BYTES_SEED_A);
    }

    for (i = 0; i < N * N; i++) {
        A[i] = LE_TO_UINT16(A[i]);
    }
    std::memcpy(out, e, NBAR * N * sizeof(word16));

    for (i = 0; i < N; i++) {                            // Matrix multiplication-addition A*s + e
        for (k = 0; k < NBAR; k++) {
            word16 sum = 0;
            for (j = 0; j < N; j++) {                                
                sum += A[j*N + i] * s[k*N + j];  
            }
            out[k*N + i] += sum;                         // Adding e. No need to reduce modulo 2^15, extra bits are taken care of during packing later on.
        }
    }
    
    return 1;
}

template<word16 T_SK_BYTES, word16 T_PK_BYTES, byte T_BYTES, word16 T_CT_BYTES, word16 T_N, byte T_LOGQ, byte T_EB>
void FrodoKEM<T_SK_BYTES, T_PK_BYTES, T_BYTES, T_CT_BYTES, T_N, T_LOGQ, T_EB>::FrodoMulBs(word16 *out, const word16 *b, const word16 *s) 
{ // Multiply by s on the right
  // Inputs: b (N_BAR x N), s (N x N_BAR)
  // Output: out = b*s (N_BAR x N_BAR)
    int i, j, k;

    for (i = 0; i < NBAR; i++) {
        for (j = 0; j < NBAR; j++) {
            out[i*NBAR + j] = 0;
            for (k = 0; k < N; k++) {
                out[i*NBAR + j] += b[i*N + k] * (int16_t)s[j*N + k];
            }
            out[i*NBAR + j] = (word32)(out[i*NBAR + j]) & ((1<<LOGQ)-1);
        }
    }
}

template<word16 T_SK_BYTES, word16 T_PK_BYTES, byte T_BYTES, word16 T_CT_BYTES, word16 T_N, byte T_LOGQ, byte T_EB>
void FrodoKEM<T_SK_BYTES, T_PK_BYTES, T_BYTES, T_CT_BYTES, T_N, T_LOGQ, T_EB>::FrodoMulAddSbPlusE(word16 *out, const word16 *b, const word16 *s, const word16 *e) 
{ // Multiply by s on the left
  // Inputs: b (N x N_BAR), s (N_BAR x N), e (N_BAR x N_BAR)
  // Output: out = s*b + e (N_BAR x N_BAR)
    int i, j, k;

    for (k = 0; k < NBAR; k++) {
        for (i = 0; i < NBAR; i++) {
            out[k*NBAR + i] = e[k*NBAR + i];
            for (j = 0; j < N; j++) {
                out[k*NBAR + i] += (sword16)s[k*N + j] * b[j*NBAR + i];
            }
            out[k*NBAR + i] = (word32)(out[k*NBAR + i]) & ((1<<LOGQ)-1);
        }
    }
}

template<word16 T_SK_BYTES, word16 T_PK_BYTES, byte T_BYTES, word16 T_CT_BYTES, word16 T_N, byte T_LOGQ, byte T_EB>
void FrodoKEM<T_SK_BYTES, T_PK_BYTES, T_BYTES, T_CT_BYTES, T_N, T_LOGQ, T_EB>::FrodoAdd(word16 *out, const word16 *a, const word16 *b) 
{ // Add a and b
  // Inputs: a, b (N_BAR x N_BAR)
  // Output: c = a + b

    for (int i = 0; i < (NBAR*NBAR); i++) {
        out[i] = (a[i] + b[i]) & ((1<<LOGQ)-1);
    }
}

template<word16 T_SK_BYTES, word16 T_PK_BYTES, byte T_BYTES, word16 T_CT_BYTES, word16 T_N, byte T_LOGQ, byte T_EB>
void FrodoKEM<T_SK_BYTES, T_PK_BYTES, T_BYTES, T_CT_BYTES, T_N, T_LOGQ, T_EB>::FrodoSub(word16 *out, const word16 *a, const word16 *b) 
{ // Subtract a and b
  // Inputs: a, b (N_BAR x N_BAR)
  // Output: c = a - b

    for (int i = 0; i < (NBAR*NBAR); i++) {
        out[i] = (a[i] - b[i]) & ((1<<LOGQ)-1);
    }
}

template<word16 T_SK_BYTES, word16 T_PK_BYTES, byte T_BYTES, word16 T_CT_BYTES, word16 T_N, byte T_LOGQ, byte T_EB>
void FrodoKEM<T_SK_BYTES, T_PK_BYTES, T_BYTES, T_CT_BYTES, T_N, T_LOGQ, T_EB>::FrodoKeyEncode(word16 *out, const word16 *in) 
{ // Encoding
    unsigned int i, j, nPiecesWord = 8;
    unsigned int nAmountWords = (NBAR*NBAR)/8;
    word64 temp, mask = ((word64)1 << EXTRACTED_BITS) - 1;
    word16* pos = out;

    for (i = 0; i < nAmountWords; i++) {
        temp = 0;
        for(j = 0; j < EXTRACTED_BITS; j++) 
            temp |= ((word64)((byte*)in)[i*EXTRACTED_BITS + j]) << (8*j);
        for (j = 0; j < nPiecesWord; j++) { 
            *pos = (word16)((temp & mask) << (LOGQ - EXTRACTED_BITS));  
            temp >>= EXTRACTED_BITS;
            pos++;
        }
    }
}

template<word16 T_SK_BYTES, word16 T_PK_BYTES, byte T_BYTES, word16 T_CT_BYTES, word16 T_N, byte T_LOGQ, byte T_EB>
void FrodoKEM<T_SK_BYTES, T_PK_BYTES, T_BYTES, T_CT_BYTES, T_N, T_LOGQ, T_EB>::FrodoKeyDecode(word16 *out, const word16 *in)
{ // Decoding
    unsigned int i, j, index = 0, nPiecesWord = 8;
    unsigned int nAmountWords = (NBAR * NBAR) / 8;
    word16 temp, maskex=((word16)1 << EXTRACTED_BITS) -1, maskq =((word16)1 << LOGQ) -1;
    byte  *pos = (byte*)out;
    word64 templong;

    for (i = 0; i < nAmountWords; i++) {
        templong = 0;
        for (j = 0; j < nPiecesWord; j++) {  // temp = floor(in*2^{-11}+0.5)
            temp = ((in[index] & maskq) + (1 << (LOGQ - EXTRACTED_BITS - 1))) >> (LOGQ - EXTRACTED_BITS);
            templong |= ((word64)(temp & maskex)) << (EXTRACTED_BITS * j);
            index++;
        }
	for(j = 0; j < EXTRACTED_BITS; j++) 
	    pos[i*EXTRACTED_BITS + j] = (templong >> (8*j)) & 0xFF;
    }
}

template<word16 T_SK_BYTES, word16 T_PK_BYTES, byte T_BYTES, word16 T_CT_BYTES, word16 T_N, byte T_LOGQ, byte T_EB>
void FrodoKEM<T_SK_BYTES, T_PK_BYTES, T_BYTES, T_CT_BYTES, T_N, T_LOGQ, T_EB>::SampleN(word16 *s, const size_t n, word16 version) 
{ // Fills vector s with n samples from the noise distribution which requires 16 bits to sample. 
  // The distribution is specified by its CDF.
  // Input: pseudo-random values (2*n bytes) passed in s. The input is overwritten by the output.
    word32 i, j;
    //Variables depend on Frodo version
    std::vector<word16> cdfTable;
    word16 cdfTableLen;

    if (version == 640) {
        cdfTable = {4643, 13363, 20579, 25843, 29227, 31145, 32103, 32525, 32689, 32745, 32762, 32766, 32767};
        cdfTableLen = 13;
    } else if (version == 976) {
        cdfTable = {5638, 15915, 23689, 28571, 31116, 32217, 32613, 32731, 32760, 32766, 32767};
        cdfTableLen = 11;
    } else {
        cdfTable = {9142, 23462, 30338, 32361, 32725, 32765, 32767};
        cdfTableLen = 7;
    };

    for (i = 0; i < n; ++i) {
        word16 sample = 0;
        word16 prnd = s[i] >> 1;    // Drop the least significant bit
        word16 sign = s[i] & 0x1;    // Pick the least significant bit

        // No need to compare with the last value.
        for (j = 0; j < (word32)(cdfTableLen - 1); j++) {
            // Constant time comparison: 1 if CDF_TABLE[j] < s, 0 otherwise. Uses the fact that CDF_TABLE[j] and s fit in 15 bits.
            sample += (word16)(cdfTable.at(j) - prnd) >> 15;
        }
        // Assuming that sign is either 0 or 1, flips sample iff sign = 1
        s[i] = ((-sign) ^ sample) + sign;
    }
}

#define min(x, y) (((x) < (y)) ? (x) : (y))

template<word16 T_SK_BYTES, word16 T_PK_BYTES, byte T_BYTES, word16 T_CT_BYTES, word16 T_N, byte T_LOGQ, byte T_EB>
void FrodoKEM<T_SK_BYTES, T_PK_BYTES, T_BYTES, T_CT_BYTES, T_N, T_LOGQ, T_EB>::Pack(byte *out, const size_t outLen, const word16 *in, const size_t inLen, const byte lsb) 
{ // Pack the input uint16 vector into a char output vector, copying lsb bits from each input element. 
  // If inlen * lsb / 8 > outlen, only outlen * 8 bits are copied.
    std::memset(out, 0, outLen);

    size_t i = 0;            // whole bytes already filled in
    size_t j = 0;            // whole uint16_t already copied
    word16 w = 0;          // the leftover, not yet copied
    byte bits = 0;  // the number of lsb in w

    while (i < outLen && (j < inLen || ((j == inLen) && (bits > 0)))) {
        /*
        in: |        |        |********|********|
                              ^
                              j
        w : |   ****|
                ^
               bits
        out:|**|**|**|**|**|**|**|**|* |
                                    ^^
                                    ib
        */
        byte b = 0;  // bits in out[i] already filled in
        while (b < 8) {
            int nbits = min(8 - b, bits);
            word16 mask = (1 << nbits) - 1;
            byte t = (w >> (bits - nbits)) & mask;  // the bits to copy from w to out
            out[i] = out[i] + (t << (8 - b - nbits));
            b += nbits;
            bits -= nbits;
            w &= ~(mask << bits);  // not strictly necessary; mostly for debugging

            if (bits == 0) {
                if (j < inLen) {
                    w = in[j];
                    bits = lsb;
                    j++;
                } else {
                    break;  // the input vector is exhausted
                }
            }
        }
        if (b == 8) {  // out[i] is filled in
            i++;
        }
    }
}


template<word16 T_SK_BYTES, word16 T_PK_BYTES, byte T_BYTES, word16 T_CT_BYTES, word16 T_N, byte T_LOGQ, byte T_EB>
void FrodoKEM<T_SK_BYTES, T_PK_BYTES, T_BYTES, T_CT_BYTES, T_N, T_LOGQ, T_EB>::Unpack(word16 *out, const size_t outLen, const byte *in, const size_t inLen, const byte lsb) 
{ // Unpack the input char vector into a uint16_t output vector, copying lsb bits
  // for each output element from input. outlen must be at least ceil(inlen * 8 / lsb).
    std::memset(out, 0, outLen * sizeof(word16));

    size_t i = 0;            // whole uint16_t already filled in
    size_t j = 0;            // whole bytes already copied
    byte w = 0;     // the leftover, not yet copied
    byte bits = 0;  // the number of lsb bits of w

    while (i < outLen && (j < inLen || ((j == inLen) && (bits > 0)))) {
        /*
        in: |  |  |  |  |  |  |**|**|...
                              ^
                              j
        w : | *|
              ^
              bits
        out:|   *****|   *****|   ***  |        |...
                              ^   ^
                              i   b
        */
        byte b = 0;  // bits in out[i] already filled in
        while (b < lsb) {
            int nbits = min(lsb - b, bits);
            word16 mask = (1 << nbits) - 1;
            byte t = (w >> (bits - nbits)) & mask;  // the bits to copy from w to out
            out[i] = out[i] + (t << (lsb - b - nbits));
            b += nbits;
            bits -= nbits;
            w &= ~(mask << bits);  // not strictly necessary; mostly for debugging

            if (bits == 0) {
                if (j < inLen) {
                    w = in[j];
                    bits = 8;
                    j++;
                } else {
                    break;  // the input vector is exhausted
                }
            }
        }
        if (b == lsb) {  // out[i] is filled in
            i++;
        }
    }
}

template<word16 T_SK_BYTES, word16 T_PK_BYTES, byte T_BYTES, word16 T_CT_BYTES, word16 T_N, byte T_LOGQ, byte T_EB>
void FrodoKEM<T_SK_BYTES, T_PK_BYTES, T_BYTES, T_CT_BYTES, T_N, T_LOGQ, T_EB>::CtSelect(byte *r, const byte *a, const byte *b, size_t len, sbyte selector) 
{ // Select one of the two input arrays to be moved to r
  // If (selector == 0) then load r with a, else if (selector == -1) load r with b
    for (size_t i = 0; i < len; i++) {
        r[i] = (~selector & a[i]) | (selector & b[i]);
    }
}

template<word16 T_SK_BYTES, word16 T_PK_BYTES, byte T_BYTES, word16 T_CT_BYTES, word16 T_N, byte T_LOGQ, byte T_EB>
void FrodoKEM<T_SK_BYTES, T_PK_BYTES, T_BYTES, T_CT_BYTES, T_N, T_LOGQ, T_EB>::ClearBytes(byte *mem, size_t n)
{ // Clear 8-bit bytes from memory. "n" indicates the number of bytes to be zeroed.
  // This function uses the volatile type qualifier to inform the compiler not to optimize out the memory clearing.
    volatile byte *v = mem; 

    for (size_t i = 0; i < n; i++) {
        v[i] = 0;
    }
}

template<word16 T_SK_BYTES, word16 T_PK_BYTES, byte T_BYTES, word16 T_CT_BYTES, word16 T_N, byte T_LOGQ, byte T_EB>
void FrodoKEM<T_SK_BYTES, T_PK_BYTES, T_BYTES, T_CT_BYTES, T_N, T_LOGQ, T_EB>::RandomBytes(byte *out, size_t outLen) {
  AutoSeededRandomPool rng;
  rng.GenerateBlock(out, outLen);
}


template<word16 T_SK_BYTES, word16 T_PK_BYTES, byte T_BYTES, word16 T_CT_BYTES, word16 T_N, byte T_LOGQ, byte T_EB>
sbyte FrodoKEM<T_SK_BYTES, T_PK_BYTES, T_BYTES, T_CT_BYTES, T_N, T_LOGQ, T_EB>::CtVerify(const word16 *a, const word16 *b, size_t len) 
{ // Compare two arrays in constant time.
  // Returns 0 if the byte arrays are equal, -1 otherwise.
    word16 r = 0;

    for (size_t i = 0; i < len; i++) {
        r |= a[i] ^ b[i];
    }

    r = (-(sword16)(r >> 1) | -(sword16)(r & 1)) >> (8*sizeof(word16)-1);
    return (sbyte)r;
}




template class FrodoKEM<19888, 9616, 16, 9720, 640, 15, 2>;
template class FrodoKEM<31296, 15632, 24, 15744, 976, 16, 3>; 
template class FrodoKEM<43088, 21520, 32, 21632, 1344, 16, 4>;




NAMESPACE_END
