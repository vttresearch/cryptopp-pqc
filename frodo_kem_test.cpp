/**
 * The C++ file for FrodoKEM tests. Adapted
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
 * 
 */

#include <iostream>
#include "frodo_kem_test.h"
#include "frodo_kem.h"

NAMESPACE_BEGIN(CryptoPP)
NAMESPACE_BEGIN(Test)


#define NTESTS 100

/**
 * Run Frodo640 tests
 * 
 * @return int 
 */
int TestFrodo640() 
{
    byte pk[Frodo640::PUBLICKEYBYTES];
    byte sk[Frodo640::SECRETKEYBYTES];
    byte ss_encap[Frodo640::BYTES], ss_decap[Frodo640::BYTES];
    byte ct[Frodo640::CIPHERTEXTBYTES];
    byte bytes[4];
    word32* pos = (word32*)bytes;
    byte Fin[Frodo640::CIPHERTEXTBYTES + Frodo640::BYTES];
    Frodo640 frodo = Frodo640();
    std::cout << "Running tests for FrodoKEM-640" << std::endl;

    for (int i = 0; i < NTESTS; i++) {
        frodo.KemKeypair(pk, sk);
        frodo.KemEnc(ct, ss_encap, pk);
        frodo.KemDec(ss_decap, ct, sk);

        if (!VerifyBufsEqual(ss_encap, ss_decap, Frodo640::BYTES)) {
            std::cout << "ERROR -- encapsulation/decapsulation mechanism failed!" << std::endl;
	        return false; 
        }

        // Testing decapsulation after changing random bits of a random 16-bit digit of ct
        frodo.RandomBytes(bytes, 4);
        *pos %= Frodo640::CIPHERTEXTBYTES/2;
        if (*pos == 0) {
            *pos = 1;
        }

        ((word16*)ct)[*pos] ^= *pos;
  
        frodo.KemDec(ss_decap, ct, sk);


        // Compute ss = F(ct || s) with modified ct
        std::memcpy(Fin, ct, Frodo640::CIPHERTEXTBYTES);
        std::memcpy(&Fin[Frodo640::CIPHERTEXTBYTES], sk, Frodo640::BYTES);
        frodo.Shake128(ss_encap, Frodo640::BYTES, Fin, Frodo640::CIPHERTEXTBYTES + Frodo640::BYTES);

        
        if (!VerifyBufsEqual(ss_encap, ss_decap, Frodo640::BYTES)) {
            std::cout << "ERROR -- changing random bits of the ciphertext should cause a failure!" << std::endl;
	          return false;
        }
    }
    std::cout << "Tests PASSED. All session keys matched." << std::endl;

    return true;
}

/**
 * Run Frodo976 tests
 * 
 * @return int 
 */
int TestFrodo976() 
{
    byte pk[Frodo976::PUBLICKEYBYTES];
    byte sk[Frodo976::SECRETKEYBYTES];
    byte ss_encap[Frodo976::BYTES], ss_decap[Frodo976::BYTES];
    byte ct[Frodo976::CIPHERTEXTBYTES];
    byte bytes[4];
    word32* pos = (word32*)bytes;
    byte Fin[Frodo976::CIPHERTEXTBYTES + Frodo976::BYTES];
    Frodo976 frodo = Frodo976();
    std::cout << "Running tests for FrodoKEM-976" << std::endl;

    for (int i = 0; i < NTESTS; i++) {
        frodo.KemKeypair(pk, sk);
        frodo.KemEnc(ct, ss_encap, pk);
        frodo.KemDec(ss_decap, ct, sk);

        if (!VerifyBufsEqual(ss_encap, ss_decap, Frodo976::BYTES)) {
            std::cout << "ERROR -- encapsulation/decapsulation mechanism failed!" << std::endl;
	        return false; 
        }

        // Testing decapsulation after changing random bits of a random 16-bit digit of ct
        frodo.RandomBytes(bytes, 4);
        *pos %= Frodo976::CIPHERTEXTBYTES/2;
        if (*pos == 0) {
            *pos = 1;
        }

        ((word16*)ct)[*pos] ^= *pos;
  
        frodo.KemDec(ss_decap, ct, sk);


        // Compute ss = F(ct || s) with modified ct
        std::memcpy(Fin, ct, Frodo976::CIPHERTEXTBYTES);
        std::memcpy(&Fin[Frodo976::CIPHERTEXTBYTES], sk, Frodo976::BYTES);
        frodo.Shake256(ss_encap, Frodo976::BYTES, Fin, Frodo976::CIPHERTEXTBYTES + Frodo976::BYTES);

        
        if (!VerifyBufsEqual(ss_encap, ss_decap, Frodo976::BYTES)) {
            std::cout << "ERROR -- changing random bits of the ciphertext should cause a failure!" << std::endl;
	          return false;
        }
    }
    std::cout << "Tests PASSED. All session keys matched." << std::endl;

    return true;
}

/**
 * Run Frodo1344 tests
 * 
 * @return int 
 */
int TestFrodo1344() 
{
    byte pk[Frodo1344::PUBLICKEYBYTES];
    byte sk[Frodo1344::SECRETKEYBYTES];
    byte ss_encap[Frodo1344::BYTES], ss_decap[Frodo1344::BYTES];
    byte ct[Frodo1344::CIPHERTEXTBYTES];
    byte bytes[4];
    word32* pos = (word32*)bytes;
    byte Fin[Frodo1344::CIPHERTEXTBYTES + Frodo1344::BYTES];
    Frodo1344 frodo = Frodo1344();
    std::cout << "Running tests for FrodoKEM-1344" << std::endl;

    for (int i = 0; i < NTESTS; i++) {
        frodo.KemKeypair(pk, sk);
        frodo.KemEnc(ct, ss_encap, pk);
        frodo.KemDec(ss_decap, ct, sk);

        if (!VerifyBufsEqual(ss_encap, ss_decap, Frodo1344::BYTES)) {
            std::cout << "ERROR -- encapsulation/decapsulation mechanism failed!" << std::endl;
	        return false; 
        }

        // Testing decapsulation after changing random bits of a random 16-bit digit of ct
        frodo.RandomBytes(bytes, 4);
        *pos %= Frodo1344::CIPHERTEXTBYTES/2;
        if (*pos == 0) {
            *pos = 1;
        }

        ((word16*)ct)[*pos] ^= *pos;
  
        frodo.KemDec(ss_decap, ct, sk);


        // Compute ss = F(ct || s) with modified ct
        std::memcpy(Fin, ct, Frodo1344::CIPHERTEXTBYTES);
        std::memcpy(&Fin[Frodo1344::CIPHERTEXTBYTES], sk, Frodo1344::BYTES);
        frodo.Shake256(ss_encap, Frodo1344::BYTES, Fin, Frodo1344::CIPHERTEXTBYTES + Frodo1344::BYTES);

        if (!VerifyBufsEqual(ss_encap, ss_decap, Frodo1344::BYTES)) {
            std::cout << "ERROR -- changing random bits of the ciphertext should cause a failure!" << std::endl;
	          return false;
        }
    }
    std::cout << "Tests PASSED. All session keys matched." << std::endl;

    return true;
}


/**
 * Run FrodoKEM tests for all parameter sets.
 */
void RunFrodoKemTests() {
    TestFrodo640();
    TestFrodo976();
    TestFrodo1344();
}

NAMESPACE_END
NAMESPACE_END

