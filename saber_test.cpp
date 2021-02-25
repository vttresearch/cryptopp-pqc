/** The C++ file for the tests of the post-quantum KEM algorithm SABER. 
Adapted by Julius Hekkala from the reference
implementation of SABER by the SABER team 
(https://github.com/KULeuven-COSIC/SABER) */

#include "saber_test.h"
#include "saber.h"
#include "misc.h"
#include "osrng.h"
#include <iostream>

NAMESPACE_BEGIN(CryptoPP)
NAMESPACE_BEGIN(Test)

#define NTESTS 1000


/**
 * Test the LightSaber version
 * 
 */
int TestLightSaber() {
    std::cout << "Testing LightSaber..." << std::endl;
    int res = 0;

    LightSaber saber = LightSaber();

    for (int i=0; i < NTESTS; i++) {
        byte pk[LightSaber::PUBLICKEYBYTES];
        byte sk[LightSaber::SECRETKEYBYTES];
        byte ct[LightSaber::CIPHERTEXTBYTES];	
        byte ssA[LightSaber::SHAREDSECRETBYTES], ssB[LightSaber::SHAREDSECRETBYTES];
            
                    
        //Generate public key pk and secret key sk
        saber.KemKeypair(pk, sk);

        //Key encapsulation
        saber.KemEnc(ct, ssA, pk);

        
        //Key decapsulation
        saber.KemDec(ssB, ct, sk);
        if (memcmp(ssA, ssB, LightSaber::SHAREDSECRETBYTES)) {
            res = 1;
        }

    }

    if (res == 0) {
        std::cout << "LightSaber tests passed." << std::endl;
    } else {
        std::cout << "LightSaber tests failed." << std::endl;
    }
    return res;
    

}


/**
 * 
 * Test the standard Saber version
 */
int TestSaber() {
    std::cout << "Testing Saber..." << std::endl;
    int res = 0;

    Saber saber = Saber();

    for (int i=0; i < NTESTS; i++) {
        byte pk[Saber::PUBLICKEYBYTES];
        byte sk[Saber::SECRETKEYBYTES];
        byte ct[Saber::CIPHERTEXTBYTES];	
        byte ssA[Saber::SHAREDSECRETBYTES], ssB[Saber::SHAREDSECRETBYTES];
            
                    
        //Generate public key pk and secret key sk
        saber.KemKeypair(pk, sk);

        //Key encapsulation
        saber.KemEnc(ct, ssA, pk);

        
        //Key decapsulation
        saber.KemDec(ssB, ct, sk);
        if (memcmp(ssA, ssB, Saber::SHAREDSECRETBYTES)) {
            res = 1;
        }

    }

    if (res == 0) {
        std::cout << "Saber tests passed." << std::endl;
    } else {
        std::cout << "Saber tests failed." << std::endl;
    }
    return res;
    

}


/**
 * Test FireSaber version
 * 
 */
int TestFireSaber() {
    std::cout << "Testing FireSaber..." << std::endl;
    int res = 0;

    FireSaber saber = FireSaber();

    for (int i=0; i < NTESTS; i++) {
        byte pk[FireSaber::PUBLICKEYBYTES];
        byte sk[FireSaber::SECRETKEYBYTES];
        byte ct[FireSaber::CIPHERTEXTBYTES];	
        byte ssA[FireSaber::SHAREDSECRETBYTES], ssB[FireSaber::SHAREDSECRETBYTES];
            
                    
        //Generate public key pk and secret key sk
        saber.KemKeypair(pk, sk);

        //Key encapsulation
        saber.KemEnc(ct, ssA, pk);

        
        //Key decapsulation
        saber.KemDec(ssB, ct, sk);
        if (memcmp(ssA, ssB, FireSaber::SHAREDSECRETBYTES)) {
            res = 1;
        }

    }

    if (res == 0) {
        std::cout << "FireSaber tests passed." << std::endl;
    } else {
        std::cout << "FireSaber tests failed." << std::endl;
    }
    return res;
    

}

int TestInvalidSkSaber()
{
  std::cout << "Testing invalid secret keys" << std::endl;
  unsigned int i;
  byte pk[Saber::PUBLICKEYBYTES];
  byte sk[Saber::SECRETKEYBYTES];
  byte ct[Saber::CIPHERTEXTBYTES];
  byte key_a[Saber::SHAREDSECRETBYTES];
  byte key_b[Saber::SHAREDSECRETBYTES];

  Saber saber = Saber();
  bool pass = true;

  for(i=0;i<NTESTS;i++) {
    //Alice generates a public key
    saber.KemKeypair(pk, sk);
    //Bob derives a secret key and creates a response
    saber.KemEnc(ct, key_b, pk);

    //Replace secret key with random values
    AutoSeededRandomPool rng;
    rng.GenerateBlock(sk, Saber::SECRETKEYBYTES);

    //Alice uses Bobs response to get her shared key
    saber.KemDec(key_a, ct, sk);

    if(!memcmp(key_a, key_b, Saber::SHAREDSECRETBYTES)) {
      std::cout << "ERROR invalid sk" << std::endl;
      pass = false;
    }
      
  }
  if (pass) {
    return 0;
  } else {
    return 1;
  }
  
}

static int TestInvalidCiphertext()
{

  unsigned int i;
  byte pk[Saber::PUBLICKEYBYTES];
  byte sk[Saber::SECRETKEYBYTES];
  byte ct[Saber::CIPHERTEXTBYTES];
  byte key_a[Saber::SHAREDSECRETBYTES];
  byte key_b[Saber::SHAREDSECRETBYTES];
  size_t pos;

  Saber saber = Saber();

  std::cout << "Testing invalid ciphertext..." << std::endl;


  bool pass = true;

  for(i=0;i<NTESTS;i++) {
    AutoSeededRandomPool rng;
    rng.GenerateBlock((byte *) &pos, sizeof(size_t));

    //Alice generates a public key
    saber.KemKeypair(pk, sk);

    //Bob derives a secret key and creates a response
    saber.KemEnc(ct, key_b, pk);

    //Change some byte in the ciphertext (i.e., encapsulated key)
    ct[pos % Saber::CIPHERTEXTBYTES] ^= 23;

    //Alice uses Bobs response to get her shared key
    saber.KemDec(key_a, ct, sk);

    if(!memcmp(key_a, key_b, Saber::SHAREDSECRETBYTES)) {
      std::cout << "ERROR invalid ciphertext\n" << std::endl;
      pass = false;
    }
      
  }
  if (pass) {
    return 0;
  } else {
    return 1;
  }
  
}



int RunSaberTests() {
    std::cout << "Testing that different versions of SABER work correctly" << std::endl;

    int ret1 = TestLightSaber();
    int ret2 = TestSaber();
    int ret3 = TestFireSaber();
    int ret4 = TestInvalidSkSaber();
    int ret5 = TestInvalidCiphertext();

    int res = ret1 + ret2 + ret3 + ret4 + ret5;
    if (res > 0) {
        std::cout << "SABER tests failed." << std::endl;
    } else {
        std::cout << "SABER tests completed succesfully." << std::endl;
    }

    return res;
}




NAMESPACE_END
NAMESPACE_END