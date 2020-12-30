/** The C++ file including the tests for CRYSTALS-Kyber. 
Adapted from the reference
implementation of kyber by the CRYSTALS team 
(https://github.com/pq-crystals/kyber) */


#include <iostream>
#include "kyber_test.h"
#include "kyber.h"
#include "sha3.h"


NAMESPACE_BEGIN(CryptoPP)
NAMESPACE_BEGIN(Test)


#define NTESTS 1000

int TestKeys512()
{
  unsigned int i;
  unsigned char pk[Kyber512::PUBLICKEYBYTES];
  unsigned char sk[Kyber512::SECRETKEYBYTES];
  unsigned char ct[Kyber512::CIPHERTEXTBYTES];
  unsigned char key_a[Kyber512::SHAREDSECRETBYTES];
  unsigned char key_b[Kyber512::SHAREDSECRETBYTES];
  Kyber512 kyber = Kyber512();
  bool pass = true;
  
  for(i=0;i<NTESTS;i++) {
    //Alice generates a public key
    kyber.KemKeypair(pk, sk);

    //Bob derives a secret key and creates a response
    kyber.KemEnc(ct, key_b, pk);

    //Alice uses Bobs response to get her shared key
    kyber.KemDec(key_a, ct, sk);

    if(memcmp(key_a, key_b, Kyber512::SHAREDSECRETBYTES)) {

        pass = false;
    }
      
  }
  if (pass) {
    return 0;
  } else {
    std::cout << "Key error in Kyber-512" << std::endl;
    return 1;
  }
}

int TestKeys768()
{
  unsigned int i;
  unsigned char pk[Kyber768::PUBLICKEYBYTES];
  unsigned char sk[Kyber768::SECRETKEYBYTES];
  unsigned char ct[Kyber768::CIPHERTEXTBYTES];
  unsigned char key_a[Kyber768::SHAREDSECRETBYTES];
  unsigned char key_b[Kyber768::SHAREDSECRETBYTES];
  Kyber768 kyber = Kyber768();
  bool pass = true;

  for(i=0;i<NTESTS;i++) {
    //Alice generates a public key
    kyber.KemKeypair(pk, sk);

    //Bob derives a secret key and creates a response
    kyber.KemEnc(ct, key_b, pk);

    //Alice uses Bobs response to get her shared key
    kyber.KemDec(key_a, ct, sk);

    if(memcmp(key_a, key_b, Kyber768::SHAREDSECRETBYTES)) {

        pass = false;
    }
      
  }
  if (pass) {
    return 0;
  } else {
    std::cout << "Key error in Kyber-768" << std::endl;
    return 1;
  }
}

int TestKeys1024()
{
  unsigned int i;
  unsigned char pk[Kyber1024::PUBLICKEYBYTES];
  unsigned char sk[Kyber1024::SECRETKEYBYTES];
  unsigned char ct[Kyber1024::CIPHERTEXTBYTES];
  unsigned char key_a[Kyber1024::SHAREDSECRETBYTES];
  unsigned char key_b[Kyber1024::SHAREDSECRETBYTES];
  Kyber1024 kyber = Kyber1024();
  bool pass = true;

  for(i=0;i<NTESTS;i++) {
    //Alice generates a public key
    kyber.KemKeypair(pk, sk);

    //Bob derives a secret key and creates a response
    kyber.KemEnc(ct, key_b, pk);

    //Alice uses Bobs response to get her shared key
    kyber.KemDec(key_a, ct, sk);

    if(memcmp(key_a, key_b, Kyber1024::SHAREDSECRETBYTES)) {

        pass = false;
    }
      
  }
  if (pass) {
    return 0;
  } else {
    std::cout << "Key error in Kyber-1024" << std::endl;
    return 1;
  }
}

static int TestInvalidSkA()
{
  unsigned int i;
  unsigned char pk[Kyber768::PUBLICKEYBYTES];
  unsigned char sk[Kyber768::SECRETKEYBYTES];
  unsigned char ct[Kyber768::CIPHERTEXTBYTES];
  unsigned char key_a[Kyber768::SHAREDSECRETBYTES];
  unsigned char key_b[Kyber768::SHAREDSECRETBYTES];

  Kyber768 kyber = Kyber768();
  bool pass = true;

  for(i=0;i<NTESTS;i++) {
    //Alice generates a public key
    kyber.KemKeypair(pk, sk);
    //Bob derives a secret key and creates a response
    kyber.KemEnc(ct, key_b, pk);

    //Replace secret key with random values
    
    kyber.randombytes(sk, Kyber768::SECRETKEYBYTES);
    //Alice uses Bobs response to get her shared key
    kyber.KemDec(key_a, ct, sk);

    if(!memcmp(key_a, key_b, Kyber768::SHAREDSECRETBYTES)) {
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
  unsigned char pk[Kyber768::PUBLICKEYBYTES];
  unsigned char sk[Kyber768::SECRETKEYBYTES];
  unsigned char ct[Kyber768::CIPHERTEXTBYTES];
  unsigned char key_a[Kyber768::SHAREDSECRETBYTES];
  unsigned char key_b[Kyber768::SHAREDSECRETBYTES];
  size_t pos;

  Kyber768 kyber = Kyber768();


  bool pass = true;

  for(i=0;i<NTESTS;i++) {
    kyber.randombytes((unsigned char *)&pos, sizeof(size_t));

    //Alice generates a public key
    kyber.KemKeypair(pk, sk);

    //Bob derives a secret key and creates a response
    kyber.KemEnc(ct, key_b, pk);

    //Change some byte in the ciphertext (i.e., encapsulated key)
    ct[pos % Kyber768::CIPHERTEXTBYTES] ^= 23;

    //Alice uses Bobs response to get her shared key
    kyber.KemDec(key_a, ct, sk);

    if(!memcmp(key_a, key_b, Kyber768::SHAREDSECRETBYTES)) {
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




int RunKyberTests()
{
  std::cout << "Testing that Kyber works correctly" << std::endl;
  int fail = TestKeys512();
  int fail2 = TestKeys768();
  int fail3 = TestKeys1024();
  int fail4 = TestInvalidSkA();
  int fail5 = TestInvalidCiphertext();
  int failure = fail + fail2 + fail3 + fail4 + fail5;
  if (failure > 0){
    std::cout << "Kyber tests failed" << std::endl;
    return 1;
  }


  std::cout << "All Kyber tests passed" << std::endl;
  
    

  return 0;
}



NAMESPACE_END
NAMESPACE_END