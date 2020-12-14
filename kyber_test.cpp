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

int TestKeys()
{
  unsigned int i;
  unsigned char pk[KYBER_PUBLICKEYBYTES];
  unsigned char sk[KYBER_SECRETKEYBYTES];
  unsigned char ct[KYBER_CIPHERTEXTBYTES];
  unsigned char key_a[KYBER_SSBYTES];
  unsigned char key_b[KYBER_SSBYTES];
  Kyber768 kyber = Kyber768();
  bool pass = true;

  for(i=0;i<NTESTS;i++) {
    //Alice generates a public key
    kyber.KemKeypair(pk, sk);

    //Bob derives a secret key and creates a response
    kyber.KemEnc(ct, key_b, pk);

    //Alice uses Bobs response to get her shared key
    kyber.KemDec(key_a, ct, sk);

    if(memcmp(key_a, key_b, KYBER_SSBYTES)) {
        pass = false;
    }
      
  }
  if (pass) {
    return 0;
  } else {
    std::cout << "Key error" << std::endl;
    return 1;
  }
}

static int TestInvalidSkA()
{
  unsigned int i;
  unsigned char pk[KYBER_PUBLICKEYBYTES];
  unsigned char sk[KYBER_SECRETKEYBYTES];
  unsigned char ct[KYBER_CIPHERTEXTBYTES];
  unsigned char key_a[KYBER_SSBYTES];
  unsigned char key_b[KYBER_SSBYTES];

  Kyber768 kyber = Kyber768();
  bool pass = true;

  for(i=0;i<NTESTS;i++) {
    //Alice generates a public key
    kyber.KemKeypair(pk, sk);
    //Bob derives a secret key and creates a response
    kyber.KemEnc(ct, key_b, pk);

    //Replace secret key with random values
    
    kyber.randombytes(sk, KYBER_SECRETKEYBYTES);
    //Alice uses Bobs response to get her shared key
    kyber.KemDec(key_a, ct, sk);

    if(!memcmp(key_a, key_b, KYBER_SSBYTES)) {
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
  unsigned char pk[KYBER_PUBLICKEYBYTES];
  unsigned char sk[KYBER_SECRETKEYBYTES];
  unsigned char ct[KYBER_CIPHERTEXTBYTES];
  unsigned char key_a[KYBER_SSBYTES];
  unsigned char key_b[KYBER_SSBYTES];
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
    ct[pos % KYBER_CIPHERTEXTBYTES] ^= 23;

    //Alice uses Bobs response to get her shared key
    kyber.KemDec(key_a, ct, sk);

    if(!memcmp(key_a, key_b, KYBER_SSBYTES)) {
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
  int fail = TestKeys();
  int fail2 = TestInvalidSkA();
  int fail3 = TestInvalidCiphertext();
  int failure = fail + fail2 + fail3;
  if (failure > 0){
    std::cout << "Kyber tests failed" << std::endl;
    return 1;
  }


  std::cout << "All Kyber tests passed" << std::endl;
  
    

  return 0;
}

NAMESPACE_END
NAMESPACE_END