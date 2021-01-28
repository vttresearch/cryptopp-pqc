/** The C++ file for the tests of CRYSTALS-Dilithium. 
Adapted by Julius Hekkala from the public-domain reference
implementation of Dilithium by the CRYSTALS team 
(https://github.com/pq-crystals/dilithium) */

#include <iostream>

#include "dilithium_test.h"
#include "dilithium.h"


NAMESPACE_BEGIN(CryptoPP)
NAMESPACE_BEGIN(Test)

#define MLEN 59
#define NTESTS 1000

int Dilithium1Test()
{
  std::cout << "Running Dilithium1 tests" << std::endl;
  word32 i, j;
  int ret;
  size_t mLen, smLen;
  byte m[MLEN];
  byte sm[MLEN + Dilithium1::BYTES];
  byte m2[MLEN + Dilithium1::BYTES];
  byte pk[Dilithium1::PUBLICKEYBYTES];
  byte sk[Dilithium1::SECRETKEYBYTES];

  Dilithium1 dilithium = Dilithium1();
  for(i = 0; i < NTESTS; ++i) {
    dilithium.RandomBytes(m, MLEN);
    
    dilithium.Keypair(pk, sk);
    
    dilithium.Sign(sm, &smLen, m, MLEN, sk);
    
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);

    if(ret) {
      std::cout << "Dilithium1: Vefication failed" << std::endl;
      return 1;
    }

    if(mLen != MLEN) {
      std::cout << "Dilithium1: Message lengths don't match" << std::endl;
      return 1;
    }

    for(j = 0; j < mLen; ++j) {
      if(m[j] != m2[j]) {
        std::cout << "Dilithium1: Messages don't match" << std::endl;
        return 1;
      }
    }

    dilithium.RandomBytes((byte *)&j, sizeof(j));
    do {
      dilithium.RandomBytes(m2, 1);
    } while(!m2[0]);
    sm[j % Dilithium1::BYTES] += m2[0];
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);
    if(!ret) {
      std::cout << "Dilithium1: Trivial forgeries possible" << std::endl;
      return 1;
    }
  }
  std::cout << "Dilithium1 tests passed" << std::endl;
  return 0;
}

int Dilithium2Test()
{
  std::cout << "Running Dilithium2 tests" << std::endl;
  word32 i, j;
  int ret;
  size_t mLen, smLen;
  byte m[MLEN];
  byte sm[MLEN + Dilithium2::BYTES];
  byte m2[MLEN + Dilithium2::BYTES];
  byte pk[Dilithium2::PUBLICKEYBYTES];
  byte sk[Dilithium2::SECRETKEYBYTES];

  Dilithium2 dilithium = Dilithium2();
  for(i = 0; i < NTESTS; ++i) {
    dilithium.RandomBytes(m, MLEN);
    
    dilithium.Keypair(pk, sk);
    
    dilithium.Sign(sm, &smLen, m, MLEN, sk);
    
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);

    if(ret) {
      std::cout << "Dilithium2: Vefication failed" << std::endl;
      return 1;
    }

    if(mLen != MLEN) {
      std::cout << "Dilithium2: Message lengths don't match" << std::endl;
      return 1;
    }

    for(j = 0; j < mLen; ++j) {
      if(m[j] != m2[j]) {
        std::cout << "Dilithium2: Messages don't match" << std::endl;
        return 1;
      }
    }

    dilithium.RandomBytes((byte *)&j, sizeof(j));
    do {
      dilithium.RandomBytes(m2, 1);
    } while(!m2[0]);
    sm[j % Dilithium2::BYTES] += m2[0];
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);
    if(!ret) {
      std::cout << "Dilithium2: Trivial forgeries possible" << std::endl;
      return 1;
    }
  }
  std::cout << "Dilithium2 tests passed" << std::endl;
  return 0;
}

int Dilithium3Test()
{
  std::cout << "Running Dilithium3 tests" << std::endl;
  word32 i, j;
  int ret;
  size_t mLen, smLen;
  byte m[MLEN];
  byte sm[MLEN + Dilithium3::BYTES];
  byte m2[MLEN + Dilithium3::BYTES];
  byte pk[Dilithium3::PUBLICKEYBYTES];
  byte sk[Dilithium3::SECRETKEYBYTES];

  Dilithium3 dilithium = Dilithium3();
  for(i = 0; i < NTESTS; ++i) {
    dilithium.RandomBytes(m, MLEN);
    
    dilithium.Keypair(pk, sk);
    
    dilithium.Sign(sm, &smLen, m, MLEN, sk);
    
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);

    if(ret) {
      std::cout << "Dilithium3: Vefication failed" << std::endl;
      return 1;
    }

    if(mLen != MLEN) {
      std::cout << "Dilithium3: Message lengths don't match" << std::endl;
      return 1;
    }

    for(j = 0; j < mLen; ++j) {
      if(m[j] != m2[j]) {
        std::cout << "Dilithium3: Messages don't match" << std::endl;
        return 1;
      }
    }

    dilithium.RandomBytes((byte *)&j, sizeof(j));
    do {
      dilithium.RandomBytes(m2, 1);
    } while(!m2[0]);
    sm[j % Dilithium3::BYTES] += m2[0];
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);
    if(!ret) {
      std::cout << "Dilithium3: Trivial forgeries possible" << std::endl;
      return 1;
    }
  }
  std::cout << "Dilithium3 tests passed" << std::endl;
  return 0;
}

int Dilithium4Test()
{
  std::cout << "Running Dilithium4 tests" << std::endl;
  word32 i, j;
  int ret;
  size_t mLen, smLen;
  byte m[MLEN];
  byte sm[MLEN + Dilithium4::BYTES];
  byte m2[MLEN + Dilithium4::BYTES];
  byte pk[Dilithium4::PUBLICKEYBYTES];
  byte sk[Dilithium4::SECRETKEYBYTES];

  Dilithium4 dilithium = Dilithium4();
  for(i = 0; i < NTESTS; ++i) {
    dilithium.RandomBytes(m, MLEN);
    
    dilithium.Keypair(pk, sk);
    
    dilithium.Sign(sm, &smLen, m, MLEN, sk);
    
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);

    if(ret) {
      std::cout << "Dilithium4: Vefication failed" << std::endl;
      return 1;
    }

    if(mLen != MLEN) {
      std::cout << "Dilithium4: Message lengths don't match" << std::endl;
      return 1;
    }

    for(j = 0; j < mLen; ++j) {
      if(m[j] != m2[j]) {
        std::cout << "Dilithium4: Messages don't match" << std::endl;
        return 1;
      }
    }

    dilithium.RandomBytes((byte *)&j, sizeof(j));
    do {
      dilithium.RandomBytes(m2, 1);
    } while(!m2[0]);
    sm[j % Dilithium4::BYTES] += m2[0];
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);
    if(!ret) {
      std::cout << "Dilithium4: Trivial forgeries possible" << std::endl;
      return 1;
    }
  }
  std::cout << "Dilithium4 tests passed" << std::endl;
  return 0;
}




//Test different modes of Dilithium
int RunDilithiumTests() {
    std::cout << "Testing that Dilithium works correctly" << std::endl;
    int dilithium1 = Dilithium1Test();
    int dilithium2 = Dilithium2Test();
    int dilithium3 = Dilithium3Test();
    int dilithium4 = Dilithium4Test();
    int ret = dilithium1 + dilithium2 + dilithium3 + dilithium4;
    if (ret != 0) {  
        std::cout << "Dilithium tests failed." << std::endl;
        return 1;
    }
    std::cout << "All Dilithium tests passed" << std::endl;
    return 0;
}

NAMESPACE_END
NAMESPACE_END