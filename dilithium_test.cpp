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



int Dilithium2Test()
{
  std::cout << "Running Dilithium2 tests" << std::endl;
  word32 i, j;
  int ret;
  size_t mLen, smLen;
  byte m[MLEN + Dilithium2::BYTES];
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
    sm[j % (MLEN + Dilithium2::BYTES)] += m2[0];
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
  byte m[MLEN + Dilithium3::BYTES];
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
    sm[j % (MLEN + Dilithium3::BYTES)] += m2[0];
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);
    if(!ret) {
      std::cout << "Dilithium3: Trivial forgeries possible" << std::endl;
      return 1;
    }
  }
  std::cout << "Dilithium3 tests passed" << std::endl;
  return 0;
}

int Dilithium5Test()
{
  std::cout << "Running Dilithium5 tests" << std::endl;
  word32 i, j;
  int ret;
  size_t mLen, smLen;
  byte m[MLEN + Dilithium5::BYTES];
  byte sm[MLEN + Dilithium5::BYTES];
  byte m2[MLEN + Dilithium5::BYTES];
  byte pk[Dilithium5::PUBLICKEYBYTES];
  byte sk[Dilithium5::SECRETKEYBYTES];

  Dilithium5 dilithium = Dilithium5();
  for(i = 0; i < NTESTS; ++i) {
    dilithium.RandomBytes(m, MLEN);
    
    dilithium.Keypair(pk, sk);
    
    dilithium.Sign(sm, &smLen, m, MLEN, sk);
    
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);

    if(ret) {
      std::cout << "Dilithium5: Vefication failed" << std::endl;
      return 1;
    }

    if(mLen != MLEN) {
      std::cout << "Dilithium5: Message lengths don't match" << std::endl;
      return 1;
    }

    for(j = 0; j < mLen; ++j) {
      if(m[j] != m2[j]) {
        std::cout << "Dilithium5: Messages don't match" << std::endl;
        return 1;
      }
    }

    dilithium.RandomBytes((byte *)&j, sizeof(j));
    do {
      dilithium.RandomBytes(m2, 1);
    } while(!m2[0]);
    sm[j % (MLEN + Dilithium5::BYTES)] += m2[0];
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);
    if(!ret) {
      std::cout << "Dilithium5: Trivial forgeries possible" << std::endl;
      return 1;
    }
  }
  std::cout << "Dilithium5 tests passed" << std::endl;
  return 0;
}

/**
 * Test the randomized signature mode of Dilithium2
 */
int Dilithium2TestRandomized() {
  std::cout << "Running Dilithium2 randomized tests" << std::endl;
  word32 i, j;
  int ret;
  size_t mLen, smLen;
  byte m[MLEN + Dilithium2::BYTES];
  byte sm[MLEN + Dilithium2::BYTES];
  byte m2[MLEN + Dilithium2::BYTES];
  byte pk[Dilithium2::PUBLICKEYBYTES];
  byte sk[Dilithium2::SECRETKEYBYTES];

  Dilithium2 dilithium = Dilithium2();
  dilithium.setRandomizedSignature();
  dilithium.RandomBytes(m, MLEN);
  for(i = 0; i < NTESTS; ++i) {
    dilithium.Keypair(pk, sk);
    
    dilithium.Sign(sm, &smLen, m, MLEN, sk);
    
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);

    if(ret) {
      std::cout << "Randomized Dilithium2: Vefication failed" << std::endl;
      return 1;
    }

    if(mLen != MLEN) {
      std::cout << "Randomized Dilithium2: Message lengths don't match" << std::endl;
      return 1;

    for(j = 0; j < mLen; ++j) {
      if(m[j] != m2[j]) {
    }
        std::cout << "Randomized Dilithium2: Messages don't match" << std::endl;
        return 1;
      }
    }

    dilithium.RandomBytes((byte *)&j, sizeof(j));
    do {
      dilithium.RandomBytes(m2, 1);
    } while(!m2[0]);
    sm[j % (MLEN + Dilithium2::BYTES)] += m2[0];
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);
    if(!ret) {
      std::cout << "Randomized Dilithium2: Trivial forgeries possible" << std::endl;
      return 1;
    }
  }
  std::cout << "Randomized Dilithium2 tests passed" << std::endl;
  return 0;
}


/**
 * Test the randomized signature mode of Dilithium3
 */
int Dilithium3TestRandomized() {
  std::cout << "Running Dilithium3 randomized tests" << std::endl;
  word32 i, j;
  int ret;
  size_t mLen, smLen;
  byte m[MLEN + Dilithium3::BYTES];
  byte sm[MLEN + Dilithium3::BYTES];
  byte m2[MLEN + Dilithium3::BYTES];
  byte pk[Dilithium3::PUBLICKEYBYTES];
  byte sk[Dilithium3::SECRETKEYBYTES];

  Dilithium3 dilithium = Dilithium3();
  dilithium.setRandomizedSignature();
  dilithium.RandomBytes(m, MLEN);
  for(i = 0; i < NTESTS; ++i) {
    dilithium.Keypair(pk, sk);
    
    dilithium.Sign(sm, &smLen, m, MLEN, sk);
    
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);

    if(ret) {
      std::cout << "Randomized Dilithium3: Vefication failed" << std::endl;
      return 1;
    }

    if(mLen != MLEN) {
      std::cout << "Randomized Dilithium3: Message lengths don't match" << std::endl;
      return 1;
    }

    for(j = 0; j < mLen; ++j) {
      if(m[j] != m2[j]) {
        std::cout << "Randomized Dilithium3: Messages don't match" << std::endl;
        return 1;
      }
    }

    dilithium.RandomBytes((byte *)&j, sizeof(j));
    do {
      dilithium.RandomBytes(m2, 1);
    } while(!m2[0]);
    sm[j % (MLEN + Dilithium3::BYTES)] += m2[0];
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);
    if(!ret) {
      std::cout << "Randomized Dilithium3: Trivial forgeries possible" << std::endl;
      return 1;
    }
  }
  std::cout << "Randomized Dilithium3 tests passed" << std::endl;
  return 0;
}

/**
 * Test the randomized signature mode of Dilithium5
 */
int Dilithium5TestRandomized() {
  std::cout << "Running Dilithium5 randomized tests" << std::endl;
  word32 i, j;
  int ret;
  size_t mLen, smLen;
  byte m[MLEN + Dilithium5::BYTES];
  byte sm[MLEN + Dilithium5::BYTES];
  byte m2[MLEN + Dilithium5::BYTES];
  byte pk[Dilithium5::PUBLICKEYBYTES];
  byte sk[Dilithium5::SECRETKEYBYTES];

  Dilithium5 dilithium = Dilithium5();
  dilithium.setRandomizedSignature();
  dilithium.RandomBytes(m, MLEN);
  for(i = 0; i < NTESTS; ++i) {
    dilithium.Keypair(pk, sk);
    
    dilithium.Sign(sm, &smLen, m, MLEN, sk);
    
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);

    if(ret) {
      std::cout << "Randomized Dilithium5: Vefication failed" << std::endl;
      return 1;
    }

    if(mLen != MLEN) {
      std::cout << "Randomized Dilithium5: Message lengths don't match" << std::endl;
      return 1;
    }

    for(j = 0; j < mLen; ++j) {
      if(m[j] != m2[j]) {
        std::cout << "Randomized Dilithium5: Messages don't match" << std::endl;
        return 1;
      }
    }

    dilithium.RandomBytes((byte *)&j, sizeof(j));
    do {
      dilithium.RandomBytes(m2, 1);
    } while(!m2[0]);
    sm[j % (MLEN + Dilithium5::BYTES)] += m2[0];
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);
    if(!ret) {
      std::cout << "Randomized Dilithium5: Trivial forgeries possible" << std::endl;
      return 1;
    }
  }
  std::cout << "Randomized Dilithium5 tests passed" << std::endl;
  return 0;
}

//Test different modes of Dilithium
int RunDilithiumTests() {
    std::cout << "Testing that Dilithium works correctly" << std::endl;
    int dilithium2 = Dilithium2Test();
    int dilithium3 = Dilithium3Test();
    int dilithium5 = Dilithium5Test();
    int dilithium2random = Dilithium2TestRandomized();
    int dilithium3random = Dilithium3TestRandomized();
    int dilithium5random = Dilithium5TestRandomized();

    int ret = dilithium2 + dilithium3 + dilithium5 + dilithium2random + dilithium3random + dilithium5random;
    if (ret != 0) {  
        std::cout << "Dilithium tests failed." << std::endl;
        return 1;
    }
    std::cout << "All Dilithium tests passed" << std::endl;
    return 0;
}

NAMESPACE_END
NAMESPACE_END