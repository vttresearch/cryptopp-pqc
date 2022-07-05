/*
 C++ file for benchmarking PQC implementations
 Implemented by Julius Hekkala. The benchmarks are based on the tests
 of the algorithms. The tests are based on the reference implementations of 
 kyber by the CRYSTALS team (https://github.com/pq-crystals/kyber),
 Dilithium by the CRYSTALS team (https://github.com/pq-crystals/dilithium)
 and SABER by the SABER team (https://github.com/KULeuven-COSIC/SABER)
*/

#include "pch.h"
#include "pqc_benchmark.h"
#include "dilithium.h"
#include "kyber.h"
#include "saber.h"
#include <iostream>
#include <fstream>

#ifdef _WIN32
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

NAMESPACE_BEGIN(CryptoPP)
NAMESPACE_BEGIN(Test)

#define NTESTS 10000


/**
 * Dilithium benchmarks
 */
void BenchmarkDilithiumVersions() {
  std::cout << "Benchmarking Dilithium.\n";
  std::cout << "Results will be saved in benchmark_dilithium.txt" << std::endl;
  BenchmarkDilithium2();
  BenchmarkDilithium3();
  BenchmarkDilithium5();
}

/**
 * Kyber benchmarks
 */
void BenchmarkKyberVersions() {
  std::cout << "Benchmarking Kyber.\n";
  std::cout << "Results will be saved in benchmark_kyber.txt" << std::endl;
  BenchmarkKyber512();
  BenchmarkKyber768();
  BenchmarkKyber1024();

}

/**
 * Saber benchmarks
 */
void BenchmarkSaberVersions() {
  std::cout << "Benchmarking SABER.\n";
  std::cout << "Results will be saved in benchmark_saber.txt" << std::endl;
  BenchmarkLightSaber();
  BenchmarkSaber();
  BenchmarkFireSaber();
}


void BenchmarkLightSaber() {
    std::ofstream saberStream("benchmark_saber.txt");
    if (saberStream.is_open()) {
      saberStream << "Benchmarking LightSaber...\n";
      saberStream.close();
    }
    
    word64 clockKeyPair, clockKemEnc, clockKemDec; 
    clockKeyPair = 0;
    clockKemEnc = 0;
    clockKemDec = 0;

    LightSaber saber = LightSaber();

    for (int i=0; i < NTESTS; i++) {
        byte pk[LightSaber::PUBLICKEYBYTES];
        byte sk[LightSaber::SECRETKEYBYTES];
        byte ct[LightSaber::CIPHERTEXTBYTES];	
        byte ssA[LightSaber::SHAREDSECRETBYTES], ssB[LightSaber::SHAREDSECRETBYTES];

        word64 clockBefore;
        word64 clockAfter;
        

        //Generate public key pk and secret key sk
        clockBefore = __rdtsc();
        saber.KemKeypair(pk, sk);
        clockAfter = __rdtsc();
        clockKeyPair = clockKeyPair + (clockAfter - clockBefore);

        //Key encapsulation
        clockBefore = __rdtsc();
        saber.KemEnc(ct, ssA, pk);
        clockAfter = __rdtsc();
        clockKemEnc = clockKemEnc + (clockAfter - clockBefore);


        
        //Key decapsulation
        clockBefore = __rdtsc();
        saber.KemDec(ssB, ct, sk);
        clockAfter = __rdtsc();
        clockKemDec = clockKemDec + (clockAfter - clockBefore);
        if (memcmp(ssA, ssB, LightSaber::SHAREDSECRETBYTES)) {
            std::cout << "Error in LightSaber" << std::endl;
        }

    }
    saberStream.open("benchmark_saber.txt", std::ios::app);
    if (saberStream.is_open()) {
      saberStream << "Number of tests: ";
      saberStream << NTESTS << "\n";
      saberStream << "Average keypair time: ";
      saberStream << clockKeyPair / NTESTS << "\n";
      saberStream << "Average KemEnc time: ";
      saberStream << clockKemEnc / NTESTS << "\n";
      saberStream << "Average KemDec time: ";
      saberStream << clockKemDec / NTESTS << "\n";
      saberStream << "\n";
      saberStream.close();
    

    }
    

}

void BenchmarkSaber() {
    std::ofstream saberStream("benchmark_saber.txt", std::ios::app);
    if (saberStream.is_open()) {
      saberStream << "Benchmarking Saber...\n";
      saberStream.close();
    }
    
    word64 clockKeyPair, clockKemEnc, clockKemDec; 
    Saber saber = Saber();
    clockKeyPair = 0;
    clockKemEnc = 0;
    clockKemDec = 0;

    for (int i=0; i < NTESTS; i++) {
        byte pk[Saber::PUBLICKEYBYTES];
        byte sk[Saber::SECRETKEYBYTES];
        byte ct[Saber::CIPHERTEXTBYTES];	
        byte ssA[Saber::SHAREDSECRETBYTES], ssB[Saber::SHAREDSECRETBYTES];

        word64 clockBefore;
        word64 clockAfter;
        

        //Generate public key pk and secret key sk
        clockBefore = __rdtsc();
        saber.KemKeypair(pk, sk);
        clockAfter = __rdtsc();
        clockKeyPair = clockKeyPair + (clockAfter - clockBefore);

        //Key encapsulation
        clockBefore = __rdtsc();
        saber.KemEnc(ct, ssA, pk);
        clockAfter = __rdtsc();
        clockKemEnc = clockKemEnc + (clockAfter - clockBefore);


        
        //Key decapsulation
        clockBefore = __rdtsc();
        saber.KemDec(ssB, ct, sk);
        clockAfter = __rdtsc();
        clockKemDec = clockKemDec + (clockAfter - clockBefore);
        if (memcmp(ssA, ssB, Saber::SHAREDSECRETBYTES)) {
            std::cout << "Error in Saber" << std::endl;
        }

    }
    saberStream.open("benchmark_saber.txt", std::ios::app);
    if (saberStream.is_open()) {
      saberStream << "Number of tests: ";
      saberStream << NTESTS << "\n";
      saberStream << "Average keypair time: ";
      saberStream << clockKeyPair / NTESTS << "\n";
      saberStream << "Average KemEnc time: ";
      saberStream << clockKemEnc / NTESTS << "\n";
      saberStream << "Average KemDec time: ";
      saberStream << clockKemDec / NTESTS << "\n";
      saberStream << "\n";
      saberStream.close();
    }
    
    

}

void BenchmarkFireSaber() {
    std::ofstream saberStream("benchmark_saber.txt", std::ios::app);
    if (saberStream.is_open()) {
      saberStream << "Benchmarking FireSaber...\n";
      saberStream.close();
    }
    
    word64 clockKeyPair, clockKemEnc, clockKemDec; 
    clockKeyPair = 0;
    clockKemEnc = 0;
    clockKemDec = 0;

    FireSaber saber = FireSaber();

    for (int i=0; i < NTESTS; i++) {
        byte pk[FireSaber::PUBLICKEYBYTES];
        byte sk[FireSaber::SECRETKEYBYTES];
        byte ct[FireSaber::CIPHERTEXTBYTES];	
        byte ssA[FireSaber::SHAREDSECRETBYTES], ssB[FireSaber::SHAREDSECRETBYTES];

        word64 clockBefore;
        word64 clockAfter;
        

        //Generate public key pk and secret key sk
        clockBefore = __rdtsc();
        saber.KemKeypair(pk, sk);
        clockAfter = __rdtsc();
        clockKeyPair = clockKeyPair + (clockAfter - clockBefore);

        //Key encapsulation
        clockBefore = __rdtsc();
        saber.KemEnc(ct, ssA, pk);
        clockAfter = __rdtsc();
        clockKemEnc = clockKemEnc + (clockAfter - clockBefore);


        
        //Key decapsulation
        clockBefore = __rdtsc();
        saber.KemDec(ssB, ct, sk);
        clockAfter = __rdtsc();
        clockKemDec = clockKemDec + (clockAfter - clockBefore);
        if (memcmp(ssA, ssB, FireSaber::SHAREDSECRETBYTES)) {
            std::cout << "Error in FireSaber" << std::endl;
        }

    }
    saberStream.open("benchmark_saber.txt", std::ios::app);
    if (saberStream.is_open()) {
      saberStream << "Number of tests: ";
      saberStream << NTESTS << "\n";
      saberStream << "Average keypair time: ";
      saberStream << clockKeyPair / NTESTS << "\n";
      saberStream << "Average KemEnc time: ";
      saberStream << clockKemEnc / NTESTS << "\n";
      saberStream << "Average KemDec time: ";
      saberStream << clockKemDec / NTESTS << "\n";
      saberStream << "\n";  
      saberStream.close();
    }
      

}

#define MLEN 59


void BenchmarkDilithium2()
{
  std::ofstream dilithiumStream("benchmark_dilithium.txt");
  if (dilithiumStream.is_open()) {
    dilithiumStream << "Benchmarking Dilithium2\n";
    dilithiumStream.close();
  }
  word32 i, j;
  int ret;
  size_t mLen, smLen;
  byte m[MLEN];
  byte sm[MLEN + Dilithium2::BYTES];
  byte m2[MLEN + Dilithium2::BYTES];
  byte pk[Dilithium2::PUBLICKEYBYTES];
  byte sk[Dilithium2::SECRETKEYBYTES];
  word64 clockKeypair, clockSign, clockOpen;
  clockKeypair = 0;
  clockSign = 0;
  clockOpen = 0;

  Dilithium2 dilithium = Dilithium2();
  for(i = 0; i < NTESTS; ++i) {
    word64 clockBefore, clockAfter;

    dilithium.RandomBytes(m, MLEN);
    
    clockBefore = __rdtsc();
    dilithium.Keypair(pk, sk);
    clockAfter = __rdtsc();
    clockKeypair = clockKeypair + (clockAfter - clockBefore);

    clockBefore = __rdtsc();
    dilithium.Sign(sm, &smLen, m, MLEN, sk);
    clockAfter = __rdtsc();
    clockSign = clockSign + (clockAfter - clockBefore);
    
    clockBefore = __rdtsc();
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);
    clockAfter = __rdtsc();
    clockOpen = clockOpen + (clockAfter - clockBefore);

    if(ret) {
      std::cout << "Dilithium2: Vefication failed" << std::endl;
    }

    if(mLen != MLEN) {
      std::cout << "Dilithium2: Message lengths don't match" << std::endl;
    }

    for(j = 0; j < mLen; ++j) {
      if(m[j] != m2[j]) {
        std::cout << "Dilithium2: Messages don't match" << std::endl;
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
    }
  }
  dilithiumStream.open("benchmark_dilithium.txt", std::ios::app);
  if (dilithiumStream.is_open()) {
    dilithiumStream << "Number of tests: ";
    dilithiumStream << NTESTS << "\n";
    dilithiumStream << "Average keypair time: ";
    dilithiumStream << clockKeypair / NTESTS << "\n";
    dilithiumStream << "Average sign time: ";
    dilithiumStream << clockSign / NTESTS << "\n";
    dilithiumStream << "Average open time: ";
    dilithiumStream << clockOpen / NTESTS << "\n";
    dilithiumStream << "\n";
    dilithiumStream.close();
  }
  
  
}


void BenchmarkDilithium3()
{
  std::ofstream dilithiumStream("benchmark_dilithium.txt", std::ios::app);
  if (dilithiumStream.is_open()) {
    dilithiumStream << "Benchmarking Dilithium3\n";
    dilithiumStream.close();
  }
  word32 i, j;
  int ret;
  size_t mLen, smLen;
  byte m[MLEN];
  byte sm[MLEN + Dilithium3::BYTES];
  byte m2[MLEN + Dilithium3::BYTES];
  byte pk[Dilithium3::PUBLICKEYBYTES];
  byte sk[Dilithium3::SECRETKEYBYTES];
  word64 clockKeypair, clockSign, clockOpen;
  clockKeypair = 0;
  clockSign = 0;
  clockOpen = 0;

  Dilithium3 dilithium = Dilithium3();
  for(i = 0; i < NTESTS; ++i) {
    word64 clockBefore, clockAfter;

    dilithium.RandomBytes(m, MLEN);
    
    clockBefore = __rdtsc();
    dilithium.Keypair(pk, sk);
    clockAfter = __rdtsc();
    clockKeypair = clockKeypair + (clockAfter - clockBefore);

    clockBefore = __rdtsc();
    dilithium.Sign(sm, &smLen, m, MLEN, sk);
    clockAfter = __rdtsc();
    clockSign = clockSign + (clockAfter - clockBefore);
    
    clockBefore = __rdtsc();
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);
    clockAfter = __rdtsc();
    clockOpen = clockOpen + (clockAfter - clockBefore);

    if(ret) {
      std::cout << "Dilithium3: Vefication failed" << std::endl;
    }

    if(mLen != MLEN) {
      std::cout << "Dilithium3: Message lengths don't match" << std::endl;
    }

    for(j = 0; j < mLen; ++j) {
      if(m[j] != m2[j]) {
        std::cout << "Dilithium3: Messages don't match" << std::endl;
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
    }
  }
  dilithiumStream.open("benchmark_dilithium.txt", std::ios::app);
  if (dilithiumStream.is_open()) {
    dilithiumStream << "Number of tests: ";
    dilithiumStream << NTESTS << "\n";
    dilithiumStream << "Average keypair time: ";
    dilithiumStream << clockKeypair / NTESTS << "\n";
    dilithiumStream << "Average sign time: ";
    dilithiumStream << clockSign / NTESTS << "\n";
    dilithiumStream << "Average open time: ";
    dilithiumStream << clockOpen / NTESTS << "\n";
    dilithiumStream << "\n";
    dilithiumStream.close();
  }
  
    
  
}

void BenchmarkDilithium5()
{
  std::ofstream dilithiumStream("benchmark_dilithium.txt", std::ios::app);
  if (dilithiumStream.is_open()) {
    dilithiumStream << "Benchmarking Dilithium5\n";
    dilithiumStream.close();
  }
  
  word32 i, j;
  int ret;
  size_t mLen, smLen;
  byte m[MLEN];
  byte sm[MLEN + Dilithium5::BYTES];
  byte m2[MLEN + Dilithium5::BYTES];
  byte pk[Dilithium5::PUBLICKEYBYTES];
  byte sk[Dilithium5::SECRETKEYBYTES];
  word64 clockKeypair, clockSign, clockOpen;
  clockKeypair = 0;
  clockSign = 0;
  clockOpen = 0;

  Dilithium5 dilithium = Dilithium5();
  for(i = 0; i < NTESTS; ++i) {
    word64 clockBefore, clockAfter;

    dilithium.RandomBytes(m, MLEN);
    
    clockBefore = __rdtsc();
    dilithium.Keypair(pk, sk);
    clockAfter = __rdtsc();
    clockKeypair = clockKeypair + (clockAfter - clockBefore);

    clockBefore = __rdtsc();
    dilithium.Sign(sm, &smLen, m, MLEN, sk);
    clockAfter = __rdtsc();
    clockSign = clockSign + (clockAfter - clockBefore);
    
    clockBefore = __rdtsc();
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);
    clockAfter = __rdtsc();
    clockOpen = clockOpen + (clockAfter - clockBefore);

    if(ret) {
      std::cout << "Dilithium5: Vefication failed" << std::endl;
    }

    if(mLen != MLEN) {
      std::cout << "Dilithium5: Message lengths don't match" << std::endl;
    }

    for(j = 0; j < mLen; ++j) {
      if(m[j] != m2[j]) {
        std::cout << "Dilithium5: Messages don't match" << std::endl;
      }
    }

    dilithium.RandomBytes((byte *)&j, sizeof(j));
    do {
      dilithium.RandomBytes(m2, 1);
    } while(!m2[0]);
    sm[j % Dilithium5::BYTES] += m2[0];
    ret = dilithium.Open(m2, &mLen, sm, smLen, pk);
    if(!ret) {
      std::cout << "Dilithium5: Trivial forgeries possible" << std::endl;
    }
  }
  dilithiumStream.open("benchmark_dilithium.txt", std::ios::app);
  if (dilithiumStream.is_open()) {
    dilithiumStream << "Number of tests: ";
    dilithiumStream << NTESTS << "\n";
    dilithiumStream << "Average keypair time: ";
    dilithiumStream << clockKeypair / NTESTS << "\n";
    dilithiumStream << "Average sign time: ";
    dilithiumStream << clockSign / NTESTS << "\n";
    dilithiumStream << "Average open time: ";
    dilithiumStream << clockOpen / NTESTS << "\n";
    dilithiumStream << "\n";
    dilithiumStream.close();
  }

  
}

void BenchmarkKyber512()
{
  std::ofstream kyberStream("benchmark_kyber.txt");
  if (kyberStream.is_open()) {
    kyberStream << "Benchmarking Kyber512...\n";
    kyberStream.close();
  }
  
  unsigned int i;
  unsigned char pk[Kyber512::PUBLICKEYBYTES];
  unsigned char sk[Kyber512::SECRETKEYBYTES];
  unsigned char ct[Kyber512::CIPHERTEXTBYTES];
  unsigned char key_a[Kyber512::SHAREDSECRETBYTES];
  unsigned char key_b[Kyber512::SHAREDSECRETBYTES];
  Kyber512 kyber = Kyber512();

  word64 clockKeypair, clockKemEnc, clockKemDec;
  clockKeypair = 0;
  clockKemEnc = 0;
  clockKemDec = 0;
  
  for(i=0;i<NTESTS;i++) {
    word64 clockBefore, clockAfter;
    //Alice generates a public key
    clockBefore = __rdtsc();
    kyber.KemKeypair(pk, sk);
    clockAfter = __rdtsc();
    clockKeypair = clockKeypair + (clockAfter - clockBefore);

    //Bob derives a secret key and creates a response
    clockBefore = __rdtsc();
    kyber.KemEnc(ct, key_b, pk);
    clockAfter = __rdtsc();
    clockKemEnc = clockKemEnc + (clockAfter - clockBefore);

    //Alice uses Bobs response to get her shared key
    clockBefore = __rdtsc();
    kyber.KemDec(key_a, ct, sk);
    clockAfter = __rdtsc();
    clockKemDec = clockKemDec + (clockAfter - clockBefore);

    if(memcmp(key_a, key_b, Kyber512::SHAREDSECRETBYTES)) {
      std::cout << "Kyber512: Key error" << std::endl;
    }
      
  }

  kyberStream.open("benchmark_kyber.txt", std::ios::app);
  if (kyberStream.is_open()) {
    kyberStream << "Number of tests: ";
    kyberStream << NTESTS << "\n";
    kyberStream << "Average keypair time: ";
    kyberStream << clockKeypair / NTESTS << "\n";
    kyberStream << "Average KemEnc time: ";
    kyberStream << clockKemEnc / NTESTS << "\n";
    kyberStream << "Average KemDec time: ";
    kyberStream << clockKemDec / NTESTS << "\n";
    kyberStream << "\n";
    kyberStream.close();
  }
  
  
}

void BenchmarkKyber768()
{
  std::ofstream kyberStream("benchmark_kyber.txt", std::ios::app);
  if (kyberStream.is_open()) {
    kyberStream << "Benchmarking Kyber768...\n";
    kyberStream.close();
  }
  unsigned int i;
  unsigned char pk[Kyber768::PUBLICKEYBYTES];
  unsigned char sk[Kyber768::SECRETKEYBYTES];
  unsigned char ct[Kyber768::CIPHERTEXTBYTES];
  unsigned char key_a[Kyber768::SHAREDSECRETBYTES];
  unsigned char key_b[Kyber768::SHAREDSECRETBYTES];
  Kyber768 kyber = Kyber768();

  word64 clockKeypair, clockKemEnc, clockKemDec;
  clockKeypair = 0;
  clockKemEnc = 0;
  clockKemDec = 0;
  
  for(i=0;i<NTESTS;i++) {
    word64 clockBefore, clockAfter;
    //Alice generates a public key
    clockBefore = __rdtsc();
    kyber.KemKeypair(pk, sk);
    clockAfter = __rdtsc();
    clockKeypair = clockKeypair + (clockAfter - clockBefore);

    //Bob derives a secret key and creates a response
    clockBefore = __rdtsc();
    kyber.KemEnc(ct, key_b, pk);
    clockAfter = __rdtsc();
    clockKemEnc = clockKemEnc + (clockAfter - clockBefore);

    //Alice uses Bobs response to get her shared key
    clockBefore = __rdtsc();
    kyber.KemDec(key_a, ct, sk);
    clockAfter = __rdtsc();
    clockKemDec = clockKemDec + (clockAfter - clockBefore);

    if(memcmp(key_a, key_b, Kyber768::SHAREDSECRETBYTES)) {
      std::cout << "Kyber768: Key error" << std::endl;
    }
      
  }

  kyberStream.open("benchmark_kyber.txt", std::ios::app);
  if (kyberStream.is_open()) {
    kyberStream << "Number of tests: ";
    kyberStream << NTESTS << "\n";
    kyberStream << "Average keypair time: ";
    kyberStream << clockKeypair / NTESTS << "\n";
    kyberStream << "Average KemEnc time: ";
    kyberStream << clockKemEnc / NTESTS << "\n";
    kyberStream << "Average KemDec time: ";
    kyberStream << clockKemDec / NTESTS << "\n";
    kyberStream << "\n";
    kyberStream.close();
  }
  
  
}

void BenchmarkKyber1024()
{
  std::ofstream kyberStream("benchmark_kyber.txt", std::ios::app);
  if (kyberStream.is_open()) {
    kyberStream << "Benchmarking Kyber1024...\n";
    kyberStream.close();
  }
  
  unsigned int i;
  unsigned char pk[Kyber1024::PUBLICKEYBYTES];
  unsigned char sk[Kyber1024::SECRETKEYBYTES];
  unsigned char ct[Kyber1024::CIPHERTEXTBYTES];
  unsigned char key_a[Kyber1024::SHAREDSECRETBYTES];
  unsigned char key_b[Kyber1024::SHAREDSECRETBYTES];
  Kyber1024 kyber = Kyber1024();

  word64 clockKeypair, clockKemEnc, clockKemDec;
  clockKeypair = 0;
  clockKemEnc = 0;
  clockKemDec = 0;
  
  for(i=0;i<NTESTS;i++) {
    word64 clockBefore, clockAfter;
    //Alice generates a public key
    clockBefore = __rdtsc();
    kyber.KemKeypair(pk, sk);
    clockAfter = __rdtsc();
    clockKeypair = clockKeypair + (clockAfter - clockBefore);

    //Bob derives a secret key and creates a response
    clockBefore = __rdtsc();
    kyber.KemEnc(ct, key_b, pk);
    clockAfter = __rdtsc();
    clockKemEnc = clockKemEnc + (clockAfter - clockBefore);

    //Alice uses Bobs response to get her shared key
    clockBefore = __rdtsc();
    kyber.KemDec(key_a, ct, sk);
    clockAfter = __rdtsc();
    clockKemDec = clockKemDec + (clockAfter - clockBefore);

    if(memcmp(key_a, key_b, Kyber1024::SHAREDSECRETBYTES)) {
      std::cout << "Kyber1024: Key error" << std::endl;
    }
      
  }
  kyberStream.open("benchmark_kyber.txt", std::ios::app);
  if (kyberStream.is_open()) {
    kyberStream << "Number of tests: ";
    kyberStream << NTESTS << "\n";
    kyberStream << "Average keypair time: ";
    kyberStream << clockKeypair / NTESTS << "\n";
    kyberStream << "Average KemEnc time: ";
    kyberStream << clockKemEnc / NTESTS << "\n";
    kyberStream << "Average KemDec time: ";
    kyberStream << clockKemDec / NTESTS << "\n";
    kyberStream << "\n";
    kyberStream.close();
  }
}









NAMESPACE_END
NAMESPACE_END