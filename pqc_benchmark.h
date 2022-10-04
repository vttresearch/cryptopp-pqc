/*
 Header file for benchmarking PQC implementations
 Implemented by Julius Hekkala
*/

#ifndef CRYPTOPP_PQC_BENCHMARK_H
#define CRYPTOPP_PQC_BENCHMARK_H

#include "cryptlib.h"

NAMESPACE_BEGIN(CryptoPP)
NAMESPACE_BEGIN(Test)


void BenchmarkDilithiumVersions();

void BenchmarkKyberVersions();

void BenchmarkSaberVersions();

void BenchmarkFrodoVersions();

void BenchmarkLightSaber();
void BenchmarkSaber();
void BenchmarkFireSaber();


void BenchmarkDilithium2();
void BenchmarkDilithium3();
void BenchmarkDilithium5();

void BenchmarkKyber512();
void BenchmarkKyber768();
void BenchmarkKyber1024();

void BenchmarkFrodo640();
void BenchmarkFrodo976();
void BenchmarkFrodo1344();




NAMESPACE_END
NAMESPACE_END

#endif