/* The C++ file for the post quantum KEM algorithm SABER. Adapted by Julius Hekkala 
from the public-domain reference
implementation of SABER by the SABER team 
(https://github.com/KULeuven-COSIC/SABER) 
 */

#include "pch.h"
#include "saber.h"
#include "sha3.h"
#include "shake.h"
#include "misc.h"

#include <cstring>

NAMESPACE_BEGIN(CryptoPP)



/**
 * Generate a keypair consisting of public key pk and secret key sk.
 * pk and sk are byte arrays provided by the user. 
 * Length of pk: PUBLICKEYBYTES
 * Length of sk: SECRETKEYBYTES
 */ 
template<byte T_L, byte T_Et>
int SABER<T_L, T_Et>::KemKeypair(byte *pk, byte *sk)
{
  int i;

  IndCpaKemKeypair(pk, sk); // sk[0:INDCPA_SECRETKEYBYTES-1] <-- sk
  for (i = 0; i < INDCPA_PUBLICKEYBYTES; i++)
    sk[i + INDCPA_SECRETKEYBYTES] = pk[i]; // sk[INDCPA_SECRETKEYBYTES:INDCPA_SECRETKEYBYTES+INDCPA_SECRETKEYBYTES-1] <-- pk

  SHA3_256 sha3 = SHA3_256();
  sha3.Update(pk, INDCPA_PUBLICKEYBYTES);
  sha3.Final(sk + SECRETKEYBYTES - 64); // Then hash(pk) is appended.

  RandomBytes(sk + SECRETKEYBYTES - KEYBYTES, KEYBYTES); // Remaining part of sk contains a pseudo-random number.
                                                                           // This is output when check in KemDec() fails.
  return 0;
}


/**
* Generates cipher text c and shared
* secret k for given public key pk
* Length of pk: PUBLICKEYBYTES
* Length of k: SHAREDSECRETBYTES
* Lenght of c: CIPHERTEXTBYTES
*/
template<byte T_L, byte T_Et>
int SABER<T_L, T_Et>::KemEnc(byte *c, byte *k, const byte *pk)
{

  byte kr[64]; // Will contain key, coins
  byte buf[64];

  RandomBytes(buf, 32);

  SHA3_256 sha3_256 = SHA3_256();

  sha3_256.Update(buf, 32); // BUF[0:31] <-- random message (will be used as the key for client) Note: hash doesnot release system RNG output
  sha3_256.Final(buf);

  sha3_256.Update(pk, INDCPA_PUBLICKEYBYTES);
  sha3_256.Final(buf + 32);            // BUF[32:63] <-- Hash(public key);  Multitarget countermeasure for coins + contributory KEM

  SHA3_512 sha3_512 = SHA3_512();

  sha3_512.Update(buf, 64);            // kr[0:63] <-- Hash(buf[0:63]);
  sha3_512.Final(kr);                  // K^ <-- kr[0:31]
                                       // noiseseed (r) <-- kr[32:63];
  IndCpaKemEnc(buf, kr + 32, pk, c);   // buf[0:31] contains message; kr[32:63] contains randomness r;

  sha3_256.Update(c, CCA_DEC);
  sha3_256.Final(kr+32);

  sha3_256.Update(kr, 64); // hash concatenation of pre-k and h(c) to k
  sha3_256.Final(k);

  return 0;
}


/*
* Generates shared secret for given
* cipher text and private key
* On failure, ss will contain a pseudo-random value.
* Length of sk: SECRETKEYBYTES
* Length of c: CIPHERTEXTBYTES
* Length of k: SHAREDSECRETBYTES
*/
template<byte T_L, byte T_Et>
int SABER<T_L, T_Et>::KemDec(byte *k, const byte *c, const byte *sk)
{
  int i;
  bool fail;
  byte cmp[CCA_DEC];
  byte buf[64];
  byte kr[64]; // Will contain key, coins
  const byte *pk = sk + INDCPA_SECRETKEYBYTES;

  IndCpaKemDec(sk, c, buf); // buf[0:31] <-- message

  // Multitarget countermeasure for coins + contributory KEM
  for (i = 0; i < 32; i++) // Save hash by storing h(pk) in sk
    buf[32 + i] = sk[SECRETKEYBYTES - 64 + i];

  SHA3_512 sha3_512 = SHA3_512(); 

  sha3_512.Update(buf, 64);
  sha3_512.Final(kr);

  IndCpaKemEnc(buf, kr + 32, pk, cmp);

  fail = !VerifyBufsEqual(c, cmp, CCA_DEC);

  SHA3_256 sha3_256 = SHA3_256();

  sha3_256.Update(c, CCA_DEC); // overwrite coins in kr with h(c)
  sha3_256.Final(kr + 32);

  Cmov(kr, sk + SECRETKEYBYTES - KEYBYTES, KEYBYTES, fail);

  sha3_256.Update(kr, 64); // hash concatenation of pre-k and h(c) to k
  sha3_256.Final(k);

  return (0);
}

/**
 * Generates public and private key for the CPA-secure
 * public-key encryption scheme underlying Saber
 */
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::IndCpaKemKeypair(byte *pk, byte *sk)
{
	word16 A[L][L][N];
	word16 s[L][N];
	word16 b[L][N] = {0};

	byte seedA[SEEDBYTES];
    byte seedS[NOISE_SEEDBYTES];
	int i, j;

	RandomBytes(seedA, SEEDBYTES);
    SHAKE128 shake128 = SHAKE128(SEEDBYTES);

    shake128.Update(seedA, SEEDBYTES);
	shake128.Final(seedA); // for not revealing system RNG state
	RandomBytes(seedS, NOISE_SEEDBYTES);

	GenMatrix(A, seedA);
	GenSecret(s, seedS);
	MatrixVectorMul(A, s, b, 1);

	for (i = 0; i < L; i++)
	{
		for (j = 0; j < N; j++)
		{
			b[i][j] = (b[i][j] + h1) >> (EQ - EP);
		}
	}

	POLVECq2BS(sk, s);
	POLVECp2BS(pk, b);
	std::memcpy(pk + POLYVECCOMPRESSEDBYTES, seedA, sizeof(seedA));
}
/**
* Encryption function of the CPA-secure
* public-key encryption scheme underlying Saber.
*/
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::IndCpaKemEnc(const byte *m, const byte *seedSp, const byte *pk, byte *cipherText)
{
	word16 A[L][L][N];
	word16 sp[L][N];
	word16 bp[L][N] = {0};
	word16 vp[N] = {0};
	word16 mp[N];
	word16 b[L][N];
	int i, j;
	const byte *seedA = pk + POLYVECCOMPRESSEDBYTES;
 
	GenMatrix(A, seedA);
	GenSecret(sp, seedSp);
	MatrixVectorMul(A, sp, bp, 0);

	for (i = 0; i < L; i++)
	{
		for (j = 0; j < N; j++)
		{
			bp[i][j] = (bp[i][j] + h1) >> (EQ - EP);
		}
	}

	POLVECp2BS(cipherText, bp);
	BS2POLVECp(pk, b);
	InnerProd(b, sp, vp);

	BS2POLmsg(m, mp);

	for (j = 0; j < N; j++)
	{
		vp[j] = (vp[j] - (mp[j] << (EP - 1)) + h1) >> (EP - ET);
	}

	POLT2BS(cipherText + POLYVECCOMPRESSEDBYTES, vp);
}

/**
 * Decryption function of the CPA-secure
 * public-key encryption scheme underlying Saber.
 */
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::IndCpaKemDec(const byte *sk, const byte *cipherText, byte *m)
{

	word16 s[L][N];
	word16 b[L][N];
	word16 v[N] = {0};
	word16 cm[N];
	int i;

	BS2POLVECq(sk, s);
	BS2POLVECp(cipherText, b);
	InnerProd(b, s, v);
	BS2POLT(cipherText + POLYVECCOMPRESSEDBYTES, cm);

	for (i = 0; i < N; i++)
	{
		v[i] = (v[i] + h2 - (cm[i] << (EP - ET))) >> (EP - 1);
	}

	POLmsg2BS(m, v);
}


template class SABER<2,3>;
template class SABER<3,4>;
template class SABER<4,6>;



NAMESPACE_END