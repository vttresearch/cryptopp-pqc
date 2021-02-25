/* The utility file for the post quantum KEM algorithm SABER. Adapted by Julius Hekkala 
from the public-domain reference
implementation of SABER by the SABER team 
(https://github.com/KULeuven-COSIC/SABER) 
 */

#include "pch.h"
#include "saber.h"
#include "osrng.h"

NAMESPACE_BEGIN(CryptoPP)


/**
 * Random bytes implementation. Uses AutoSeededRandomPool
 * of Crypto++.
 */
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::RandomBytes(byte *out, word64 outLen)
{
    AutoSeededRandomPool rng;
    rng.GenerateBlock(out, outLen);
}


/* b = 1 means mov, b = 0 means don't mov*/
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::Cmov(byte *r, const byte *x, size_t len, byte b)
{
  size_t i;

  b = -b;
  for (i = 0; i < len; i++)
    r[i] ^= b & (x[i] ^ r[i]);
}

/**
 * T = 2 ^ ET
 * Takes a polynomial from R_T and transforms it into a bytestring of length ET × 256/8
 */
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::POLT2BS(byte *bytes, const word16 *data)
{
	size_t j, offsetByte, offsetData;
  if (ET == 3) {
    for (j = 0; j < N / 8; j++)
    {
      offsetByte = 3 * j;
      offsetData = 8 * j;
      bytes[offsetByte + 0] = (data[offsetData + 0] & 0x7) | ((data[offsetData + 1] & 0x7) << 3) | ((data[offsetData + 2] & 0x3) << 6);
      bytes[offsetByte + 1] = ((data[offsetData + 2] >> 2) & 0x01) | ((data[offsetData + 3] & 0x7) << 1) | ((data[offsetData + 4] & 0x7) << 4) | (((data[offsetData + 5]) & 0x01) << 7);
      bytes[offsetByte + 2] = ((data[offsetData + 5] >> 1) & 0x03) | ((data[offsetData + 6] & 0x7) << 2) | ((data[offsetData + 7] & 0x7) << 5);
    }

  }
    
  else if (ET == 4) {
    for (j = 0; j < N / 2; j++)
    {
      offsetByte = j;
      offsetData = 2 * j;
      bytes[offsetByte] = (data[offsetData] & 0x0f) | ((data[offsetData + 1] & 0x0f) << 4);
    }

  }
    
  else if (ET == 6) {
    for (j = 0; j < N / 4; j++)
    {
      offsetByte = 3 * j;
      offsetData = 4 * j;
      bytes[offsetByte + 0] = (data[offsetData + 0] & 0x3f) | ((data[offsetData + 1] & 0x03) << 6);
      bytes[offsetByte + 1] = ((data[offsetData + 1] >> 2) & 0x0f) | ((data[offsetData + 2] & 0x0f) << 4);
      bytes[offsetByte + 2] = ((data[offsetData + 2] >> 4) & 0x03) | ((data[offsetData + 3] & 0x3f) << 2);
    }

  }
	
}

/**
 * T = 2 ^ ET
 * Takes a byte string of length ET × 256/8 and
 * transforms it into a polynomial in R_T
 */ 
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::BS2POLT(const byte *bytes, word16 *data)
{
	size_t j, offsetByte, offsetData;
	if (ET == 3) {
		for (j = 0; j < N / 8; j++)
		{
			offsetByte = 3 * j;
			offsetData = 8 * j;
			data[offsetData + 0] = (bytes[offsetByte + 0]) & 0x07;
			data[offsetData + 1] = ((bytes[offsetByte + 0]) >> 3) & 0x07;
			data[offsetData + 2] = (((bytes[offsetByte + 0]) >> 6) & 0x03) | (((bytes[offsetByte + 1]) & 0x01) << 2);
			data[offsetData + 3] = ((bytes[offsetByte + 1]) >> 1) & 0x07;
			data[offsetData + 4] = ((bytes[offsetByte + 1]) >> 4) & 0x07;
			data[offsetData + 5] = (((bytes[offsetByte + 1]) >> 7) & 0x01) | (((bytes[offsetByte + 2]) & 0x03) << 1);
			data[offsetData + 6] = ((bytes[offsetByte + 2] >> 2) & 0x07);
			data[offsetData + 7] = ((bytes[offsetByte + 2] >> 5) & 0x07);
		}
	}
		
	else if (ET == 4) {
		for (j = 0; j < N / 2; j++)
		{
			offsetByte = j;
			offsetData = 2 * j;
			data[offsetData] = bytes[offsetByte] & 0x0f;
			data[offsetData + 1] = (bytes[offsetByte] >> 4) & 0x0f;
		}

	}
		
	else if (ET == 6) {
		for (j = 0; j < N / 4; j++)
		{
			offsetByte = 3 * j;
			offsetData = 4 * j;
			data[offsetData + 0] = bytes[offsetByte + 0] & 0x3f;
			data[offsetData + 1] = ((bytes[offsetByte + 0] >> 6) & 0x03) | ((bytes[offsetByte + 1] & 0x0f) << 2);
			data[offsetData + 2] = ((bytes[offsetByte + 1] & 0xff) >> 4) | ((bytes[offsetByte + 2] & 0x03) << 4);
			data[offsetData + 3] = ((bytes[offsetByte + 2] & 0xff) >> 2);
		}
	}
		
	else {
		throw Exception(Exception::OTHER_ERROR, "SABER: Unsupported SABER parameter.");
	}
}

/**
 * q = 2 ^ EQ
 * Takes a polynomial from R_q and transforms it into a bytestring of length EQ × 256/8.
 */  
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::POLq2BS(byte *bytes, const word16 *data)
{
	size_t j, offsetByte, offsetData;
	for (j = 0; j < N / 8; j++)
	{
		offsetByte = 13 * j;
		offsetData = 8 * j;
		bytes[offsetByte + 0] = (data[offsetData + 0] & (0xff));
		bytes[offsetByte + 1] = ((data[offsetData + 0] >> 8) & 0x1f) | ((data[offsetData + 1] & 0x07) << 5);
		bytes[offsetByte + 2] = ((data[offsetData + 1] >> 3) & 0xff);
		bytes[offsetByte + 3] = ((data[offsetData + 1] >> 11) & 0x03) | ((data[offsetData + 2] & 0x3f) << 2);
		bytes[offsetByte + 4] = ((data[offsetData + 2] >> 6) & 0x7f) | ((data[offsetData + 3] & 0x01) << 7);
		bytes[offsetByte + 5] = ((data[offsetData + 3] >> 1) & 0xff);
		bytes[offsetByte + 6] = ((data[offsetData + 3] >> 9) & 0x0f) | ((data[offsetData + 4] & 0x0f) << 4);
		bytes[offsetByte + 7] = ((data[offsetData + 4] >> 4) & 0xff);
		bytes[offsetByte + 8] = ((data[offsetData + 4] >> 12) & 0x01) | ((data[offsetData + 5] & 0x7f) << 1);
		bytes[offsetByte + 9] = ((data[offsetData + 5] >> 7) & 0x3f) | ((data[offsetData + 6] & 0x03) << 6);
		bytes[offsetByte + 10] = ((data[offsetData + 6] >> 2) & 0xff);
		bytes[offsetByte + 11] = ((data[offsetData + 6] >> 10) & 0x07) | ((data[offsetData + 7] & 0x1f) << 3);
		bytes[offsetByte + 12] = ((data[offsetData + 7] >> 5) & 0xff);
	}
}

/**
 * q = 2 ^ EQ
 * Takes a byte string of length EQ × 256/8 and
 * transforms it into a polynomial in R_q
 */
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::BS2POLq(const byte *bytes, word16 *data)
{
	size_t j, offsetByte, offsetData;
	for (j = 0; j < N / 8; j++)
	{
		offsetByte = 13 * j;
		offsetData = 8 * j;
		data[offsetData + 0] = (bytes[offsetByte + 0] & (0xff)) | ((bytes[offsetByte + 1] & 0x1f) << 8);
		data[offsetData + 1] = (bytes[offsetByte + 1] >> 5 & (0x07)) | ((bytes[offsetByte + 2] & 0xff) << 3) | ((bytes[offsetByte + 3] & 0x03) << 11);
		data[offsetData + 2] = (bytes[offsetByte + 3] >> 2 & (0x3f)) | ((bytes[offsetByte + 4] & 0x7f) << 6);
		data[offsetData + 3] = (bytes[offsetByte + 4] >> 7 & (0x01)) | ((bytes[offsetByte + 5] & 0xff) << 1) | ((bytes[offsetByte + 6] & 0x0f) << 9);
		data[offsetData + 4] = (bytes[offsetByte + 6] >> 4 & (0x0f)) | ((bytes[offsetByte + 7] & 0xff) << 4) | ((bytes[offsetByte + 8] & 0x01) << 12);
		data[offsetData + 5] = (bytes[offsetByte + 8] >> 1 & (0x7f)) | ((bytes[offsetByte + 9] & 0x3f) << 7);
		data[offsetData + 6] = (bytes[offsetByte + 9] >> 6 & (0x03)) | ((bytes[offsetByte + 10] & 0xff) << 2) | ((bytes[offsetByte + 11] & 0x07) << 10);
		data[offsetData + 7] = (bytes[offsetByte + 11] >> 3 & (0x1f)) | ((bytes[offsetByte + 12] & 0xff) << 5);
	}
}

/**
 * p = 2 ^ EP
 * Takes a polynomial from R_p and transforms it into a bytestring of length EP × 256/8.
 */
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::POLp2BS(byte *bytes, const word16 *data)
{
	size_t j, offsetByte, offsetData;
	for (j = 0; j < N / 4; j++)
	{
		offsetByte = 5 * j;
		offsetData = 4 * j;
		bytes[offsetByte + 0] = (data[offsetData + 0] & (0xff));
		bytes[offsetByte + 1] = ((data[offsetData + 0] >> 8) & 0x03) | ((data[offsetData + 1] & 0x3f) << 2);
		bytes[offsetByte + 2] = ((data[offsetData + 1] >> 6) & 0x0f) | ((data[offsetData + 2] & 0x0f) << 4);
		bytes[offsetByte + 3] = ((data[offsetData + 2] >> 4) & 0x3f) | ((data[offsetData + 3] & 0x03) << 6);
		bytes[offsetByte + 4] = ((data[offsetData + 3] >> 2) & 0xff);
	}
}


/**
 * p = 2 ^ EP
 * Takes a byte string of length EP × 256/8 and
 * transforms it into a polynomial in R_p
 */
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::BS2POLp(const byte *bytes, word16 *data)
{
	size_t j, offsetByte, offsetData;
	for (j = 0; j < N / 4; j++)
	{
		offsetByte = 5 * j;
		offsetData = 4 * j;
		data[offsetData + 0] = (bytes[offsetByte + 0] & (0xff)) | ((bytes[offsetByte + 1] & 0x03) << 8);
		data[offsetData + 1] = ((bytes[offsetByte + 1] >> 2) & (0x3f)) | ((bytes[offsetByte + 2] & 0x0f) << 6);
		data[offsetData + 2] = ((bytes[offsetByte + 2] >> 4) & (0x0f)) | ((bytes[offsetByte + 3] & 0x3f) << 4);
		data[offsetData + 3] = ((bytes[offsetByte + 3] >> 6) & (0x03)) | ((bytes[offsetByte + 4] & 0xff) << 2);
	}
}

/**
 * q = 2 ^ EQ
 * Takes a vector from R_q ^ (L x 1). Transforms it into a byte string of length
 * L * EQ * 256 / 8.
 */ 
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::POLVECq2BS(byte bytes[POLYVECBYTES], const word16 data[L][N])
{
	size_t i;
	for (i = 0; i < L; i++)
	{
		POLq2BS(bytes + i * POLYBYTES, data[i]);
	}
}

/* q = 2 ^ EQ
 * Takes a byte string of length L * EQ * 256 / 8. Transforms it into a vector from R_q ^ (L x 1)
 *
 */
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::BS2POLVECq(const byte bytes[POLYVECBYTES], word16 data[L][N])
{
	size_t i;
	for (i = 0; i < L; i++)
	{
		BS2POLq(bytes + i * POLYBYTES, data[i]);
	}
}

/**
 * p = 2 ^ EP
 * Takes a vector from R_p ^ (L x 1). Transforms it into a byte string of length
 * L * EP * 256 / 8.
 */ 
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::POLVECp2BS(byte bytes[POLYVECCOMPRESSEDBYTES], const word16 data[L][N])
{
	size_t i;
	for (i = 0; i < L; i++)
	{
		POLp2BS(bytes + i * (EP * N / 8), data[i]);
	}
}

/*
 * p = 2 ^ EP
 * Takes a byte string of length L * EP * 256 / 8. Transforms it into a vector from R_p ^ (L x 1)
 *
 */
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::BS2POLVECp(const byte bytes[POLYVECCOMPRESSEDBYTES], word16 data[L][N])
{
	size_t i;
	for (i = 0; i < L; i++)
	{
		BS2POLp(bytes + i * (EP * N / 8), data[i]);
	}
}

/**
 * Takes a the message byte string and
 * transforms it into a polynomial.
 */
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::BS2POLmsg(const byte bytes[KEYBYTES], word16 data[N])
{
	size_t i, j;
	for (j = 0; j < KEYBYTES; j++)
	{
		for (i = 0; i < 8; i++)
		{
			data[j * 8 + i] = ((bytes[j] >> i) & 0x01);
		}
	}
}

/**
 * Takes the message polynomial and trasnforms it back into 
 * a byte string.
 */
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::POLmsg2BS(byte bytes[KEYBYTES], const word16 data[N])
{
	size_t i, j;
	std::memset(bytes, 0, KEYBYTES);

	for (j = 0; j < KEYBYTES; j++)
	{
		for (i = 0; i < 8; i++)
		{
			bytes[j] = bytes[j] | ((data[j * 8 + i] & 0x01) << i);
		}
	}
}

static word64 LoadLittleEndian(const byte *x, int bytes)
{
	int i;
	word64 r = x[0];
	for (i = 1; i < bytes; i++)
		r |= (word64)x[i] << (8 * i);
	return r;
}

/**
 * Compute a polynomial with coefficients distributed according to
 * a centered binomial distribution.
 */
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::Cbd(word16 *s, const byte *buf) {
	if (MU == 6) {
		word32 t, d, a[4], b[4];
		int i, j;

		for (i = 0; i < N / 4; i++)
		{
			t = LoadLittleEndian(buf + 3 * i, 3);
			d = 0;
			for (j = 0; j < 3; j++)
			d += (t >> j) & 0x249249;

			a[0] = d & 0x7;
			b[0] = (d >> 3) & 0x7;
			a[1] = (d >> 6) & 0x7;
			b[1] = (d >> 9) & 0x7;
			a[2] = (d >> 12) & 0x7;
			b[2] = (d >> 15) & 0x7;
			a[3] = (d >> 18) & 0x7;
			b[3] = (d >> 21);

			s[4 * i + 0] = (word16)(a[0] - b[0]);
			s[4 * i + 1] = (word16)(a[1] - b[1]);
			s[4 * i + 2] = (word16)(a[2] - b[2]);
			s[4 * i + 3] = (word16)(a[3] - b[3]);
		}
	}
  
	else if (MU == 8) {
		word32 t, d, a[4], b[4];
		int i, j;

		for (i = 0; i < N / 4; i++)
		{
			t = LoadLittleEndian(buf + 4 * i, 4);
			d = 0;
			for (j = 0; j < 4; j++)
			d += (t >> j) & 0x11111111;

			a[0] = d & 0xf;
			b[0] = (d >> 4) & 0xf;
			a[1] = (d >> 8) & 0xf;
			b[1] = (d >> 12) & 0xf;
			a[2] = (d >> 16) & 0xf;
			b[2] = (d >> 20) & 0xf;
			a[3] = (d >> 24) & 0xf;
			b[3] = (d >> 28);

			s[4 * i + 0] = (word16)(a[0] - b[0]);
			s[4 * i + 1] = (word16)(a[1] - b[1]);
			s[4 * i + 2] = (word16)(a[2] - b[2]);
			s[4 * i + 3] = (word16)(a[3] - b[3]);
		}

	} 
  
	else if (MU == 10) {
		word64 t, d, a[4], b[4];
		int i, j;

		for (i = 0; i < N / 4; i++)
		{
			t = LoadLittleEndian(buf + 5 * i, 5);
			d = 0;
			for (j = 0; j < 5; j++)
			d += (t >> j) & 0x0842108421UL;

			a[0] = d & 0x1f;
			b[0] = (d >> 5) & 0x1f;
			a[1] = (d >> 10) & 0x1f;
			b[1] = (d >> 15) & 0x1f;
			a[2] = (d >> 20) & 0x1f;
			b[2] = (d >> 25) & 0x1f;
			a[3] = (d >> 30) & 0x1f;
			b[3] = (d >> 35);

			s[4 * i + 0] = (word16)(a[0] - b[0]);
			s[4 * i + 1] = (word16)(a[1] - b[1]);
			s[4 * i + 2] = (word16)(a[2] - b[2]);
			s[4 * i + 3] = (word16)(a[3] - b[3]);
		}

	}
 
	else {
		throw Exception(Exception::OTHER_ERROR, "SABER: Unsupported SABER parameter.");
	}
}



template class SABER<2,3>;
template class SABER<3,4>;
template class SABER<4,6>;



NAMESPACE_END