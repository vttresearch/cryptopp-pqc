/** The C++ file for handling the polynomials of 
CRYSTALS-Kyber. Adapted from the reference
implementation of kyber by the CRYSTALS team 
(https://github.com/pq-crystals/kyber) */

#include "pch.h"
#include "kyber.h"

NAMESPACE_BEGIN(CryptoPP)



//Compression and subsequent serialization of a polynomial
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyCompress(uint8_t *r, poly *a)
{
  unsigned int i,j;
  uint8_t t[8];

  PolyCsubq(a);

  if (POLYCOMPRESSEDBYTES == 96) {
    for(i=0;i<N/8;i++) {
      for(j=0;j<8;j++)
        t[j] = ((((uint16_t)a->coeffs[8*i+j] << 3) + Q/2)/Q) & 7;

      r[0] = (t[0] >> 0) | (t[1] << 3) | (t[2] << 6);
      r[1] = (t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7);
      r[2] = (t[5] >> 1) | (t[6] << 2) | (t[7] << 5);
      r += 3;
    }
  }
  else if (POLYCOMPRESSEDBYTES == 128) {
    for(i=0;i<N/8;i++) {
      for(j=0;j<8;j++)
        t[j] = ((((uint16_t)a->coeffs[8*i+j] << 4) + Q/2)/Q) & 15;

      r[0] = t[0] | (t[1] << 4);
      r[1] = t[2] | (t[3] << 4);
      r[2] = t[4] | (t[5] << 4);
      r[3] = t[6] | (t[7] << 4);
      r += 4;
    }
  }
    
  else if (POLYCOMPRESSEDBYTES == 160) {
    for(i=0;i<N/8;i++) {
      for(j=0;j<8;j++)
        t[j] = ((((uint32_t)a->coeffs[8*i+j] << 5) + Q/2)/Q) & 31;

      r[0] = (t[0] >> 0) | (t[1] << 5);
      r[1] = (t[1] >> 3) | (t[2] << 2) | (t[3] << 7);
      r[2] = (t[3] >> 1) | (t[4] << 4);
      r[3] = (t[4] >> 4) | (t[5] << 1) | (t[6] << 6);
      r[4] = (t[6] >> 2) | (t[7] << 3);
      r += 5;
    }
  }
  
}


//De-serialization and subsequent decompression of a polynomial;
//approximate inverse of PolyCompress
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyDecompress(poly *r, const uint8_t *a)
{
  unsigned int i;

  if (POLYCOMPRESSEDBYTES == 96) {
    unsigned int j;
    uint8_t t[8];
    for(i=0;i<N/8;i++) {
      t[0] = (a[0] >> 0);
      t[1] = (a[0] >> 3);
      t[2] = (a[0] >> 6) | (a[1] << 2);
      t[3] = (a[1] >> 1);
      t[4] = (a[1] >> 4);
      t[5] = (a[1] >> 7) | (a[2] << 1);
      t[6] = (a[2] >> 2);
      t[7] = (a[2] >> 5);
      a += 3;

      for(j=0;j<8;j++)
        r->coeffs[8*i+j] = ((uint16_t)(t[j] & 7)*Q + 4) >> 3;
    }
  }
  else if (POLYCOMPRESSEDBYTES == 128) {
    for(i=0;i<N/2;i++) {
      r->coeffs[2*i+0] = (((uint16_t)(a[0] & 15)*Q) + 8) >> 4;
      r->coeffs[2*i+1] = (((uint16_t)(a[0] >> 4)*Q) + 8) >> 4;
      a += 1;
    }
  }
  else if (POLYCOMPRESSEDBYTES == 160) {
    unsigned int j;
    uint8_t t[8];
    for(i=0;i<N/8;i++) {
      t[0] = (a[0] >> 0);
      t[1] = (a[0] >> 5) | (a[1] << 3);
      t[2] = (a[1] >> 2);
      t[3] = (a[1] >> 7) | (a[2] << 1);
      t[4] = (a[2] >> 4) | (a[3] << 4);
      t[5] = (a[3] >> 1);
      t[6] = (a[3] >> 6) | (a[4] << 2);
      t[7] = (a[4] >> 3);
      a += 5;

      for(j=0;j<8;j++)
        r->coeffs[8*i+j] = ((uint32_t)(t[j] & 31)*Q + 16) >> 5;
    }
  }

}


//Serialize the polynomial
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyToBytes(uint8_t *r, poly *a)
{
  unsigned int i;
  uint16_t t0, t1;

  PolyCsubq(a);

  for(i=0;i<N/2;i++) {
    t0 = a->coeffs[2*i];
    t1 = a->coeffs[2*i+1];
    r[3*i+0] = (t0 >> 0);
    r[3*i+1] = (t0 >> 8) | (t1 << 4);
    r[3*i+2] = (t1 >> 4);
  }
}

//Deserialize the polynomial
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyFromBytes(poly* r, const uint8_t *a)
{
  unsigned int i;
  for(i=0;i<N/2;i++) {
    r->coeffs[2*i]   = ((a[3*i+0] >> 0) | ((uint16_t)a[3*i+1] << 8)) & 0xFFF;
    r->coeffs[2*i+1] = ((a[3*i+1] >> 4) | ((uint16_t)a[3*i+2] << 4)) & 0xFFF;
  }
}

//Convert a 32 byte message to polynomial
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyFromMsg(poly *r, const uint8_t *msg)
{
  unsigned int i,j;
  int16_t mask;


  for(i=0;i<N/8;i++) {
    for(j=0;j<8;j++) {
      mask = -(int16_t)((msg[i] >> j)&1);
      r->coeffs[8*i+j] = mask & ((Q+1)/2);
    }
  }
}

//Convert polynomial to a 32 byte message
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyToMsg(uint8_t *msg, poly *a)
{
  unsigned int i,j;
  uint16_t t;

  PolyCsubq(a);

  for(i=0;i<N/8;i++) {
    msg[i] = 0;
    for(j=0;j<8;j++) {
      t = ((((uint16_t)a->coeffs[8*i+j] << 1) + Q/2)/Q) & 1;
      msg[i] |= t << j;
    }
  }
}

//Sample a polynomial deterministically from a seed and a nonce,
//with output polynomial close to centered binomial distribution
//with parameter KYBER_ETA
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyGetNoise(poly *r, const uint8_t *seed, uint8_t nonce)
{
  uint8_t buf[ETA*N/4];
  Shake256Prf(buf, sizeof(buf), seed, nonce);
  Cbd(r, buf);
}

//Computes negacyclic number-theoretic transform (NTT) of
//a polynomial in place;
//inputs assumed to be in normal order, output in bitreversed order
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyNtt(poly *r)
{
  Ntt(r->coeffs);
  PolyReduce(r);
}


//Computes inverse of negacyclic number-theoretic transform (NTT)
//of a polynomial in place;
//inputs assumed to be in bitreversed order, output in normal order
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyInvNttToMont(poly *r)
{
  InvNtt(r->coeffs);
}

//Multiplication of two polynomials in NTT domain
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyBasemulMontgomery(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  for(i=0;i<N/4;i++) {
    Basemul(&r->coeffs[4*i], &a->coeffs[4*i], &b->coeffs[4*i], zetas[64+i]);
    Basemul(&r->coeffs[4*i+2], &a->coeffs[4*i+2], &b->coeffs[4*i+2],
            -zetas[64+i]);
  }
}


//Inplace conversion of all coefficients of a polynomial
//from normal domain to Montgomery domain
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyToMont(poly *r)
{
  unsigned int i;
  const int16_t f = (1ULL << 32) % Q;
  for(i=0;i<N;i++)
    r->coeffs[i] = MontgomeryReduce((int32_t)r->coeffs[i]*f);
}


//Applies Barrett reduction to all coefficients of a polynomial
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyReduce(poly *r)
{
  unsigned int i;
  for(i=0;i<N;i++)
    r->coeffs[i] = BarrettReduce(r->coeffs[i]);
}

// Applies conditional subtraction of q to each coefficient
// of a polynomial.
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyCsubq(poly *r)
{
  unsigned int i;
  for(i=0;i<N;i++)
    r->coeffs[i] = Csubq(r->coeffs[i]);
}


//Add two polynomials

template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyAdd(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  for(i=0;i<N;i++)
    r->coeffs[i] = a->coeffs[i] + b->coeffs[i];
}


//Subtract two polynomials
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolySub(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  for(i=0;i<N;i++)
    r->coeffs[i] = a->coeffs[i] - b->coeffs[i];
}



//Compress and serialize vector of polynomials
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyvecCompress(uint8_t *r, polyvec *a)
{
  unsigned int i,j,k;

  PolyvecCsubq(a);

  if (POLYVECCOMPRESSEDBYTES == (K * 352)) {
    uint16_t t[8];
    for(i=0;i<K;i++) {
      for(j=0;j<N/8;j++) {
        for(k=0;k<8;k++)
          t[k] = ((((uint32_t)a->vec[i].coeffs[8*j+k] << 11) + Q/2)
                  /Q) & 0x7ff;

        r[ 0] = (t[0] >>  0);
        r[ 1] = (t[0] >>  8) | (t[1] << 3);
        r[ 2] = (t[1] >>  5) | (t[2] << 6);
        r[ 3] = (t[2] >>  2);
        r[ 4] = (t[2] >> 10) | (t[3] << 1);
        r[ 5] = (t[3] >>  7) | (t[4] << 4);
        r[ 6] = (t[4] >>  4) | (t[5] << 7);
        r[ 7] = (t[5] >>  1);
        r[ 8] = (t[5] >>  9) | (t[6] << 2);
        r[ 9] = (t[6] >>  6) | (t[7] << 5);
        r[10] = (t[7] >>  3);
        r += 11;
      }
    }
  }
  else if (POLYVECCOMPRESSEDBYTES == (K * 320)) {
    uint16_t t[4];
    for(i=0;i<K;i++) {
      for(j=0;j<N/4;j++) {
        for(k=0;k<4;k++)
          t[k] = ((((uint32_t)a->vec[i].coeffs[4*j+k] << 10) + Q/2)
                  / Q) & 0x3ff;

        r[0] = (t[0] >> 0);
        r[1] = (t[0] >> 8) | (t[1] << 2);
        r[2] = (t[1] >> 6) | (t[2] << 4);
        r[3] = (t[2] >> 4) | (t[3] << 6);
        r[4] = (t[3] >> 2);
        r += 5;
      }
    }
  }

}


//De-serialize and decompress vector of polynomials;
//approximate inverse of polyvec_compress
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyvecDecompress(polyvec *r, const uint8_t *a)
{
  unsigned int i,j,k;

  if (POLYVECCOMPRESSEDBYTES == (K * 352)) {
    uint16_t t[8];
    for(i=0;i<K;i++) {
      for(j=0;j<N/8;j++) {
        t[0] = (a[0] >> 0) | ((uint16_t)a[ 1] << 8);
        t[1] = (a[1] >> 3) | ((uint16_t)a[ 2] << 5);
        t[2] = (a[2] >> 6) | ((uint16_t)a[ 3] << 2) | ((uint16_t)a[4] << 10);
        t[3] = (a[4] >> 1) | ((uint16_t)a[ 5] << 7);
        t[4] = (a[5] >> 4) | ((uint16_t)a[ 6] << 4);
        t[5] = (a[6] >> 7) | ((uint16_t)a[ 7] << 1) | ((uint16_t)a[8] << 9);
        t[6] = (a[8] >> 2) | ((uint16_t)a[ 9] << 6);
        t[7] = (a[9] >> 5) | ((uint16_t)a[10] << 3);
        a += 11;

        for(k=0;k<8;k++)
          r->vec[i].coeffs[8*j+k] = ((uint32_t)(t[k] & 0x7FF)*Q + 1024) >> 11;
      }
    }
  }
  else if (POLYVECCOMPRESSEDBYTES == (K * 320)) {
    uint16_t t[4];
    for(i=0;i<K;i++) {
      for(j=0;j<N/4;j++) {
        t[0] = (a[0] >> 0) | ((uint16_t)a[1] << 8);
        t[1] = (a[1] >> 2) | ((uint16_t)a[2] << 6);
        t[2] = (a[2] >> 4) | ((uint16_t)a[3] << 4);
        t[3] = (a[3] >> 6) | ((uint16_t)a[4] << 2);
        a += 5;

        for(k=0;k<4;k++)
          r->vec[i].coeffs[4*j+k] = ((uint32_t)(t[k] & 0x3FF)*Q + 512) >> 10;
      }
    }
  }
}


//Description: Serialize vector of polynomials
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyvecToBytes(uint8_t *r, polyvec *a)
{
  unsigned int i;
  for(i=0;i<K;i++)
    PolyToBytes(r+i*POLYBYTES, &a->vec[i]);
}

//De-serialize vector of polynomials;
//inverse of polyvec_tobytes

template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyvecFromBytes(polyvec *r, const uint8_t *a)
{
  unsigned int i;
  for(i=0;i<K;i++)
    PolyFromBytes(&r->vec[i], a+i*POLYBYTES);
}

//Apply forward NTT to all elements of a vector of polynomials
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyvecNtt(polyvec *r)
{
  unsigned int i;
  for(i=0;i<K;i++)
    PolyNtt(&r->vec[i]);
}

//Apply inverse NTT to all elements of a vector of polynomials
//and multiply by Montgomery factor 2^16

template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyvecInvNttToMont(polyvec *r)
{
  unsigned int i;
  for(i=0;i<K;i++)
    PolyInvNttToMont(&r->vec[i]);
}


//Pointwise multiply elements of a and b, accumulate into r,
//and multiply by 2^-16.

template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyvecPointwiseAccMontgomery(poly *r, const polyvec *a, const polyvec *b)
{
  unsigned int i;
  poly t;

  PolyBasemulMontgomery(r, &a->vec[0], &b->vec[0]);
  for(i=1;i<K;i++) {
    PolyBasemulMontgomery(&t, &a->vec[i], &b->vec[i]);
    PolyAdd(r, r, &t);
  }

  PolyReduce(r);
}


//Applies Barrett reduction to each coefficient
//of each element of a vector of polynomials
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyvecReduce(polyvec *r)
{
  unsigned int i;
  for(i=0;i<K;i++)
    PolyReduce(&r->vec[i]);
}

//Applies conditional subtraction of q to each coefficient
//of each element of a vector of polynomials
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyvecCsubq(polyvec *r)
{
  unsigned int i;
  for(i=0;i<K;i++)
    PolyCsubq(&r->vec[i]);
}

//Add vectors of polynomials
template<int T_K, unsigned int T_Compr>
void Kyber<T_K, T_Compr>::PolyvecAdd(polyvec *r, const polyvec *a, const polyvec *b)
{
  unsigned int i;
  for(i=0;i<K;i++)
    PolyAdd(&r->vec[i], &a->vec[i], &b->vec[i]);
}





template class Kyber<2, 320>;
template class Kyber<3, 320>;
template class Kyber<4, 352>;



NAMESPACE_END
