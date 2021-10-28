/** The C++ file for handling the polynomials of 
CRYSTALS-Kyber. Adapted from the reference
implementation of kyber by the CRYSTALS team 
(https://github.com/pq-crystals/kyber) */

#include "pch.h"
#include "kyber.h"

NAMESPACE_BEGIN(CryptoPP)



//Compression and subsequent serialization of a polynomial
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyCompress(byte *r, poly *a)
{
  word32 i,j;
  sword16 u;
  byte t[8];

  if (POLYCOMPRESSEDBYTES == 128) {
    for(i=0;i<N/8;i++) {
      for(j=0;j<8;j++) {
        // map to positive standard representatives
        u  = a->coeffs[8*i+j];
        u += (u >> 15) & Q;
        t[j] = ((((word16)u << 4) + Q/2)/Q) & 15;
      }
      r[0] = t[0] | (t[1] << 4);
      r[1] = t[2] | (t[3] << 4);
      r[2] = t[4] | (t[5] << 4);
      r[3] = t[6] | (t[7] << 4);
      r += 4;
    }
  }
    
  else if (POLYCOMPRESSEDBYTES == 160) {
    for(i=0;i<N/8;i++) {
      for(j=0;j<8;j++) {
        // map to positive standard representatives
        u  = a->coeffs[8*i+j];
        u += (u >> 15) & Q;
        t[j] = ((((word32)u << 5) + Q/2)/Q) & 31;
      }
        

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
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyDecompress(poly *r, const byte *a)
{
  word32 i;

  if (POLYCOMPRESSEDBYTES == 128) {
    for(i=0;i<N/2;i++) {
      r->coeffs[2*i+0] = (((word16)(a[0] & 15)*Q) + 8) >> 4;
      r->coeffs[2*i+1] = (((word16)(a[0] >> 4)*Q) + 8) >> 4;
      a += 1;
    }
  }
  else if (POLYCOMPRESSEDBYTES == 160) {
    word32 j;
    byte t[8];
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
        r->coeffs[8*i+j] = ((word32)(t[j] & 31)*Q + 16) >> 5;
    }
  }

}


//Serialize the polynomial
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyToBytes(byte *r, const poly *a)
{
  word32 i;

  word16 t0, t1;

  for(i=0;i<N/2;i++) {
    // map to positive standard representatives
    t0  = a->coeffs[2*i];
    t0 += ((sword16)t0 >> 15) & Q;
    t1 = a->coeffs[2*i+1];
    t1 += ((sword16)t1 >> 15) & Q;
    r[3*i+0] = (t0 >> 0);
    r[3*i+1] = (t0 >> 8) | (t1 << 4);
    r[3*i+2] = (t1 >> 4);
  }
}

//Deserialize the polynomial
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyFromBytes(poly* r, const byte *a)
{
  word32 i;
  for(i=0;i<N/2;i++) {
    r->coeffs[2*i]   = ((a[3*i+0] >> 0) | ((word16)a[3*i+1] << 8)) & 0xFFF;
    r->coeffs[2*i+1] = ((a[3*i+1] >> 4) | ((word16)a[3*i+2] << 4)) & 0xFFF;
  }
}

//Convert a 32 byte message to polynomial
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyFromMsg(poly *r, const byte *msg)
{
  word32 i,j;
  sword16 mask;


  for(i=0;i<N/8;i++) {
    for(j=0;j<8;j++) {
      mask = -(sword16)((msg[i] >> j)&1);
      
      r->coeffs[8*i+j] = mask & ((Q+1)/2);
    }
  }
}

//Convert polynomial to a 32 byte message
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyToMsg(byte *msg, poly *a)
{
  word32 i,j;
  word16 t;

  for(i=0;i<N/8;i++) {
    msg[i] = 0;
    for(j=0;j<8;j++) {
      t  = a->coeffs[8*i+j];
      t += ((sword16)t >> 15) & Q;
      t  = (((t << 1) + Q/2)/Q) & 1;
      msg[i] |= t << j;
    }
  }
}

//Sample a polynomial deterministically from a seed and a nonce,
//with output polynomial close to centered binomial distribution
//with parameter KYBER_ETA
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyGetNoise(poly *r, const byte *seed, byte nonce, const sword32 eta)
{
  std::vector<byte> bufVector(eta*N/4);
  std::vector<byte> *buf = &bufVector;
  KyberShake256PRF(buf->data(), eta*N/4, seed, nonce);

  Cbd(r, buf->data(), eta);
}

//Computes negacyclic number-theoretic transform (NTT) of
//a polynomial in place;
//inputs assumed to be in normal order, output in bitreversed order
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyNtt(poly *r)
{
  Ntt(r->coeffs);
  PolyReduce(r);
}


//Computes inverse of negacyclic number-theoretic transform (NTT)
//of a polynomial in place;
//inputs assumed to be in bitreversed order, output in normal order
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyInvNttToMont(poly *r)
{
  InvNtt(r->coeffs);
}

//Multiplication of two polynomials in NTT domain
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyBasemulMontgomery(poly *r, const poly *a, const poly *b)
{
  word32 i;
  for(i=0;i<N/4;i++) {
    Basemul(&r->coeffs[4*i], &a->coeffs[4*i], &b->coeffs[4*i], zetas[64+i]);
    Basemul(&r->coeffs[4*i+2], &a->coeffs[4*i+2], &b->coeffs[4*i+2],
            -zetas[64+i]);
  }
}


//Inplace conversion of all coefficients of a polynomial
//from normal domain to Montgomery domain
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyToMont(poly *r)
{
  word32 i;
  const sword16 f = (1ULL << 32) % Q;
  for(i=0;i<N;i++)
    r->coeffs[i] = MontgomeryReduce((sword32)r->coeffs[i]*f);
}


//Applies Barrett reduction to all coefficients of a polynomial
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyReduce(poly *r)
{
  word32 i;
  for(i=0;i<N;i++)
    r->coeffs[i] = BarrettReduce(r->coeffs[i]);
}


//Add two polynomials
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyAdd(poly *r, const poly *a, const poly *b)
{
  word32 i;
  for(i=0;i<N;i++)
    r->coeffs[i] = a->coeffs[i] + b->coeffs[i];
}


//Subtract two polynomials
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolySub(poly *r, const poly *a, const poly *b)
{
  word32 i;
  for(i=0;i<N;i++)
    r->coeffs[i] = a->coeffs[i] - b->coeffs[i];
}



//Compress and serialize vector of polynomials
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyvecCompress(byte *r, polyvec *a)
{
  word32 i,j,k;

  if (POLYVECCOMPRESSEDBYTES == (K * 352)) {
    word16 t[8];
    for(i=0;i<K;i++) {
      for(j=0;j<N/8;j++) {
        for(k=0;k<8;k++) {
          t[k]  = a->vec[i].coeffs[8*j+k];
          t[k] += ((sword16)t[k] >> 15) & Q;
          t[k]  = ((((word32)t[k] << 11) + Q/2)/Q) & 0x7ff;
        }
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
    word16 t[4];
    for(i=0;i<K;i++) {
      for(j=0;j<N/4;j++) {
        for(k=0;k<4;k++) {
          t[k]  = a->vec[i].coeffs[4*j+k];
          t[k] += ((sword16)t[k] >> 15) & Q;
          t[k]  = ((((word32)t[k] << 10) + Q/2)/Q) & 0x3ff;
        }

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
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyvecDecompress(polyvec *r, const byte *a)
{
  word32 i,j,k;
  if (POLYVECCOMPRESSEDBYTES == (K * 352)) {

    word16 t[8];
    for(i=0;i<K;i++) {
      for(j=0;j<N/8;j++) {
        t[0] = (a[0] >> 0) | ((word16)a[ 1] << 8);
        t[1] = (a[1] >> 3) | ((word16)a[ 2] << 5);
        t[2] = (a[2] >> 6) | ((word16)a[ 3] << 2) | ((word16)a[4] << 10);
        t[3] = (a[4] >> 1) | ((word16)a[ 5] << 7);
        t[4] = (a[5] >> 4) | ((word16)a[ 6] << 4);
        t[5] = (a[6] >> 7) | ((word16)a[ 7] << 1) | ((word16)a[8] << 9);
        t[6] = (a[8] >> 2) | ((word16)a[ 9] << 6);
        t[7] = (a[9] >> 5) | ((word16)a[10] << 3);
        a += 11;

        for(k=0;k<8;k++)
          r->vec[i].coeffs[8*j+k] = ((word32)(t[k] & 0x7FF)*Q + 1024) >> 11;
      }
    }
  }
  else if (POLYVECCOMPRESSEDBYTES == (K * 320)) {
    word16 t[4];
    for(i=0;i<K;i++) {
      for(j=0;j<N/4;j++) {
        t[0] = (a[0] >> 0) | ((word16)a[1] << 8);
        t[1] = (a[1] >> 2) | ((word16)a[2] << 6);
        t[2] = (a[2] >> 4) | ((word16)a[3] << 4);
        t[3] = (a[3] >> 6) | ((word16)a[4] << 2);
        a += 5;

        for(k=0;k<4;k++)
          r->vec[i].coeffs[4*j+k] = ((word32)(t[k] & 0x3FF)*Q + 512) >> 10;
      }
    }
  }
}


//Description: Serialize vector of polynomials
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyvecToBytes(byte *r, const polyvec *a)
{
  word32 i;
  for(i=0;i<K;i++)
    PolyToBytes(r+i*POLYBYTES, &a->vec[i]);
}

//De-serialize vector of polynomials;
//inverse of polyvec_tobytes

template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyvecFromBytes(polyvec *r, const byte *a)
{
  word32 i;
  for(i=0;i<K;i++)
    PolyFromBytes(&r->vec[i], a+i*POLYBYTES);
}

//Apply forward NTT to all elements of a vector of polynomials
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyvecNtt(polyvec *r)
{
  word32 i;
  for(i=0;i<K;i++)
    PolyNtt(&r->vec[i]);
}

//Apply inverse NTT to all elements of a vector of polynomials
//and multiply by Montgomery factor 2^16

template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyvecInvNttToMont(polyvec *r)
{
  word32 i;
  for(i=0;i<K;i++)
    PolyInvNttToMont(&r->vec[i]);
}


//Multiply elements of a and b in NTT domain, accumulate into r,
//and multiply by 2^-16.
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyvecBasemulAccMontgomery(poly *r, const polyvec *a, const polyvec *b)
{
  word32 i;
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
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyvecReduce(polyvec *r)
{
  word32 i;
  for(i=0;i<K;i++)
    PolyReduce(&r->vec[i]);
}

//Add vectors of polynomials
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::PolyvecAdd(polyvec *r, const polyvec *a, const polyvec *b)
{
  word32 i;
  for(i=0;i<K;i++)
    PolyAdd(&r->vec[i], &a->vec[i], &b->vec[i]);
}




template class Kyber<2, 320, 128, 3>;
template class Kyber<3, 320, 128, 2>;
template class Kyber<4, 352, 160, 2>;



NAMESPACE_END
