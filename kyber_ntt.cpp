/** The C++ file including the NTT stuff for CRYSTALS-Kyber. 
Adapted from the reference
implementation of kyber by the CRYSTALS team 
(https://github.com/pq-crystals/kyber) */

#include "pch.h"
#include "kyber.h"

NAMESPACE_BEGIN(CryptoPP)

template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
const sword16 Kyber<T_K, T_C, T_C2, T_ETA1>::zetas[128] = {
  -1044,  -758,  -359, -1517,  1493,  1422,   287,   202,
  -171,   622,  1577,   182,   962, -1202, -1474,  1468,
  573, -1325,   264,   383,  -829,  1458, -1602,  -130,
  -681,  1017,   732,   608, -1542,   411,  -205, -1571,
  1223,   652,  -552,  1015, -1293,  1491,  -282, -1544,
  516,    -8,  -320,  -666, -1618, -1162,   126,  1469,
  -853,   -90,  -271,   830,   107, -1421,  -247,  -951,
  -398,   961, -1508,  -725,   448, -1065,   677, -1275,
  -1103,   430,   555,   843, -1251,   871,  1550,   105,
  422,   587,   177,  -235,  -291,  -460,  1574,  1653,
  -246,   778,  1159,  -147,  -777,  1483,  -602,  1119,
  -1590,   644,  -872,   349,   418,   329,  -156,   -75,
  817,  1097,   603,   610,  1322, -1285, -1465,   384,
  -1215,  -136,  1218, -1335,  -874,   220, -1187, -1659,
  -1185, -1530, -1278,   794, -1510,  -854,  -870,   478,
  -108,  -308,   996,   991,   958, -1460,  1522,  1628
};



//Multiplication followed by Montgomery reduction
//Returns 16-bit integer congruent to a*b*R^{-1} mod q
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
sword16 Kyber<T_K, T_C, T_C2, T_ETA1>::Fqmul(sword16 a, sword16 b) {
  return MontgomeryReduce((sword32)a*b);
}

//Inplace number-theoretic transform (NTT) in Rq
//input is in standard order, output is in bitreversed order
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::Ntt(sword16 r[256]) {
  word32 len, start, j, k;
  sword16 t, zeta;

  k = 1;
  for(len = 128; len >= 2; len >>= 1) {
    for(start = 0; start < 256; start = j + len) {
      zeta = zetas[k++];
      for(j = start; j < start + len; j++) {
        t = Fqmul(zeta, r[j + len]);
        r[j + len] = r[j] - t;
        r[j] = r[j] + t;
      }
    }
  }
}

//Inplace inverse number-theoretic transform in Rq and
//multiplication by Montgomery factor 2^16.
//Input is in bitreversed order, output is in standard order
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::InvNtt(sword16 r[256]) {
  word32 start, len, j, k;
  sword16 t, zeta;
  const sword16 f = 1441; // mont^2/128

  k = 127;
  for(len = 2; len <= 128; len <<= 1) {
    for(start = 0; start < 256; start = j + len) {
      zeta = zetas[k--];
      for(j = start; j < start + len; ++j) {
        t = r[j];
        r[j] = BarrettReduce(t + r[j + len]);
        r[j + len] = r[j + len] - t;
        r[j + len] = Fqmul(zeta, r[j + len]);
      }
    }
  }

  for(j = 0; j < 256; ++j)
    r[j] = Fqmul(r[j], f);
}
//sword16 r[2]:       pointer to the output polynomial
//const sword16 b[2]: pointer to the second factor
//sword16 zeta:       integer defining the reduction polynomial
//const sword16 a[2]: pointer to the first factor
template<sword32 T_K, word32 T_C,  word32 T_C2, byte T_ETA1>
void Kyber<T_K, T_C, T_C2, T_ETA1>::Basemul(sword16 r[2],
             const sword16 a[2],
             const sword16 b[2],
             sword16 zeta)
{
  r[0]  = Fqmul(a[1], b[1]);
  r[0]  = Fqmul(r[0], zeta);
  r[0] += Fqmul(a[0], b[0]);

  r[1]  = Fqmul(a[0], b[1]);
  r[1] += Fqmul(a[1], b[0]);
}

template class Kyber<2, 320, 128, 3>;
template class Kyber<3, 320, 128, 2>;
template class Kyber<4, 352, 160, 2>;

NAMESPACE_END
