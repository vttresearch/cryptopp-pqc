/*
 * Utility functions for Kyber. Adapted from the reference
 * implementation of Kyber by the CRYSTALS team
 * (https://github.com/pq-crystals/kyber)
 */


#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include "kyber_utils.h"
#include "kyber.h"
#include "shake.h"


#ifdef _WIN32
#include <windows.h>
#include <wincrypt.h>
#else
#include <fcntl.h>
#include <errno.h>
#ifdef __linux__

#include <unistd.h>
#include <sys/syscall.h>
#else
#include <unistd.h>
#endif
#endif


NAMESPACE_BEGIN(CryptoPP)

// The Keccak core function
extern void KeccakF1600(word64 *state);

#ifdef _WIN32
void RandomBytes::randombytes(uint8_t *out, size_t outlen) {
  HCRYPTPROV ctx;
  DWORD len;

  if(!CryptAcquireContext(&ctx, NULL, NULL, PROV_RSA_FULL, CRYPT_VERIFYCONTEXT))
    abort();

  while(outlen > 0) {
    len = (outlen > 1048576) ? 1048576 : outlen;
    if(!CryptGenRandom(ctx, len, (BYTE *)out))
      abort();

    out += len;
    outlen -= len;
  }

  if(!CryptReleaseContext(ctx, 0))
    abort();
}
#elif defined(__linux__) && defined(SYS_getrandom)
void RandomBytes::randombytes(uint8_t *out, size_t outlen) {
  ssize_t ret;

  while(outlen > 0) {
    ret = syscall(SYS_getrandom, out, outlen, 0);
    if(ret == -1 && errno == EINTR)
      continue;
    else if(ret == -1)
      abort();

    out += ret;
    outlen -= ret;
  }
}
#else
void RandomBytes::randombytes(uint8_t *out, size_t outlen) {
  static int fd = -1;
  ssize_t ret;

  while(fd == -1) {
    fd = open("/dev/urandom", O_RDONLY);
    if(fd == -1 && errno == EINTR)
      continue;
    else if(fd == -1)
      abort();
  }

  while(outlen > 0) {
    ret = read(fd, out, outlen);
    if(ret == -1 && errno == EINTR)
      continue;
    else if(ret == -1)
      abort();

    out += ret;
    outlen -= ret;
  }
}
#endif

/*************************************************
* Name:        load64
*
* Description: Load 8 bytes into uint64_t in little-endian order
*
* Arguments:   - const uint8_t *x: pointer to input byte array
*
* Returns the loaded 64-bit unsigned integer
**************************************************/
static uint64_t Load64(const uint8_t x[8]) {
  unsigned int i;
  uint64_t r = 0;

  for(i=0;i<8;i++)
    r |= (uint64_t)x[i] << 8*i;

  return r;
}

/*************************************************
* Name:        store64
*
* Description: Store a 64-bit integer to array of 8 bytes in little-endian order
*
* Arguments:   - uint8_t *x: pointer to the output byte array (allocated)
*              - uint64_t u: input 64-bit unsigned integer
**************************************************/
static void Store64(uint8_t x[8], uint64_t u) {
  unsigned int i;

  for(i=0;i<8;i++)
    x[i] = u >> 8*i;
}


/*************************************************
* Name:        keccak_absorb
*
* Description: Absorb step of Keccak;
*              non-incremental, starts by zeroeing the state.
*
* Arguments:   - uint64_t *s: pointer to (uninitialized) output Keccak state
*              - unsigned int r: rate in bytes (e.g., 168 for SHAKE128)
*              - const uint8_t *m: pointer to input to be absorbed into s
*              - size_t mlen: length of input in bytes
*              - uint8_t p: domain-separation byte for different
*                           Keccak-derived functions
**************************************************/
void Kyber::KeccakAbsorb(uint64_t s[25], unsigned int r, const uint8_t *m, int mlen, uint8_t p)
{
  size_t i;
  uint8_t t[200] = {0};

  /* Zero state */
  for(i=0;i<25;i++)
    s[i] = 0;

  while(mlen >= r) {
    for(i=0;i<r/8;i++)
      s[i] ^= Load64(m + 8*i);

    KeccakF1600(s);
    mlen -= r;
    m += r;
  }

  for(i=0;i<mlen;i++)
    t[i] = m[i];
  t[i] = p;
  t[r-1] |= 128;
  for(i=0;i<r/8;i++)
    s[i] ^= Load64(t + 8*i);
}

/*************************************************
* Name:        keccak_squeezeblocks
*
* Description: Squeeze step of Keccak. Squeezes full blocks of r bytes each.
*              Modifies the state. Can be called multiple times to keep
*              squeezing, i.e., is incremental.
*
* Arguments:   - uint8_t *h: pointer to output blocks
*              - size_t nblocks: number of blocks to be squeezed (written to h)
*              - uint64_t *s: pointer to input/output Keccak state
*              - unsigned int r: rate in bytes (e.g., 168 for SHAKE128)
**************************************************/
void Kyber::KeccakSqueezeBlocks(uint8_t *out,
                                 size_t nblocks,
                                 uint64_t s[25],
                                 unsigned int r)
{
  unsigned int i;
  while(nblocks > 0) {
    KeccakF1600(s);
    for(i=0;i<r/8;i++)
      Store64(out + 8*i, s[i]);
    out += r;
    --nblocks;
  }
}

void Kyber::Shake128Absorb(keccak_state *state, const uint8_t *in, int inlen) {
  KeccakAbsorb(state->s, 168, in, inlen, 0x1F);
}



void Kyber::Shake128SqueezeBlocks(uint8_t *out, size_t nblocks, keccak_state *state) {
  KeccakSqueezeBlocks(out, nblocks, state->s, 168);
}


/*************************************************
* Name:        montgomery_reduce
*
* Description: Montgomery reduction; given a 32-bit integer a, computes
*              16-bit integer congruent to a * R^-1 mod q,
*              where R=2^16
*
* Arguments:   - int32_t a: input integer to be reduced;
*                           has to be in {-q2^15,...,q2^15-1}
*
* Returns:     integer in {-q+1,...,q-1} congruent to a * R^-1 modulo q.
**************************************************/
int16_t Kyber::MontgomeryReduce(int32_t a)
{
  int32_t t;
  int16_t u;

  u = a*QINV;
  t = (int32_t)u*KYBER_Q;
  t = a - t;
  t >>= 16;
  return t;
}

/*************************************************
* Name:        barrett_reduce
*
* Description: Barrett reduction; given a 16-bit integer a, computes
*              16-bit integer congruent to a mod q in {0,...,q}
*
* Arguments:   - int16_t a: input integer to be reduced
*
* Returns:     integer in {0,...,q} congruent to a modulo q.
**************************************************/
int16_t Kyber::BarrettReduce(int16_t a) {
  int16_t t;
  const int16_t v = ((1U << 26) + KYBER_Q/2)/KYBER_Q;

  t  = (int32_t)v*a >> 26;
  t *= KYBER_Q;
  return a - t;
}

/*************************************************
* Name:        csubq
*
* Description: Conditionallly subtract q
*
* Arguments:   - int16_t x: input integer
*
* Returns:     a - q if a >= q, else a
**************************************************/
int16_t Kyber::Csubq(int16_t a) {
  a -= KYBER_Q;
  a += (a >> 15) & KYBER_Q;
  return a;
};


/*************************************************
* Name:        kyber_shake256_prf
*
* Description: Usage of SHAKE256 as a PRF, concatenates secret and public input
*              and then generates outlen bytes of SHAKE256 output
*
* Arguments:   - uint8_t *out:       pointer to output
*              - size_t outlen:      number of requested output bytes
*              - const uint8_t *key: pointer to the key
*                                    (of length KYBER_SYMBYTES)
*              - uint8_t nonce:      single-byte nonce (public PRF input)
**************************************************/
void Kyber::Shake256Prf(uint8_t *out,
                        size_t outlen,
                        const uint8_t key[KYBER_SYMBYTES],
                        uint8_t nonce)
{
  unsigned int i;
  uint8_t extkey[KYBER_SYMBYTES+1];

  for(i=0;i<KYBER_SYMBYTES;i++)
    extkey[i] = key[i];
  extkey[i] = nonce;

  SHAKE256 shake = SHAKE256(outlen);
  shake.Update(extkey, sizeof(extkey));
  shake.Final(out);
  //shake256(out, outlen, extkey, sizeof(extkey));
}


/*************************************************
* Name:        load32_littleendian
*
* Description: load bytes into a 32-bit integer
*              in little-endian order
*
* Arguments:   - const uint8_t *x: pointer to input byte array
*
* Returns 32-bit unsigned integer loaded from x
**************************************************/
static uint32_t load32_littleendian(const uint8_t x[4])
{
  uint32_t r;
  r  = (uint32_t)x[0];
  r |= (uint32_t)x[1] << 8;
  r |= (uint32_t)x[2] << 16;
  r |= (uint32_t)x[3] << 24;
  return r;
}

/*************************************************
* Name:        cbd
*
* Description: Given an array of uniformly random bytes, compute
*              polynomial with coefficients distributed according to
*              a centered binomial distribution with parameter KYBER_ETA
*
* Arguments:   - poly *r:            pointer to output polynomial
*              - const uint8_t *buf: pointer to input byte array
**************************************************/
void Kyber::Cbd(poly *r, const uint8_t buf[KYBER_ETA*KYBER_N/4])
{
#if KYBER_ETA != 2
#error "poly_getnoise in poly.c only supports eta=2"
#endif
  unsigned int i,j;
  uint32_t t,d;
  int16_t a,b;

  for(i=0;i<KYBER_N/8;i++) {
    t  = load32_littleendian(buf+4*i);
    d  = t & 0x55555555;
    d += (t>>1) & 0x55555555;

    for(j=0;j<8;j++) {
      a = (d >> (4*j+0)) & 0x3;
      b = (d >> (4*j+2)) & 0x3;
      r->coeffs[8*i+j] = a - b;
    }
  }
}


NAMESPACE_END