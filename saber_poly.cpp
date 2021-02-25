/* The polynomial functions for the post quantum KEM algorithm SABER. Adapted by Julius Hekkala 
from the public-domain reference
implementation of SABER by the SABER team 
(https://github.com/KULeuven-COSIC/SABER) 
 */

#include "pch.h"
#include "saber.h"
#include "shake.h"

NAMESPACE_BEGIN(CryptoPP)



/**
 * Multiplies a matrix A with a vector s. Returns result res.
 * Flag transpose is used to control whether to multiply A or A^T. 
 */
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::MatrixVectorMul(const word16 A[L][L][N], const word16 s[L][N], word16 res[L][N], word16 transpose)
{
	int i, j;
	for (i = 0; i < L; i++)
	{
		for (j = 0; j < L; j++)
		{
			if (transpose == 1)
			{
				PolyMulAcc(A[j][i], s[j], res[i]);
			}
			else
			{
				PolyMulAcc(A[i][j], s[j], res[i]);
			}	
		}
	}
}

/**
 * Calculates the inner product of two vectors. 
 */
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::InnerProd(const word16 b[L][N], const word16 s[L][N], word16 *res)
{
	int j;
	for (j = 0; j < L; j++)
	{
		PolyMulAcc(b[j], s[j], res);
	}
}

/**
 * Generates a matrix from seed of length SEEDBYTES.
 */
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::GenMatrix(word16 A[L][L][N], const byte seed[SEEDBYTES])
{
	byte buf[L * POLYVECBYTES];
	int i;

    SHAKE128 shake128 = SHAKE128(sizeof(buf));
	shake128.Update(seed, SEEDBYTES);
    shake128.Final(buf);

	for (i = 0; i < L; i++)
	{
		BS2POLVECq(buf + i * POLYVECBYTES, A[i]);
	}
}

/**
 * Generates a secret vector with coefficients that are sampled from
 * centered binomial distribution from the input seed of length NOISE_SEEDBYTES.
 */
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::GenSecret(word16 s[L][N], const byte seed[NOISE_SEEDBYTES])
{
	byte buf[L * POLYCOINBYTES];
	size_t i;

    SHAKE128 shake128 = SHAKE128(sizeof(buf));
	shake128.Update(seed, NOISE_SEEDBYTES);
    shake128.Final(buf);

	for (i = 0; i < L; i++)
	{
		Cbd(s[i], buf + i * POLYCOINBYTES);
	}
}

#define SCHB_N 16
#define OVERFLOWING_MUL(X, Y) ((word16)((word32)(X) * (word32)(Y)))
#define KARATSUBA_N 64
template<byte T_L, byte T_Et>
/**
 * Karatsuba algorithm for polynomial multiplication.
 * Provided by the SABER team.
 */
void SABER<T_L, T_Et>::KaratsubaSimple(const word16 *a1, const word16 *b1, word16 *resultFinal) {
    word16 d01[KARATSUBA_N / 2 - 1];
    word16 d0123[KARATSUBA_N / 2 - 1];
    word16 d23[KARATSUBA_N / 2 - 1];
    word16 result_d01[KARATSUBA_N - 1];

    sword32 i, j;

    std::memset(result_d01, 0, (KARATSUBA_N - 1)*sizeof(word16));
    std::memset(d01, 0, (KARATSUBA_N / 2 - 1)*sizeof(word16));
    std::memset(d0123, 0, (KARATSUBA_N / 2 - 1)*sizeof(word16));
    std::memset(d23, 0, (KARATSUBA_N / 2 - 1)*sizeof(word16));
    std::memset(resultFinal, 0, (2 * KARATSUBA_N - 1)*sizeof(word16));

    word16 acc1, acc2, acc3, acc4, acc5, acc6, acc7, acc8, acc9, acc10;


    for (i = 0; i < KARATSUBA_N / 4; i++) {
        acc1 = a1[i]; //a0
        acc2 = a1[i + KARATSUBA_N / 4]; //a1
        acc3 = a1[i + 2 * KARATSUBA_N / 4]; //a2
        acc4 = a1[i + 3 * KARATSUBA_N / 4]; //a3
        for (j = 0; j < KARATSUBA_N / 4; j++) {

            acc5 = b1[j]; //b0
            acc6 = b1[j + KARATSUBA_N / 4]; //b1

            resultFinal[i + j + 0 * KARATSUBA_N / 4] =
                resultFinal[i + j + 0 * KARATSUBA_N / 4] +
                OVERFLOWING_MUL(acc1, acc5);
            resultFinal[i + j + 2 * KARATSUBA_N / 4] =
                resultFinal[i + j + 2 * KARATSUBA_N / 4] +
                OVERFLOWING_MUL(acc2, acc6);

            acc7 = acc5 + acc6; //b01
            acc8 = acc1 + acc2; //a01
            d01[i + j] = d01[i + j] + (word16)(acc7 * (word64)acc8);
            //--------------------------------------------------------

            acc7 = b1[j + 2 * KARATSUBA_N / 4]; //b2
            acc8 = b1[j + 3 * KARATSUBA_N / 4]; //b3
            resultFinal[i + j + 4 * KARATSUBA_N / 4] =
                resultFinal[i + j + 4 * KARATSUBA_N / 4] +
                OVERFLOWING_MUL(acc7, acc3);

            resultFinal[i + j + 6 * KARATSUBA_N / 4] =
                resultFinal[i + j + 6 * KARATSUBA_N / 4] +
                OVERFLOWING_MUL(acc8, acc4);

            acc9 = acc3 + acc4;
            acc10 = acc7 + acc8;
            d23[i + j] = d23[i + j] + OVERFLOWING_MUL(acc9, acc10);
            //--------------------------------------------------------

            acc5 = acc5 + acc7; //b02
            acc7 = acc1 + acc3; //a02
            result_d01[i + j + 0 * KARATSUBA_N / 4] =
                result_d01[i + j + 0 * KARATSUBA_N / 4] +
                OVERFLOWING_MUL(acc5, acc7);

            acc6 = acc6 + acc8; //b13
            acc8 = acc2 + acc4;
            result_d01[i + j + 2 * KARATSUBA_N / 4] =
                result_d01[i + j + 2 * KARATSUBA_N / 4] +
                OVERFLOWING_MUL(acc6, acc8);

            acc5 = acc5 + acc6;
            acc7 = acc7 + acc8;
            d0123[i + j] = d0123[i + j] + OVERFLOWING_MUL(acc5, acc7);
        }
    }

    // 2nd last stage

    for (i = 0; i < KARATSUBA_N / 2 - 1; i++) {
        d0123[i] = d0123[i] - result_d01[i + 0 * KARATSUBA_N / 4] - result_d01[i + 2 * KARATSUBA_N / 4];
        d01[i] = d01[i] - resultFinal[i + 0 * KARATSUBA_N / 4] - resultFinal[i + 2 * KARATSUBA_N / 4];
        d23[i] = d23[i] - resultFinal[i + 4 * KARATSUBA_N / 4] - resultFinal[i + 6 * KARATSUBA_N / 4];
    }

    for (i = 0; i < KARATSUBA_N / 2 - 1; i++) {
        result_d01[i + 1 * KARATSUBA_N / 4] = result_d01[i + 1 * KARATSUBA_N / 4] + d0123[i];
        resultFinal[i + 1 * KARATSUBA_N / 4] = resultFinal[i + 1 * KARATSUBA_N / 4] + d01[i];
        resultFinal[i + 5 * KARATSUBA_N / 4] = resultFinal[i + 5 * KARATSUBA_N / 4] + d23[i];
    }

    // Last stage
    for (i = 0; i < KARATSUBA_N - 1; i++) {
        result_d01[i] = result_d01[i] - resultFinal[i] - resultFinal[i + KARATSUBA_N];
    }

    for (i = 0; i < KARATSUBA_N - 1; i++) {
        resultFinal[i + 1 * KARATSUBA_N / 2] = resultFinal[i + 1 * KARATSUBA_N / 2] + result_d01[i];
    }

}

/**
 * Toom-Cook polynomial multiplication.
 */
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::ToomCook4Way (const word16 *a1, const word16 *b1, word16 *result) {
    word16 inv3 = 43691, inv9 = 36409, inv15 = 61167;

    word16 aw1[N_SB], aw2[N_SB], aw3[N_SB], aw4[N_SB], aw5[N_SB], aw6[N_SB], aw7[N_SB];
    word16 bw1[N_SB], bw2[N_SB], bw3[N_SB], bw4[N_SB], bw5[N_SB], bw6[N_SB], bw7[N_SB];
    word16 w1[N_SB_RES] = {0}, w2[N_SB_RES] = {0}, w3[N_SB_RES] = {0}, w4[N_SB_RES] = {0},
                            w5[N_SB_RES] = {0}, w6[N_SB_RES] = {0}, w7[N_SB_RES] = {0};
    word16 r0, r1, r2, r3, r4, r5, r6, r7;
    word16 *A0, *A1, *A2, *A3, *B0, *B1, *B2, *B3;
    A0 = (word16 *)a1;
    A1 = (word16 *)&a1[N_SB];
    A2 = (word16 *)&a1[2 * N_SB];
    A3 = (word16 *)&a1[3 * N_SB];
    B0 = (word16 *)b1;
    B1 = (word16 *)&b1[N_SB];
    B2 = (word16 *)&b1[2 * N_SB];
    B3 = (word16 *)&b1[3 * N_SB];

    word16 *C;
    C = result;

    int i, j;

    // EVALUATION
    for (j = 0; j < N_SB; ++j) {
        r0 = A0[j];
        r1 = A1[j];
        r2 = A2[j];
        r3 = A3[j];
        r4 = r0 + r2;
        r5 = r1 + r3;
        r6 = r4 + r5;
        r7 = r4 - r5;
        aw3[j] = r6;
        aw4[j] = r7;
        r4 = ((r0 << 2) + r2) << 1;
        r5 = (r1 << 2) + r3;
        r6 = r4 + r5;
        r7 = r4 - r5;
        aw5[j] = r6;
        aw6[j] = r7;
        r4 = (r3 << 3) + (r2 << 2) + (r1 << 1) + r0;
        aw2[j] = r4;
        aw7[j] = r0;
        aw1[j] = r3;
    }
    for (j = 0; j < N_SB; ++j) {
        r0 = B0[j];
        r1 = B1[j];
        r2 = B2[j];
        r3 = B3[j];
        r4 = r0 + r2;
        r5 = r1 + r3;
        r6 = r4 + r5;
        r7 = r4 - r5;
        bw3[j] = r6;
        bw4[j] = r7;
        r4 = ((r0 << 2) + r2) << 1;
        r5 = (r1 << 2) + r3;
        r6 = r4 + r5;
        r7 = r4 - r5;
        bw5[j] = r6;
        bw6[j] = r7;
        r4 = (r3 << 3) + (r2 << 2) + (r1 << 1) + r0;
        bw2[j] = r4;
        bw7[j] = r0;
        bw1[j] = r3;
    }

    // MULTIPLICATION

    KaratsubaSimple(aw1, bw1, w1);
    KaratsubaSimple(aw2, bw2, w2);
    KaratsubaSimple(aw3, bw3, w3);
    KaratsubaSimple(aw4, bw4, w4);
    KaratsubaSimple(aw5, bw5, w5);
    KaratsubaSimple(aw6, bw6, w6);
    KaratsubaSimple(aw7, bw7, w7);

    // INTERPOLATION
    for (i = 0; i < N_SB_RES; ++i) {
        r0 = w1[i];
        r1 = w2[i];
        r2 = w3[i];
        r3 = w4[i];
        r4 = w5[i];
        r5 = w6[i];
        r6 = w7[i];

        r1 = r1 + r4;
        r5 = r5 - r4;
        r3 = ((r3 - r2) >> 1);
        r4 = r4 - r0;
        r4 = r4 - (r6 << 6);
        r4 = (r4 << 1) + r5;
        r2 = r2 + r3;
        r1 = r1 - (r2 << 6) - r2;
        r2 = r2 - r6;
        r2 = r2 - r0;
        r1 = r1 + 45 * r2;
        r4 = (word16)(((r4 - (r2 << 3)) * (word32)inv3) >> 3);
        r5 = r5 + r1;
        r1 = (word16)(((r1 + (r3 << 4)) * (word32)inv9) >> 1);
        r3 = -(r3 + r1);
        r5 = (word16)(((30 * r1 - r5) * (word32)inv15) >> 2);
        r2 = r2 - r4;
        r1 = r1 - r5;

        C[i]     += r6;
        C[i + 64]  += r5;
        C[i + 128] += r4;
        C[i + 192] += r3;
        C[i + 256] += r2;
        C[i + 320] += r1;
        C[i + 384] += r0;
    }
}

/**
 * Multiplies polynomials a and b and adds them to res.
 * Returns result res += a*b
 */
template<byte T_L, byte T_Et>
void SABER<T_L, T_Et>::PolyMulAcc(const word16 *a, const word16 *b, word16 *res)
{
	word16 c[2 * N] = {0};
	int i;

	ToomCook4Way(a, b, c);

	/* reduction */
	for (i = N; i < 2 * N; i++)
	{
		res[i - N] += (c[i - N] - c[i]);
	}
}

template class SABER<2,3>;
template class SABER<3,4>;
template class SABER<4,6>;


NAMESPACE_END