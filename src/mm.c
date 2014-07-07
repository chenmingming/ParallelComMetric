#include <stdio.h>
#include<stdlib.h>
#include<unistd.h>

#ifndef MM_C_
#define MM_C_

#if defined(__i386__)

static __inline__ unsigned long long rdtsc(void) {
	unsigned long long int x;
	__asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
	return x;
}
#elif defined(__x86_64__)

static __inline__ unsigned long long rdtsc(void)
{
	unsigned hi, lo;
	__asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
	return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

#elif defined(__powerpc__)
static __inline__ unsigned long long rdtsc(void)
{
	unsigned long long int result=0;
	unsigned long int upper, lower,tmp;
	__asm__ volatile(
			"0:                  \n"
			"\tmftbu   %0           \n"
			"\tmftb    %1           \n"
			"\tmftbu   %2           \n"
			"\tcmpw    %2,%0        \n"
			"\tbne     0b         \n"
			: "=r"(upper),"=r"(lower),"=r"(tmp)
	);
	result = upper;
	result = result<<32;
	result = result|lower;

	return(result);
}
#endif

/***********************************************************************/
/* START: MT 19937******************************************************/
/***********************************************************************/

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti = N + 1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s) {
	mt[0] = s & 0xffffffffUL;
	for (mti = 1; mti < N; mti++) {
		mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
		/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
		/* In the previous versions, MSBs of the seed affect   */
		/* only MSBs of the array mt[].                        */
		/* 2002/01/09 modified by Makoto Matsumoto             */
		mt[mti] &= 0xffffffffUL;
		/* for >32 bit machines */
	}
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length) {
	int i, j, k;
	init_genrand(19650218UL);
	i = 1;
	j = 0;
	k = (N > key_length ? N : key_length);
	for (; k; k--) {
		mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525UL))
				+ init_key[j] + j; /* non linear */
		mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
		i++;
		j++;
		if (i >= N) {
			mt[0] = mt[N - 1];
			i = 1;
		}
		if (j >= key_length)
			j = 0;
	}
	for (k = N - 1; k; k--) {
		mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941UL)) - i; /* non linear */
		mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
		i++;
		if (i >= N) {
			mt[0] = mt[N - 1];
			i = 1;
		}
	}

	mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void) {
	unsigned long y;
	static unsigned long mag01[2] = { 0x0UL, MATRIX_A };
	/* mag01[x] = x * MATRIX_A  for x=0,1 */

	if (mti >= N) { /* generate N words at one time */
		int kk;

		if (mti == N + 1) /* if init_genrand() has not been called, */
			init_genrand(5489UL); /* a default initial seed is used */

		for (kk = 0; kk < N - M; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		for (; kk < N - 1; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

		mti = 0;
	}

	y = mt[mti++];

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void) {
	return (long) (genrand_int32() >> 1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void) {
	return genrand_int32() * (1.0 / 4294967295.0);
	/* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void) {
	return genrand_int32() * (1.0 / 4294967296.0);
	/* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void) {
	return (((double) genrand_int32()) + 0.5) * (1.0 / 4294967296.0);
	/* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) {
	unsigned long a = genrand_int32() >> 5, b = genrand_int32() >> 6;
	return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
}
/* These real versions are due to Isaku Wada, 2002/01/09 added */

/***********************************************************************/
/* END: MT 19937 *******************************************************/
/***********************************************************************/

/* Standard matrix multiplication */
/* Arrays start at 0 */

double **A = NULL;
double **B = NULL;
double **C = NULL;
unsigned int matrix_size = 4096;
unsigned long rng_init_seeds[6] = { 0x0, 0x123, 0x234, 0x345, 0x456, 0x789 };
unsigned long rng_init_length = 6;
double clock_rate = 2666700000.0; // change to 1600000000.0 for Blue Gene/Q

void matrix_multiply(double **A, double **B, double **C, int size) {
	int i = 0, j = 0, k = 0;
	unsigned long long last = rdtsc();
	unsigned long long current = rdtsc();

	for (i = 0; i < size; i++) {
		if (i != 0 && i % 16 == 0) {
			last = current;
			current = rdtsc();
			printf("Completed %d in time %lf secs \n", i,
					((double) current - (double) last) / clock_rate);
		}
		for (j = 0; j < size; j++)
			for (k = 0; k < size; k++)
				C[i][j] += A[i][k] * B[k][j];
	}
}

void main_mm() {
	int i, j;

	// WHEN USING MPI DO: rng_init_seeds[0] = my_rank;
	init_by_array(rng_init_seeds, rng_init_length);

	A = (double **) calloc(matrix_size, sizeof(double*));
	for (i = 0; i < matrix_size; i++)
		A[i] = (double *) calloc(matrix_size, sizeof(double));

	B = (double **) calloc(matrix_size, sizeof(double*));
	for (i = 0; i < matrix_size; i++)
		B[i] = (double *) calloc(matrix_size, sizeof(double));

	C = (double **) calloc(matrix_size, sizeof(double*));
	for (i = 0; i < matrix_size; i++)
		C[i] = (double *) calloc(matrix_size, sizeof(double));

	for (i = 0; i < matrix_size; i++) /* Initialization */
		for (j = 0; j < matrix_size; j++) {
			A[i][j] = genrand_res53(); // (double)i*j;
			B[i][j] = genrand_res53(); //A[i][j] * A[i][j];
		}

	matrix_multiply(A, B, C, matrix_size);
}

#endif
