#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <string>
#include <time.h>
#include "picosha2.h"
/******************************************/
#include <openssl/opensslv.h>
#include <openssl/err.h>
#include <openssl/evp.h>
#include <openssl/rsa.h>
#include <openssl/rand.h>
/******************************************/

#define MODULUS_SIZE 3072
#define lamda1 1536
#define lamda2 138 
#define s 43
#define num 10
#define n 4096

gmp_randstate_t state;
RSA *RSAPSS_KEY = NULL;
std::string index0,index1;

typedef struct
{
	mpz_t N; /* Modulus */
	mpz_t e; /* Public Exponent */
} public_key;

typedef struct
{
	mpz_t d; /* Private Exponent */
} private_key;

/************************** Key generation algorithm *************************/
void KeyGen(private_key *ku, public_key *kp)
{
	gmp_randinit_mt(state);
	mpz_t phi;
	mpz_init(phi);
	mpz_t temp;
	mpz_init(temp);
	mpz_t temp1;
	mpz_init(temp1);
	mpz_t temp2;
	mpz_init(temp2);
	mpz_t max1;
	mpz_init_set_ui(max1, 1);
	mpz_mul_2exp(max1, max1, lamda1);
	mpz_t max2;
	mpz_init_set_ui(max2, 1);
	mpz_mul_2exp(max2, max2, 2 * lamda1 + lamda2);
	do
	{
		do
		{
			mpz_rrandomb(temp, state, 2 * lamda1 + lamda2);
			mpz_nextprime(kp->e, temp);
		} while (mpz_cmp(kp->e, max2) >= 0);

		do
		{
			mpz_rrandomb(temp1, state, lamda1);
			mpz_nextprime(temp1, temp1);
		} while (mpz_cmp(temp1, max1) >= 0);

		do
		{
			mpz_rrandomb(temp2, state, lamda1);
			mpz_nextprime(temp2, temp2);
		} while (mpz_cmp(temp2, max1) >= 0);

		mpz_mul(kp->N, temp1, temp2);

		mpz_sub_ui(temp1, temp1, 1);
		mpz_sub_ui(temp2, temp2, 1);
		mpz_mul(phi, temp1, temp2);
	} while (!mpz_invert(ku->d, kp->e, phi));

	RSAPSS_KEY = RSA_generate_key(3072, 0x010001, NULL, NULL); //Generating public-private key pair of RSAPSS signature scheme.
}

/************************** File-storing algorithm *************************/
void Store(mpz_t E[], mpz_t **M, public_key kp, private_key ku, int L, mpz_t u[], std::string &T, mpz_t name, unsigned char pSignature[], unsigned char EM[])
{
	gmp_randinit_mt(state);
	mpz_t h1;
	mpz_init(h1);
	mpz_t temp;
	mpz_init(temp);
	mpz_t temp_um;
	mpz_init(temp_um);
	mpz_urandomb(name, state, 128); //Randomly select a 128-bit file identifier
	std::string String_name, TN, TU, index, H_tag, hash0, hash256;
	char cname[32], cu[768];
	mpz_get_str(cname, 16, name);
	String_name = cname;
	for (int i = 0; i < s; i++)
	{
		do
		{
			mpz_urandomm(u[i], state, kp.N); //Randomly select s auxiliary elements
		} while (mpz_cmp_ui(u[i], 0) == 0);
	}

	//Concatenate file parameters as a string name||L||u_1||...||u_s.
	TN = std::to_string(L);
	T = String_name + TN;
	for (int i = 0; i < s; i++)
	{
		mpz_get_str(cu, 16, u[i]);
		T += cu;
	}

	unsigned char pDigest[32];
	unsigned int uDigestLen = 32;
	EVP_MD_CTX *md_ctx;
	int status = 0;
	/* hash the string T */
	md_ctx = EVP_MD_CTX_new();
	EVP_DigestInit(md_ctx, EVP_sha256());
	EVP_DigestUpdate(md_ctx, (const unsigned char *)T.c_str(), T.size());
	EVP_DigestFinal(md_ctx, pDigest, &uDigestLen);
	EVP_MD_CTX_free(md_ctx);

	/* compute the PSS padded data */
	status = RSA_padding_add_PKCS1_PSS_mgf1(RSAPSS_KEY, EM, pDigest, EVP_sha256(), EVP_sha256(), -1);

	/* perform digital signature */
	status = RSA_private_encrypt(384, EM, pSignature, RSAPSS_KEY, 3);

	std::string all_hash;
	for (int i = 0; i < L; i++)
	{
		//Compute H(name||i) by using SHA-256 in counter mode.
		index = std::to_string(i + 1);
		H_tag = index0 + String_name + index;
		picosha2::hash256_hex_string(H_tag, hash0);
		all_hash = "";
		for (int j = 0; j < 12; j++)
		{
			index = std::to_string((int)j + 1);
			H_tag = index + hash0;
			picosha2::hash256_hex_string(H_tag, hash256);
			all_hash = all_hash + hash256;
		}
		mpz_set_str(h1, (char *)all_hash.data(), 16);

		//Compute block tags
		for (int j = 0; j < s; j++)
		{
			mpz_powm(temp_um, u[j], M[i][j], kp.N);
			mpz_mul(h1, h1, temp_um);
			mpz_mod(h1, h1, kp.N);
		}
		mpz_powm(E[i], h1, ku.d, kp.N);
	}

}

/************************** Audit algorithm *************************/
void Audit(int c, int L, int I[], mpz_t Q[])
{
	//Generate a challenge
	int temp[L] = {0};
	srand(time(NULL));
	gmp_randinit_mt(state);
	for (int i = 0; i < c; i++)
	{
		do
		{
			I[i] = rand() % L + 1;
		} while (temp[I[i] - 1] == 1);
		temp[I[i] - 1] = 1;
		mpz_urandomb(Q[i], state, 128);
		mpz_add_ui(Q[i], Q[i], 1);
	}
}

/************************** Prove algorithm *************************/
void Prove(mpz_t **M, mpz_t E[], int c, int I[], mpz_t Q[], mpz_t Tag, mpz_t mu[], public_key kp)
{
	mpz_t temp;
	mpz_init(temp);
	
	//compute homomorphic data block
	for (int i = 0; i < s; i++)
	{
		mpz_set_ui(mu[i], 0);
		for (int j = 0; j < c; j++)
		{
			mpz_mul(temp, Q[j], M[I[j] - 1][i]);
			mpz_add(mu[i], mu[i], temp);
		}
	}

	//compute homomorphic tag
	mpz_set_ui(Tag, 1);
	for (int i = 0; i < c; i++)
	{
		mpz_powm(temp, E[I[i] - 1], Q[i], kp.N);
		mpz_mul(Tag, Tag, temp);
		mpz_mod(Tag, Tag, kp.N);
	}
}

/************************** Verify algorithm *************************/
void Verify(int c, int I[], mpz_t Q[], mpz_t Tag, mpz_t mu[], public_key kp, mpz_t u[], mpz_t name, std::string T, unsigned char pSignature[], unsigned char EM[])
{

	std::string String_name, index, H_tag, hash256, hash0;
	mpz_t h1;
	mpz_init(h1);
	mpz_t temp;
	mpz_init(temp);
	mpz_t temp1;
	mpz_init(temp1);
	mpz_t Tag_e;
	mpz_init(Tag_e);
	char cname[768], cu[768];
	mpz_t B;
	mpz_init_set_ui(B, 1);
	mpz_mul_2exp(B, B, lamda2);
	mpz_mul(B, B, kp.N);
	mpz_mul_ui(B, B, (uint32_t)c);

	mpz_get_str(cname, 16, name);
	String_name = cname;
	
	unsigned char pDigest[32];
	unsigned int uDigestLen = 32;
	EVP_MD_CTX *md_ctx;
	unsigned char pDecrypted[384];
	int status = 0;

	/* hash the string T */
	md_ctx = EVP_MD_CTX_new();
	EVP_DigestInit(md_ctx, EVP_sha256());
	EVP_DigestUpdate(md_ctx, (const unsigned char *)T.c_str(), T.size());
	EVP_DigestFinal(md_ctx, pDigest, &uDigestLen);
	EVP_MD_CTX_free(md_ctx);

	/* Verify the signature. */
	status = RSA_public_decrypt(384, pSignature, pDecrypted, RSAPSS_KEY, RSA_NO_PADDING);
	if (status == -1)
		printf("Signature verification failed with error %s\n", ERR_error_string(ERR_get_error(), NULL));

	/* Verify the data. */
	status = RSA_verify_PKCS1_PSS_mgf1(RSAPSS_KEY, pDigest, EVP_sha256(), EVP_sha256(), pDecrypted, -1);
	if (status != 1)
		printf("RSA_verify_PKCS1_PSS failed with error %s\n", ERR_error_string(ERR_get_error(), NULL));
	
	
	//Compute the combination result of hashes
	std::string all_hash;
	mpz_set_ui(temp, 1);
	for (int i = 0; i < c; i++)
	{
		index = std::to_string(I[i]);
		H_tag = index0+ String_name + index;
		picosha2::hash256_hex_string(H_tag, hash0);
		all_hash = "";
		for (int j = 0; j < 12; j++)
		{
			index = std::to_string((int)j + 1);
			H_tag = index + hash0;
			picosha2::hash256_hex_string(H_tag, hash256);
			all_hash += hash256;
		}
		mpz_set_str(h1, (char *)all_hash.data(), 16);
		mpz_powm(temp1, h1, Q[i], kp.N);
		mpz_mul(temp, temp, temp1);
		mpz_mod(temp, temp, kp.N);
	}

	//Verify the validity of the proof
	for (int j = 0; j < s; j++)
	{
		if (mpz_cmp(u[j], B) > 0 && mpz_cmp_si(u[j], 0) < 0)
			printf("size error!\n");
		mpz_powm(temp1, u[j], mu[j], kp.N);
		mpz_mul(temp, temp, temp1);
		mpz_mod(temp, temp, kp.N);
	}
	mpz_powm(Tag_e, Tag, kp.e, kp.N);
	if (mpz_cmp(Tag_e, temp) != 0)
		printf("Proof verification failed!\n");
}

int main()
{
	//Initialization scheme parameters
    index0 = std::to_string((uint32_t)0);
	index1 = std::to_string((uint32_t)1);
	private_key ku;
	public_key kp;
	mpz_t *E = new mpz_t[n];
	mpz_t **M = new mpz_t *[n];
	for (uint32_t i = 0; i < n; i++)
		M[i] = new mpz_t[s];
	for (uint32_t i = 0; i < n; i++){
		mpz_init(E[i]);
		for (uint32_t j = 0; j < s; j++)
		{
			mpz_init(M[i][j]);
		}
	}
	mpz_t u[s];
	mpz_t mu[s];
	for (uint32_t i = 0; i < s; i++)
	{
		mpz_init(u[i]);
		mpz_init(mu[i]);
	}
	
	uint32_t over=s*3072%8192; //compute the bit length that the total length of s elements in Z_N exceeds 16KB.
	mpz_t over_mpz;
	mpz_init_set_ui(over_mpz,1);
	mpz_mul_2exp(over_mpz,over_mpz,over);
    
	std::string T;
	int I[460];
	mpz_t Q[460];
	for (uint32_t i = 0; i < 460; i++)
		mpz_init(Q[i]);
	unsigned char pSignature[384], EM[384];
	mpz_t name;
	mpz_init(name);
	mpz_t Tag;
	mpz_init(Tag);
	clock_t times[10];
	double times0=0, times1=0, times2=0, times3=0, times4=0;
	mpz_init(kp.N); /* Modulus */
	mpz_init(kp.e); /* Public Exponent */
	mpz_init(ku.d); /* Private Exponent */

	/* openssl initialization */
	ERR_load_crypto_strings();
	OpenSSL_add_all_algorithms();
	RAND_poll();

	printf("------------- Simulation experiment of the RSA-based POR scheme. --------------\n");
	
	printf("【KeyGen:......】\n");
	
	for (int k = 0; k < num; k++){
        times[0] = clock(); 
		KeyGen(&ku, &kp);
	    times[1] = clock();
	    times0 += (double)(times[1] - times[0]) / CLOCKS_PER_SEC;
	    printf("The KeyGen time is %f seconds.\n", (double)(times[1] - times[0]) /CLOCKS_PER_SEC);}
	   printf("Average KeyGen time is %f seconds.\n", times0/num);

	//Random generation of data block array. The last sector of each data block is |N|-over bit long to satisfy that each data block is 16 KB long.
	for (int i = 0; i < n; i++){
		for (int j = 0; j < s; j++)
			mpz_urandomm(M[i][j], state, kp.N);
        mpz_cdiv_q(M[i][s-1],M[i][s-1], over_mpz);
	}

	printf("【Store:......】\n");
	int L = n;
	do
	{
		//L=n;
		times1=0.0;
		for (int k = 0; k < num; k++){
			times[2] = clock();
			Store(E, M, kp, ku, L, u, T, name, pSignature, EM);
		    times[3] = clock();
		    times1+= (double)(times[3] - times[2]) / CLOCKS_PER_SEC;
			printf("The Store time is %f seconds.\n", (double)(times[3] - times[2]) /CLOCKS_PER_SEC);}
		printf("For %d data blocks containing %d sectors, average Store time is %f seconds.\n", L, s, times1/num);
	} while (L < n);

	int c;
	int list[3] = {230, 300, 460};
	for (int l = 0; l < 3; l++)
	{
		c = list[l];
		printf("\n\n......【Challenge %d blocks......】\n\n", c);
		printf("【TPA generate challenge:......】\n");
		times2=0.0;
		for (int k = 0; k < num; k++){
            times[4] = clock();
			Audit(c, L, I, Q);
		    times[5] = clock();
		    times2 += (double)(times[5] - times[4])  / CLOCKS_PER_SEC;
			printf("The generation time of challenge is %f seconds.\n", (double)(times[5] - times[4]) /CLOCKS_PER_SEC);}
		printf("Average time of challenge generation is %f seconds.\n", times2/num);

		printf("【Prove algorithm:......】\n");

		times3=0.0;
		for (int k = 0; k < num; k++){
			times[6] = clock();
			Prove(M, E, c, I, Q, Tag, mu, kp);
		    times[7] = clock();
		    times3 += (double)(times[7] - times[6]) / CLOCKS_PER_SEC;
			printf("The Prove time is %f seconds.\n", (double)(times[7] - times[6]) /CLOCKS_PER_SEC);}
		printf("Average Prove time is %f seconds.\n", times3/num);

		printf("【Verify algorithm:......】\n");
		times4=0.0;
		for (int k = 0; k < num; k++){
			times[8] = clock();
			Verify(c, I, Q, Tag, mu, kp, u, name, T, pSignature, EM);
		    times[9] = clock();
		    times4 += (double)(times[9] - times[8]) / CLOCKS_PER_SEC;
			printf("The Verify time is %f seconds.\n", (double)(times[9] - times[8]) /CLOCKS_PER_SEC);}
		printf("Average Verify time is %f seconds.\n", times4/num);
	}
	return 0;
}
