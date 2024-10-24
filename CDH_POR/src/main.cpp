#include <gmp.h>
#include </usr/local/include/pbc/pbc.h>
#include </usr/local/include/pbc/pbc_test.h>
#include <time.h>
#include <memory.h>
#include <string.h>
#include <stdio.h>
#include <string>
#include <openssl/evp.h>
#include <unistd.h>
#define cmax 460
#define num 10

/*********************** Related parameters **********************/
uint32_t s;
uint32_t n;
clock_t times[10];
pairing_t pairing;
gmp_randstate_t state;
element_t g;
mpz_t name;
element_t spk, ssk;
element_t v, alpha;
element_t sigt;
int I[cmax];
element_t Q[cmax];
element_t E;
element_t h, h1, h2;
element_t tempG1[2];
element_t tempG2;
element_t tempGT[2];
element_t tempZ;
mpz_t Q_temp;
int compress_sizeG1;
int compress_sizeG2;
unsigned char fdata[32];
unsigned char pDigest[32];
unsigned int uDigestLen = 32;
std::string TN, hash256;
std::string T, TU, H_tag, indexn, String_name, index0, index1;
int sizeZ, sizeG1, sizeG2, sizeGT, integrity;
/*****************************************************************/

void CDH_POR(void)
{
	//Initialization and Random generation of data block array
	element_t **m = new element_t *[n]; 
	for (int i = 0; i < n; i++)
		m[i] = new element_t[s];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < s; j++){
			element_init_Zr(m[i][j], pairing);
		    element_random(m[i][j]);
		}
	
	//Initialization of auxiliary elements
	element_t u[s];
	for (int i = 0; i < s; i++)
		element_init_G1(u[i], pairing);
	
	//Initialization of block tags
	element_t *e = new element_t[n];
	for (int i = 0; i < n; i++)
		element_init_G1(e[i], pairing);
	
	//Initialization of homomorphic data blocks
	element_t M[s];
	for (int i = 0; i < s; i++)
		element_init_Zr(M[i], pairing);
	
	//Initialization of strings of related elements
	unsigned char CE[compress_sizeG1];
	unsigned char Cspk[compress_sizeG2];
	unsigned char Cv[compress_sizeG2];
	unsigned char Csigt[compress_sizeG1];
	unsigned char Cu[s][compress_sizeG1];
	unsigned char Ce[n][compress_sizeG1];
    char cc[sizeZ], cu[sizeG1], cname[32];;
	
    srand(time(NULL));
	
	printf("--------------- Simulation experiment of the CDH-based POR scheme. -----------------\n");
	//----------------------------------------------------------------------------
	
	int L = 0;
	int c; //Number of indexes included in the challenge
	int list[3] = {230, 300, 460}; //Numbers of indexes included in the challenge

	do
	{
		L = n;
		double times0 = 0.0, times1 = 0.0, times2[3] = {0.0}, times3[3] = {0.0}, times4[3] = {0.0};
		
		for (int k = 0; k < num; k++)
		{
			/************************** Key generation algorithm *************************/
			printf("【KeyGen:...............】\n");
			times[0] = clock();
			
			element_random(g); //Random generation of generators
			element_random(ssk); //Random generation of private key in signature scheme
			element_pow_zn(spk, g, ssk); //Calculation of public key in signature scheme
			element_to_bytes_compressed(Cspk, spk); //The compression of elements in group G
			element_random(alpha); //Random generation of private key in POR scheme
			element_pow_zn(v, g, alpha); //Calculation of public key in POR scheme
			element_to_bytes_compressed(Cv, v);
			
			times[1] = clock();
			times0 += (double)(times[1] - times[0]) / CLOCKS_PER_SEC;
			printf("%d-th time of KeyGen algorithm is %f seconds.\n", k+1, (double)(times[1] - times[0]) / CLOCKS_PER_SEC);

			/************************** File-storing algorithm *************************/
			printf("【Store:...............】\n");
			times[2] = clock();
			
			mpz_urandomb(name, state, 128); //Randomly select a 128-bit file identifier
			for (int i = 0; i < s; i++)
			{
				element_random(u[i]); //Randomly select s auxiliary elements
				element_to_bytes_compressed(Cu[i], u[i]);
			}
			//Concatenate file parameters as a string name||L||u_1||...||u_s.
			mpz_get_str(cname, 16, name);
			T = cname;
			TN = std::to_string(L);
			String_name = cname;
			T += TN;
			for (int i = 0; i < s; i++)
			{
				element_snprint(cu, sizeG1, u[i]);
				T += cu;
			}
			
			//Hash the string T with SHA256
            T += index0;
			EVP_Digest((char *)T.data(), T.size(), pDigest, &uDigestLen, EVP_sha256(), NULL);
			element_from_hash(h1, pDigest, 32);
			
			//Generate signature using BLS signature scheme
			element_pow_zn(sigt, h1, ssk);
			element_to_bytes_compressed(Csigt, sigt);
			
			//Generate L block tags
			for (int i = 0; i < L; i++)
			{
				//Compute Hash(name||i)
				indexn = std::to_string((int)i + 1);
				H_tag = String_name + indexn + index1;
				EVP_Digest((char *)H_tag.data(), H_tag.size(), pDigest, &uDigestLen, EVP_sha256(), NULL);
				element_from_hash(h1, pDigest, 32);
				
				//Compute i-th block tag e[i].
				for (int j = 0; j < s; j++)
				{
					element_pow_zn(tempG1[0], u[j], m[i][j]);
					element_mul(h1, h1, tempG1[0]);
				}
				element_pow_zn(e[i], h1, alpha);
				element_to_bytes_compressed(Ce[i], e[i]);
			}
			times[3] = clock();
			times1 += (double)(times[3] - times[2]) / CLOCKS_PER_SEC;

			printf("For %d data blocks containing %d sectors, the %d-th Store time is %f seconds.\n", L,s, k+1, (double)(times[3] - times[2]) / CLOCKS_PER_SEC);

			/****************************** Audit algorithm ********************************/

			for (int l = 0; l < 3; l++)
			{
				c = list[l];
				printf("\n【Challenge %d blocks...............】\n", c);
				printf("【TPA generate challenge:......】\n");
				times[4] = clock();
				
				//Randomly select c tuples {I[i],Q[i]}.
				int temp[L] = {0};
				for (int i = 0; i < c; i++)
				{
					do
					{
						I[i] = rand() % L + 1;
					} while (temp[I[i] - 1] == 1);//Random selection without putting back
					temp[I[i] - 1] = 1;
					mpz_urandomb(Q_temp, state, 128);
					element_set_mpz(Q[i], Q_temp);
				}

				times[5] = clock();
				times2[l] += (double)(times[5] - times[4]) / CLOCKS_PER_SEC;
				printf("%d-th time of generating a challenge is %f seconds.\n", k+1, (double)(times[5] - times[4]) / CLOCKS_PER_SEC);

				/****************************** Prove algorithm ********************************/
				printf("【Prove algorithm:..........】\n");

				times[6] = clock();

				//compute homomorphic data block
				for (int i = 0; i < s; i++)
				{
					element_set0(M[i]);
					for (int j = 0; j < c; j++)
					{
						element_mul(tempZ, Q[j], m[I[j] - 1][i]);
						element_add(M[i], M[i], tempZ);
					}
				}
				
				//compute homomorphic tag
				element_set1(E);
				for (int i = 0; i < c; i++)
				{
					element_from_bytes_compressed(e[I[i] - 1], Ce[I[i] - 1]);
					element_pow_zn(tempG1[0], e[I[i] - 1], Q[i]);
					element_mul(E, E, tempG1[0]);
				}
				element_to_bytes_compressed(CE, E);

				times[7] = clock();
				times3[l] += (double)(times[7] - times[6]) / CLOCKS_PER_SEC;
				printf("%d-th Prove time is %f seconds.\n", k+1, (double)(times[7] - times[6]) / CLOCKS_PER_SEC);
			
                /***************************** Verify algorithm ******************************/
				printf("【Verify algorithm:.........】\n");
				times[8] = clock();
				integrity = 1; //If integrity = 1, validation passed; otherwise, validation failed.
				element_from_bytes_compressed(spk, Cspk);
				element_from_bytes_compressed(sigt, Csigt);
				element_from_bytes_compressed(E, CE);
				element_from_bytes_compressed(v, Cv);
				
				//Verify the validity of the signature
				EVP_Digest((char *)T.data(), T.size(), pDigest, &uDigestLen, EVP_sha256(), NULL);
				element_from_hash(h1, pDigest, 32);
				pairing_apply(tempGT[0], sigt, g, pairing);
				pairing_apply(tempGT[1], h1, spk, pairing);
				if (element_cmp(tempGT[0], tempGT[1]))
				{
					printf("The signature verification failed. \n");
					integrity = 0;
				}
				
				//Compute the combination result of hashes
				element_set1(tempG1[0]);
				for (int i = 0; i < c; i++)
				{
					indexn = std::to_string(I[i]);
					H_tag = String_name + indexn + index1;
					EVP_Digest((char *)H_tag.data(), H_tag.size(), pDigest, &uDigestLen, EVP_sha256(), NULL);
					element_from_hash(h1, pDigest, 32);
					element_pow_zn(tempG1[1], h1, Q[i]);
					element_mul(tempG1[0], tempG1[0], tempG1[1]);
				}
				
				//Verify the validity of the proof
				for (int i = 0; i < s; i++)
				{
					element_from_bytes_compressed(u[i], Cu[i]);
					element_pow_zn(tempG1[1], u[i], M[i]);
					element_mul(tempG1[0], tempG1[0], tempG1[1]);
				}
				pairing_apply(tempGT[0], E, g, pairing);
				pairing_apply(tempGT[1], tempG1[0], v, pairing);
				if (element_cmp(tempGT[0], tempGT[1]))
				{
					integrity = 0;
					printf("The proof verification failed\n");
				}

				times[9] = clock();
				times4[l] += (double)(times[9] - times[8]) / CLOCKS_PER_SEC;
				printf("%d-th Verify time is %f seconds, integrity = %u\n", k+1,(double)(times[9] - times[8]) / CLOCKS_PER_SEC, integrity);
			}
		}
		printf("Average time of %d key generations is %f seconds.\n", num, times0 / num);
		printf("For %d data blocks containing %d sectors, average time of %d times Store is %f seconds.\n", L, s, num, times1 / num);
		printf("For c=230, average time of %d times challenge generation is %f seconds.\n", num, times2[0] / num);
		printf("For c=230, average time of %d times Prove is %f seconds.\n", num, times3[0] / num);
		printf("For c=230, average time of %d times Verify is %f seconds.\n", num, times4[0] / num);
		printf("For c=300, average time of %d times challenge generation is %f seconds.\n", num, times2[1] / num);
		printf("For c=300, average time of %d times Prove is %f seconds.\n", num, times3[1] / num);
		printf("For c=300, average time of %d times Verify is %f seconds.\n", num, times4[1] / num);
		printf("For c=460, average time of %d times challenge generation is %f seconds.\n", num, times2[2] / num);
		printf("For c=460, average time of %d times Prove is %f seconds.\n", num, times3[2] / num);
		printf("For c=460, average time of %d times Verify is %f seconds.\n", num, times4[2] / num);
	} while (L < n);

	for (int i = 0; i < n; i++)
		delete[] * (m + i);
	delete[] m;
	delete[] e;
}

int main(int argc, char **argv)
{
	/************** Initialize main parameters *******************/
	//Parameter generation of BN curve
	/*pbc_param_t param;
	pbc_param_init_f_gen(param,256);
	FILE *curve;
	curve=fopen("F_128.param","wb");
	pbc_param_out_str(curve, param);
	fclose(curve);*/

	pbc_demo_pairing_init(pairing, argc, argv); //Initializes bilinear pairing
	mpz_init(name);
	element_init_G2(g, pairing);
	element_init_Zr(ssk, pairing);
	element_init_G2(spk, pairing);
	element_init_Zr(alpha, pairing); 
	element_init_G2(v, pairing);	 
	element_init_G1(sigt, pairing);
	element_init_Zr(h, pairing);
	element_init_Zr(tempZ, pairing);
	element_init_G2(tempG2, pairing);
	for (int i = 0; i < 2; i++)
		element_init_G1(tempG1[i], pairing);
	for (int i = 0; i < 2; i++)
		element_init_GT(tempGT[i], pairing);
	element_init_G1(h1, pairing);
	element_init_Zr(h2, pairing);
	element_init_G1(E, pairing);
	for (int i = 0; i < cmax; i++)
		element_init_Zr(Q[i], pairing);
	mpz_init(Q_temp);
	/****************************************************************/
	
	//Compute the byte length of each type of element
	element_random(tempZ);
	element_random(tempG1[0]);
	element_random(tempG2);
	compress_sizeG1 = pairing_length_in_bytes_compressed_G1(pairing);
	compress_sizeG2 = pairing_length_in_bytes_compressed_G2(pairing);
	sizeZ = element_length_in_bytes(tempZ);
	sizeG1 = element_length_in_bytes(tempG1[0]);
	sizeG2 = element_length_in_bytes(tempG2);
	sizeGT = element_length_in_bytes(tempGT[0]);
	printf("sizeZ=%d,  sizeG1=%d,  sizeG2=%d, sizeGT=%d,  size_compressedG1=%d,  size_compressedG2=%d. \n", sizeZ, sizeG1, sizeG2, sizeGT, compress_sizeG1, compress_sizeG2);
	
	//Used to add different suffixes to different hash functions
	index0 = std::to_string(0L);
	index1 = std::to_string(1L);
	
	gmp_randinit_mt(state); //Initialize a seed of the random function
	
	s = 512, n = 4096;
	CDH_POR();
	sleep(600);

	s = 128, n = 16384;
	CDH_POR();
	sleep(600);
	
	s = 256, n = 8192;
	CDH_POR();
	sleep(600);
	
	s = 1024, n = 2048;
	CDH_POR();
	sleep(600);

	s = 2048, n = 1024;
	CDH_POR();
	sleep(600);
	return 0;
}
