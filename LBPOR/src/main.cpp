/*
The code evaluates the average running time of each algorithm of our lattice-based POR scheme under the security level of 128-bit. These results are the average of 10 tests. The default block size of the current code is 16 KB.
*/

#include <iostream>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <complex>
#include "POR/LBPOR.hpp"
#include <fstream>
#include <random>
#include "param.hpp"
#define num 10 
#define L 64//Total number of data blocks. L*d*n*log2(p)/8/1024/1024 MB
using namespace POR;
int main(void)
{
	//Gauss parameters of the POR scheme

	double_t sigma = 230;
	double_t s1R = sigma * O * (sqrt(w)*sqrt(n) + sqrt(k)*sqrt(n) + r);
	double_t stag = sqrt(s1R * s1R + 1.0) * sqrt((base*base+1.0)+2) * r;

	//Gauss parameters of the signature scheme
	double_t sigma_sign = 230;
	double_t s1R_sign =  sigma_sign* O * (sqrt(w)*sqrt(n)  + sqrt(k_sign)*sqrt(n)  + r_sign);
	double_t ssign = sqrt(s1R_sign * s1R_sign + 1.0) * sqrt((base_sign*base_sign+1.0)+2) * r_sign;

	std::cout << "Parameter: n = " << n << ", q = " << q << ", w = " << w << std::endl;
	std::cout << "sigma: " << std::setprecision(16) << sigma << std::endl;
	std::cout << "stag: " << std::setprecision(16) << stag << std::endl;
	std::cout << "ssign: " << std::setprecision(16) << ssign << std::endl;
	std::cout << "base: "<< base<<", base of signature: "<< base_sign << std::endl;

	//Initialization of elements

	int64_t X_sign[m_sign][n];
	std::string T, name;
	Poly_t U[d];
	uint32_t l, c;
	double_t time1, time2,time3=0;

	int64_t ***X = new int64_t **[L]; //block tags

	for (uint32_t i = 0; i < L; i++)
	    X[i] = new int64_t *[m];
	for (uint32_t i = 0; i < L; i++)
		for (uint32_t j = 0; j < m; j++)
			X[i][j] = new int64_t[n];

	uint16_t ***M = new uint16_t **[L]; //data blocks
	for (uint32_t i = 0; i < L; i++)
		M[i] = new uint16_t *[d];
	for (uint32_t i = 0; i < L; i++)
		for (uint32_t j = 0; j < d; j++)
			M[i][j] = new uint16_t[n];

	//Random generation of data block array
	srand(time(NULL));
	for (uint32_t i = 0; i < L; i++)
		for (uint32_t j = 0; j < d; j++)
			for (uint32_t in = 0; in < n; in++)
				M[i][j][in] = rand() % p+1;

	LBPOR LBPOR{sigma, stag, sigma_sign, ssign, 128}; //Initialization scheme parameters
	
	printf("--------------- Simulation experiment of the lattice-based POR scheme. -----------------\n");
	
	/************************** Setup algorithm *************************/
	std::cout << "Setup:" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();
	
	for (int32_t i = 0; i < num; i++)
		LBPOR.setup();

	auto end = std::chrono::high_resolution_clock::now();
	auto timing = std::chrono::duration<double, std::milli>(end - start);
	std::cout << "Average Setup Time: " << std::setprecision(16) << timing.count() / num << " ms " << std::endl;

	LBPOR.setGaussian(); //Initialize SamplePre algorithm

	/************************** File-storing algorithm *************************/
	std::cout << "Store:" << std::endl;
	l = L;
	do
	{
		//l =L;
		std::cout << "The number of data blocks: " << l << std::endl;
		time1 = 0.0, time2 = 0.0;
		for (uint32_t i = 0; i < num; ++i)
		{
			start = std::chrono::high_resolution_clock::now();
			LBPOR.Store(name, U, M, X, l, T, X_sign);
			end = std::chrono::high_resolution_clock::now();
			timing = std::chrono::duration<double, std::milli>(end - start);
			std::cout << "Store Time of " << l << " blocks: " << std::setprecision(16) << timing.count() << " ms " << std::endl;
			time2 += timing.count();
		}
		std::cout << "Average Store Time of " << l << " blocks: " << std::setprecision(16) << time2 / num << " ms " << std::endl;
	} while (l < L);
	
	
	/************************* Audit algorithm ***********************/
	int32_t integrity; //The result returned by the Verify algorithm
	Poly_t P_x[m], P_m[d]; //proof:(P_x, P_m)
	uint32_t list[1] = {16}; //Number of challenge blocks
	for (uint32_t j = 0; j < 1; j++)
	{
		c = list[j];

		std::cout << "\n \n \nThe number of challenged blocks is" << c <<".................."<< std::endl;
		uint32_t *I = new uint32_t[c];
		int16_t **V = new int16_t *[c]; //Challenge weight
		for (uint32_t i = 0; i < c; i++)
			V[i] = new int16_t[2 * knum]; 

		time1 = 0.0;
	
		std::cout << "Generation of challenge: " << std::endl;
		for (uint32_t i = 0; i < num; i++)
		{
			start = std::chrono::high_resolution_clock::now();
			LBPOR.Audit(I, V, c, l);
			end = std::chrono::high_resolution_clock::now();
			timing = std::chrono::duration<double, std::milli>(end - start);
			std::cout << "Generation time of challenge: " << std::setprecision(16) << timing.count() << " ms " << std::endl;
			time1 += timing.count();
		}
		std::cout << "Average generation time of challenge with " << c << " challenged blocks: " << std::setprecision(16) << time1 / num << " ms " << std::endl;
		time3+=time1;
		/************************ Prove algorithm ***********************/
		time1 = 0.0;
		std::cout << "\n Prove:" << std::endl;
		for (uint32_t i = 0; i < num; i++)
		{
			start = std::chrono::high_resolution_clock::now();
			LBPOR.Prove(P_x, P_m, I, V, c, M, X);
			end = std::chrono::high_resolution_clock::now();
			timing = std::chrono::duration<double, std::milli>(end - start);
			std::cout << "Prove Time: " << std::setprecision(16) << timing.count() << " ms " << std::endl;
			time1 += timing.count();
		}
		std::cout << "Average Prove Time of " << c << " blocks: " << std::setprecision(16) << time1 / num << " ms " << std::endl;
		time3+=time1;
		/************************** Verify algorithm *************************/
		time1 = 0.0;
		std::cout << "\n Verify: " << std::endl;
		for (uint32_t i = 0; i < num; i++)
		{
			start = std::chrono::high_resolution_clock::now();
			integrity = LBPOR.Verify(name, T, U, P_x, P_m, c, I, V, X_sign);
			end = std::chrono::high_resolution_clock::now();
			timing = std::chrono::duration<double, std::milli>(end - start);
			std::cout << "Verify Time: " << std::setprecision(16) << timing.count() << " ms " << std::endl;
			time1 += timing.count();
		}
		std::cout << "Average Verify Time of " << c << " blocks: " << std::setprecision(16) << time1/num << " ms. And integrity result is " << integrity << std::endl;
		time3+=time1;
		std::cout << "************************************"<<std::endl;
		std::cout << "Average Prove+challenge+Verify Time of " << c << " blocks: " << std::setprecision(16) << (time3)/num<<"ms"<<std::endl;

		for (uint32_t i = 0; i < c; i++)
			delete[] V[i];

		delete[] I;
		delete[] V;
	}

    for (uint32_t i = 0; i < L; i++){
		for (uint32_t j = 0; j < d; j++)
			delete[] M[i][j];
		for (uint32_t j = 0; j < m; j++)
			delete[] X[i][j];}

	for (uint32_t i = 0; i < L; i++){
		delete[] X[i];
		delete[] M[i];}
	delete[] M;
	delete[] X;

	return 0;
}
