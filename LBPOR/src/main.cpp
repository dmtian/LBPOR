#include <iostream>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <complex>
#include "POR/LBPOR.hpp"
#include <fstream>
#include <random>
#include <time.h>
#include "param.hpp"
#include <unistd.h>
#include <string>
#define num 10
#define L 4096
using namespace POR;

void getTime()
{
    time_t timep;
    time (&timep);
    char tmp[64];
    strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&timep) );
		std::string T=tmp;
		std::cout<<T<<std::endl;
}
 
int main(void)
{
	double_t sigma = 254;
	double_t s1R = sigma * O * (sqrt(w_POR)*sqrt(n_POR) + sqrt(k_POR-approx)*sqrt(n_POR) + r_POR);
	double_t stag = sqrt(s1R * s1R + 1.0) * sqrt((base_POR*base_POR+1.0)+2) * r_POR;

	double_t sigma_sign = 10;
	double_t s1R_sign =  sigma_sign* O * (sqrt(w_POR)*sqrt(n_sign)  + sqrt(k_sign-approx_sign)*sqrt(n_sign)  + r_sign);
	double_t ssign = sqrt(s1R_sign * s1R_sign + 1.0) * sqrt((base_sign*base_sign+1.0)+2) * r_sign;
   getTime();
	std::cout << "Parameter: n = " << n_POR << ", q = " << q_POR << ", w = " << w_POR << std::endl;
	std::cout << "sigma: " << std::setprecision(16) << sigma << std::endl;
	std::cout << "stag: " << std::setprecision(16) << stag << std::endl;
	std::cout << "ssign: " << std::setprecision(16) << ssign << std::endl;
	std::cout << "base: "<< base_POR<<", approx: "<< approx << std::endl;
  std::cout << "block size: "<< d_POR*2<<std::endl;

	int32_t X_sign[m_sign][n_sign];
	  char name[32];

	uint32_t l, c=128;
	double_t time1, time2;
	int64_t ***X = new int64_t **[L];
	for (uint32_t i = 0; i < L; i++)
	    X[i] = new int64_t *[m_POR];
	for (uint32_t i = 0; i < L; i++)
		for (uint32_t j = 0; j < m_POR; j++)
			X[i][j] = new int64_t[n_POR];

	int8_t ***M = new int8_t **[L];
	for (uint32_t i = 0; i < L; i++)
		M[i] = new int8_t *[d_POR];
	for (uint32_t i = 0; i < L; i++)
		for (uint32_t j = 0; j < d_POR; j++)
			M[i][j] = new int8_t[n_POR];

	srand(time(NULL));
	for (uint32_t i = 0; i < L; i++)
		for (uint32_t j = 0; j < d_POR; j++)
			for (uint32_t in = 0; in < n_POR; in++){
					M[i][j][in] = rand() % p_POR;
					if(rand()%2==1)
					M[i][j][in] =-M[i][j][in]-1;}

	LBPOR LBPOR{sigma, stag, sigma_sign, ssign, 128};
	
	printf("--------------- Simulation experiment of the lattice-based POR scheme. -----------------\n");
	
	/************************** Setup algorithm *************************/
	std::cout << "Setup:" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();
	
	for (int32_t i = 0; i < num; i++)
		LBPOR.setup();

	auto end = std::chrono::high_resolution_clock::now();
	auto timing = std::chrono::duration<double, std::milli>(end - start);
	std::cout << "Average Setup Time: " << std::setprecision(16) << timing.count() / num << " ms " << std::endl;

	LBPOR.setGaussian();

	/************************** File-storing algorithm *************************/
	std::cout << "Store:" << std::endl;
	l = L>>5;
	
	do
	{
		l=l<<1;
		std::cout << "The number of data blocks: " << l << std::endl;
		time1 = 0.0, time2 = 0.0;
		for (uint32_t i = 0; i < num; ++i)
		{
			start = std::chrono::high_resolution_clock::now();
			LBPOR.Store(name, M, X, l,  X_sign);
			end = std::chrono::high_resolution_clock::now();
			timing = std::chrono::duration<double, std::milli>(end - start);
			std::cout << "Store Time of " << l << " blocks: " << std::setprecision(16) << timing.count() << " ms " << std::endl;
			time2 += timing.count();
		}
		std::cout << "Average Store Time of " << l << " blocks: " << std::setprecision(16) << time2 / num << " ms " << std::endl;
	} while (l < L);
	
	
	/************************* Audit algorithm ***********************/
	int32_t integrity;
	Poly_t P_x[m_POR], P_m[d_POR];

		std::cout << "\n \n \nThe number of challenged blocks is" << c <<".................."<< std::endl;
		uint32_t *I = new uint32_t[c];
		int16_t **V = new int16_t *[c];
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

		/************************** Verify algorithm *************************/
		time1 = 0.0;
		std::cout << "\n Verify: " << std::endl;
		for (uint32_t i = 0; i < num; i++)
		{
			start = std::chrono::high_resolution_clock::now();
			integrity = LBPOR.Verify(name,l, P_x, P_m, c, I, V, X_sign);
			end = std::chrono::high_resolution_clock::now();
			timing = std::chrono::duration<double, std::milli>(end - start);
			std::cout << "Verify Time: " << std::setprecision(16) << timing.count() << " ms " << std::endl;
			time1 += timing.count();
		}
		std::cout << "Average Verify Time of " << c << " blocks: " << std::setprecision(16) << time1/num << " ms. And integrity result is " << integrity << std::endl;

		for (uint32_t i = 0; i < c; i++)
			delete[] V[i];

		delete[] I;
		delete[] V;
	

    for (uint32_t i = 0; i < L; i++){
		for (uint32_t j = 0; j < d_POR; j++)
			delete[] M[i][j];
		for (uint32_t j = 0; j < m_POR; j++)
			delete[] X[i][j];}

	for (uint32_t i = 0; i < L; i++){
		delete[] X[i];
		delete[] M[i];}
	delete[] M;
	delete[] X;
 getTime();

	return 0;
}
