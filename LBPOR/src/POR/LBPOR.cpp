#include <cstdint>
#include <cmath>
#include <algorithm>
#include <memory>
#include <random>
#include "LBPOR.hpp"
#include "shake.h"
#include "gauss/gaussian.hpp"
#include "gauss/gaussian_sign.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <time.h>
using namespace Gaussian;



/*
 * Private fields
 */
struct POR::LBPOR::encapsulation
{
	const double_t sigma;
	const double_t stag;
	const double_t sigma_sign;
	const double_t ssign;
	
	const uint32_t security;
	
	Poly_t A[m_POR]; 
	Poly_t R[w_POR][k_POR-approx];
	Poly_s spk[m_sign];
	Poly_s ssk[w_POR][k_sign-approx_sign];
	
	Poly_t gs[k_POR-approx];
	Poly_s gs_sign[k_sign-approx_sign];

	uint64_t tempz[n_POR];
	uint64_t lim;
	uint64_t lim_sign;
	shake_context sc;
	
	std::unique_ptr<GaussianSampling> gaussianSampler;
	std::unique_ptr<GaussianSampling_sign> gaussianSampler_sign;
	std::unique_ptr<Gauss_t> gaussianNoise;
	std::unique_ptr<fastGauss_t> fastGaussNoise;
	std::unique_ptr<Gauss_s> gaussianNoise_sign;
	std::unique_ptr<fastGauss_s> fastGaussNoise_sign;

	void hash(Poly_t *output, unsigned char *name, uint32_t j) noexcept;
	void hash_sign(Poly_s *output, unsigned char *T,uint32_t T_size) noexcept;
	void PRF(Poly_t *output, unsigned char *name) noexcept;
	void Approx_SamplePre(Poly_t *x, Poly_t u) noexcept;
	void samplePre_sign(Poly_s *x, Poly_s u) noexcept;
	void RLWE_Approx_GenTrap(void) noexcept;
	void Sign_GenTrap(void) noexcept;
};

/*The hash function is implemented using SHKE128 and the mapping function HashToPoint */
void POR::LBPOR::encapsulation::hash(Poly_t *output, unsigned char *name, uint32_t j) noexcept
{
	name[32]=j&0x000000FF;
	name[33]=j>>8;
	
	uint64_t temp[n_POR];
	size_t ni = n_POR;
	shake_init(&sc, 256);
	shake_inject(&sc, name, 34);
	shake_flip(&sc); 
	while (ni > 0)
	{
		unsigned char buf[8];
		uint64_t tempZ;
		shake_extract(&sc, buf, 8);
		tempZ = ((uint64_t)buf[0] << 56) |((uint64_t)buf[1] << 48) | ((uint64_t)buf[2] << 40) | ((uint64_t)buf[3] << 32) | ((uint64_t)buf[4] << 24) | ((uint64_t)buf[5] << 16) | ((uint64_t)buf[6] << 8) | ((uint64_t)buf[7]);
		
		if (tempZ < lim)
		{
			temp[n_POR - ni] = tempZ % q_POR;
			ni--;
		}
	}
	output->set(temp, temp + n_POR, false);
}

void POR::LBPOR::encapsulation::hash_sign(Poly_s *output, unsigned char *T,uint32_t T_size) noexcept
{
	uint32_t temp[n_sign];
	size_t ni = n_sign;
	shake_init(&sc, 256);
	shake_inject(&sc,T, T_size);
	shake_flip(&sc);
	while (ni > 0)
	{
        unsigned char buf[4];
		uint64_t tempZ;
		shake_extract(&sc, buf, 4);
		tempZ = ((uint32_t)buf[0] << 24) | ((uint32_t)buf[1] << 16) | ((uint32_t)buf[2] << 8) | ((uint32_t)buf[3]);
		if (tempZ < lim_sign)
		{
			temp[n_sign - ni] = tempZ % q_sign;
			ni--;
		}
	}
	output->set(temp, temp + n_sign, false);
}
void POR::LBPOR::encapsulation::PRF(Poly_t *output, unsigned char *name) noexcept
{
	uint64_t temp[n_POR];
	size_t ni = n_POR;
	shake_init(&sc, 256);
	shake_inject(&sc, name, 32);
	shake_flip(&sc); 
	for (uint32_t i=0; i<d_POR;i++, ni=n_POR){
	while (ni > 0)
	{
		unsigned char buf[8];
		uint64_t tempZ;
		shake_extract(&sc, buf, 8);
		tempZ = ((uint64_t)buf[0] << 56) |((uint64_t)buf[1] << 48) | ((uint64_t)buf[2] << 40) | ((uint64_t)buf[3] << 32) | ((uint64_t)buf[4] << 24) | ((uint64_t)buf[5] << 16) | ((uint64_t)buf[6] << 8) | ((uint64_t)buf[7]);
		
		if (tempZ < lim)
		{
			temp[n_POR - ni] = tempZ % q_POR;
			ni--;
		}
	}
	output[i].set(temp, temp + n_POR, false);}
}
//The function is used to sampling preimages as block tags. 
void POR::LBPOR::encapsulation::Approx_SamplePre(Poly_t x[], Poly_t u) noexcept
{
	Poly_t perturbation[m_POR+1];
    gaussianSampler->sampleP(perturbation, A);
	
	Poly_t v;
	v = u - perturbation[m_POR];
	v.invntt_pow_invphi();
	
	uint64_t vCoefs[n_POR];
	uint32_t index = 0;
	for (auto &v_i : v.poly_obj())
		vCoefs[index++] = (uint64_t)v_i;
	Poly_t z[k_POR-approx];
	gaussianSampler->sampleGPoly(z, vCoefs);

	for (uint32_t i = 0; i < w_POR; i++)
	{
		x[i] = {perturbation[i]};
		for (uint32_t l = 0; l < k_POR-approx; l++)
		{
			x[i] = x[i] + (R[i][l]) * z[l];
		}
	}
	for (uint32_t i = w_POR; i < m_POR; i++)
		x[i] = perturbation[i] + z[i - w_POR];
}

//The function is used to sampling preimages as signatures. 
void POR::LBPOR::encapsulation::samplePre_sign(Poly_s x[], Poly_s u) noexcept
{
	Poly_s perturbation[m_sign+1];
	gaussianSampler_sign->sampleP(perturbation, spk);

	Poly_s v;
	v = u - perturbation[m_sign];
	v.invntt_pow_invphi();
	
	uint32_t vCoefs[n_sign];
	uint32_t index = 0;
	for (auto &v_i : v.poly_obj())
		vCoefs[index++] = (uint32_t)v_i;
	Poly_s z[k_sign-approx_sign];
	gaussianSampler_sign->sampleGPoly(z, vCoefs);
//std::cout<<"1sampleG1111"<<std::endl;

	for (uint32_t i = 0; i < w_POR; i++)
	{
		x[i] = {perturbation[i]};
		for (uint32_t l = 0; l < k_sign-approx_sign; l++)
		{
			x[i] = x[i] + (ssk[i][l]) * z[l];
		}
	}
	for (uint32_t i = w_POR; i < m_sign; i++)
		x[i] = perturbation[i] + z[i - w_POR];
}

//The function is used to generate public-private key pair for our POR scheme. The construction of trapdoor is based on the security of Ring-LWE problem.
void POR::LBPOR::encapsulation::RLWE_Approx_GenTrap(void) noexcept
{
    Gauss_t noise = *(gaussianNoise.get());
	for (uint32_t j = 0; j < w_POR; j++)
	{
		for (uint32_t i = 0; i < k_POR-approx; i++)
		{
			R[j][i].set(noise);
			R[j][i].ntt_pow_phi();
		}
	}
	A[0].set(1,false);
	A[0].ntt_pow_phi();
	A[1] = nfl::uniform();
	A[1].ntt_pow_phi();

	for (uint32_t j = w_POR; j < m_POR; j++)
	{
		A[j] = {gs[j - w_POR]};
		A[j] = A[j] - R[0][j - w_POR] - A[1] * R[1][j - w_POR];
	}
}

//The function is used to generate public-private key pair.
void POR::LBPOR::encapsulation::Sign_GenTrap(void) noexcept
{
	Gauss_s noise = *(gaussianNoise_sign.get());
	for (uint32_t j = 0; j < w_POR; j++)
	{
		for (uint32_t i = 0; i < k_sign-approx_sign; i++)
		{
			ssk[j][i].set(noise);
			ssk[j][i].ntt_pow_phi();
		}
	}

	spk[0].set(1,false);
	spk[0].ntt_pow_phi();
	spk[1] = nfl::uniform();
	spk[1].ntt_pow_phi();

	for (uint32_t j = w_POR; j < m_sign; j++)
	{
		spk[j] = {gs_sign[j - w_POR]};
		spk[j] = spk[j] - ssk[0][j - w_POR] - spk[1] * ssk[1][j - w_POR];
	}
}

/* Public Interface ---------------------------------------------------------------------------------- */
POR::LBPOR::LBPOR(const double_t sigma, const double_t stag, const double_t sigma_sign, const double_t ssign, const uint32_t lambda) noexcept : impl(new encapsulation{.sigma = sigma, .stag = stag, .sigma_sign=sigma_sign, .ssign = ssign, .security = lambda})
{
	impl->fastGaussNoise.reset(new fastGauss_t(sigma*O, lambda, n_POR));
	impl->gaussianNoise.reset(new Gauss_t(impl->fastGaussNoise.get()));
	impl->fastGaussNoise_sign.reset(new fastGauss_s(sigma_sign*O, lambda, n_sign));
	impl->gaussianNoise_sign.reset(new Gauss_s(impl->fastGaussNoise_sign.get()));

	NFL_POLY_COEF_TYPE gi = 1;
	gi<<=(logbase*approx);
	for (uint32_t i = 0; i < k_POR-approx; ++i)
	{
		impl->gs[i].set(gi,false);
		impl->gs[i].ntt_pow_phi();
		gi <<= logbase;
	}
	uint32_t gi_sign = 1;
	gi_sign<<=(logbase_sign*approx_sign);
	for (uint32_t i = 0; i < k_sign-approx_sign; ++i)
	{
		impl->gs_sign[i].set(gi_sign,false);
		impl->gs_sign[i].ntt_pow_phi();
		gi_sign <<= logbase_sign;
	}

	//Values of padding and related parameters in hash function
	//impl->index0 = 0;
	//impl->index1 = 1;
	impl->lim = 4 * q_POR;
	impl->lim_sign = 4 * q_sign;
}

POR::LBPOR::~LBPOR(void) noexcept
{
	//Destructor
}

/************************** Setup algorithm *************************/
void POR::LBPOR::setup(void) noexcept
{
	impl->RLWE_Approx_GenTrap();
	impl->Sign_GenTrap();
}

/************************** File-storing algorithm *************************/
void POR::LBPOR::Store(char *cname, int8_t ***M, int64_t ***X, const uint32_t l,  int32_t X_sign[][n_sign]) const noexcept
{
	Poly_t h, u, output[m_POR], temp,U[d_POR];
	Poly_s h_sign, output_sign[m_sign];
	uint64_t M_temp[n_POR], Ztemp;

	mpz_t name;
	mpz_init(name);
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	mpz_rrandomb(name, state, impl->security);
	unsigned char ucname[34];
	
	uint32_t T_length=32+2+d_POR*n_POR*8;
	unsigned char *T=new unsigned char [T_length];
	unsigned char *T2=T;
	mpz_get_str(cname, 16, name);
	for (int i = 0; i < 32; i++)
	    ucname[i]=cname[i];
	for (int i = 0; i < 32; i++){
	    *T2=cname[i];T2++;}
	*T2=l&0x000000FF;T2++;
	*T2=l>>8;T2++;

	impl->PRF(U, ucname);

	for (int i = 0; i < d_POR; i++)
	{
		for (auto &temp_i : U[i].poly_obj())
		{
            Ztemp=temp_i;
			*T2=Ztemp&0x00000000000000FF;T2++;
			Ztemp=Ztemp>>8;
			*T2=Ztemp&0x00000000000000FF;T2++;
			Ztemp=Ztemp>>8;
			*T2=Ztemp&0x00000000000000FF;T2++;
			Ztemp=Ztemp>>8;
			*T2=Ztemp&0x00000000000000FF;T2++;
			Ztemp=Ztemp>>8;
			*T2=Ztemp&0x00000000000000FF;T2++;
			Ztemp=Ztemp>>8;
			*T2=Ztemp&0x00000000000000FF;T2++;
			Ztemp=Ztemp>>8;
			*T2=Ztemp&0x00000000000000FF;T2++;
			Ztemp=Ztemp>>8;
			*T2=Ztemp&0x00000000000000FF;T2++;
		}
	}

	//perform digital signature
	impl->hash_sign(&h_sign, T,T_length);

	impl->samplePre_sign(output_sign, h_sign);

	for (uint32_t j = 0; j < m_sign; j++)
	{
		output_sign[j].invntt_pow_invphi();
		uint32_t nnum = 0;
		for (auto &temp_i : output_sign[j].poly_obj())
		{
			if (temp_i > q_sign / 2)
				X_sign[j][nnum++] = int32_t(temp_i) - q_sign;
			else
				X_sign[j][nnum++] = (int32_t)temp_i;
		}
	}

	//Generate l block tags 
	for (uint32_t i = 0; i < l; i++)
	{
		impl->hash(&h, ucname, i + 1);
		u = {h};
		
		for (uint32_t j = 0; j < d_POR; j++)
		{
			for(uint32_t mi=0;mi<n_POR;mi++)
			{ 
				if(M[i][j][mi]>=0)
				    M_temp[mi]=(uint64_t)M[i][j][mi];
				else
					M_temp[mi]=q_POR - (uint64_t)(-M[i][j][mi]);
			}
			temp.set(M_temp, M_temp + n_POR,false);
			temp.ntt_pow_phi();
			u = u + U[j] * temp;
		}
		impl->Approx_SamplePre(output, u);
		//std::cout<<"                  "    <<i<<std::endl;
		for (uint32_t j = 0; j < m_POR; j++)
		{
			output[j].invntt_pow_invphi();
			uint32_t nnum = 0;
			for (auto &temp_i : output[j].poly_obj())
			{
				if (temp_i > q_POR / 2)
					X[i][j][nnum++] = int64_t(temp_i) - q_POR;
				else
					X[i][j][nnum++] = (int64_t)temp_i;
			}
		}
	}
	delete[] T;
	mpz_clear(name);
}

/*  Initialize SamplePre algorithm  */
void POR::LBPOR::setGaussian(void) noexcept
{
	impl->gaussianSampler.reset(new GaussianSampling(impl->stag, impl->R));
	impl->gaussianSampler_sign.reset(new GaussianSampling_sign(impl->ssign, impl->ssk));
}

/****************************** Audit algorithm ********************************/
void POR::LBPOR::Audit(uint32_t I[], int16_t **V, const uint32_t c, const uint32_t l) noexcept
{
	srand(time(NULL));
	int32_t temp;
	int16_t block[l] = {0};
	for (uint32_t i = 0; i < c; i++)
	{
		do
		{
			I[i] = (rand() % l) + 1;
		} while (block[I[i] - 1] == 1);
		block[I[i] - 1] = 1;
		int16_t location[n_POR] = {0};

		for (uint32_t j = 0; j < 2 * knum; j++)
		{
			do
			{
				V[i][j] = rand() % n_POR;
			} while (location[V[i][j]] == 1);
			location[V[i][j++]] = 1;

			V[i][j] = rand() % 2;

			if (V[i][j] == 0)
				V[i][j] = -1;
		}
	}
}

/****************************** Prove algorithm ********************************/
void POR::LBPOR::Prove(Poly_t P_x[], Poly_t P_m[], uint32_t I[], int16_t **V, uint32_t c, int8_t ***M, int64_t ***X) noexcept
{
	Poly_32 sP_m[d_POR]{0};
	for (uint32_t j = 0; j < m_POR; j++)
		P_x[j] = {0};
	uint32_t index;
	uint64_t temp_X[n_POR];
	uint32_t M_temp[n_POR];
    Poly_t temp, V_temp;
	Poly_32 temps, VS_temp;
	
	for (uint32_t i = 0; i < c; i++)
	{
		index = I[i] - 1;
		uint64_t tempV[n_POR] = {0};
		uint32_t tempVS[n_POR] = {0};
		for (uint32_t j = 0; j < 2 * knum; j++)
		{
			if (V[i][j + 1] == -1)
			{
				tempV[V[i][j]] = q_POR - 1;
				tempVS[V[i][j++]] = q_small - 1;
			}
			else
			{
				tempV[V[i][j]] = 1;
				tempVS[V[i][j++]] = 1;
			}
		}
		V_temp.set(tempV, tempV + n_POR, false);
		V_temp.ntt_pow_phi();
		VS_temp.set(tempVS, tempVS + n_POR,false);
		VS_temp.ntt_pow_phi();

		for (uint32_t j = 0; j < m_POR; j++)
		{
			for (uint32_t nnum = 0; nnum < n_POR; nnum++)
				if (X[index][j][nnum] >= 0)
					temp_X[nnum] = (uint64_t)X[index][j][nnum];
				else
					temp_X[nnum] = q_POR - (uint64_t)(-X[index][j][nnum]);
			temp.set(temp_X, temp_X + n_POR, false);
			temp.ntt_pow_phi();
			temp = V_temp * temp;
			P_x[j] = P_x[j] + temp;
		}

		for (uint32_t j = 0; j < d_POR; j++)
		{
			for(uint32_t mi=0;mi<n_POR;mi++)
			{ 
				if(M[index][j][mi]>=0)
				    M_temp[mi]=(uint32_t)M[index][j][mi];
				else
					M_temp[mi]=q_small - (uint32_t)(-M[index][j][mi]);
			}
			temps.set(M_temp, M_temp + n_POR,false);
			temps.ntt_pow_phi();
			temps = VS_temp * temps;
			sP_m[j] = sP_m[j] + temps;
		}
	}

	for (uint32_t j = 0; j < m_POR; j++)
		P_x[j].invntt_pow_invphi();

    uint64_t temp_q= q_POR - (uint64_t)q_small;
	for (uint32_t j = 0; j < d_POR; j++)
	{
		sP_m[j].invntt_pow_invphi();
		index = 0;
		for (auto &h_i : sP_m[j].poly_obj())
			if (h_i > q_small / 2)
				temp_X[index++] = (uint64_t)h_i + temp_q;
			else
				temp_X[index++] = h_i;
		P_m[j].set(temp_X, temp_X + n_POR,false);
	}
}

/***************************** Verify algorithm ******************************/
int32_t POR::LBPOR::Verify(char *name, const uint32_t l,  Poly_t P_x[], Poly_t P_m[], uint32_t c, uint32_t I[], int16_t **V, int32_t X_sign[][n_sign]) noexcept
{
	Poly_t h{0}, hh{0}, error{0},U[d_POR];
	Poly_s h_sign{0}, spkx{0}, temp_s, error_sign{0};
	uint32_t temp_sign[n_sign];
    double_t beta = impl->ssign * sqrt(n_sign * m_sign);
	double_t x_norm = 0.0;
	int32_t integrity = 1;
	
	unsigned char ucname[34];
	uint64_t Ztemp;
	uint32_t T_length=32+2+d_POR*n_POR*8;
	unsigned char *T=new unsigned char [T_length];
	unsigned char *T2=T;
	
	for (int i = 0; i < 32; i++)
	    ucname[i]=name[i];
	for (int i = 0; i < 32; i++){
	    *T2=name[i];T2++;}
	*T2=l&0x000000FF;T2++;
	*T2=l>>8;T2++;

	impl->PRF(U, ucname);

	for (int i = 0; i < d_POR; i++)
	{
		for (auto &temp_i : U[i].poly_obj())
		{
            Ztemp=temp_i;
			*T2=Ztemp&0x00000000000000FF;T2++;
			Ztemp=Ztemp>>8;
			*T2=Ztemp&0x00000000000000FF;T2++;
			Ztemp=Ztemp>>8;
			*T2=Ztemp&0x00000000000000FF;T2++;
			Ztemp=Ztemp>>8;
			*T2=Ztemp&0x00000000000000FF;T2++;
			Ztemp=Ztemp>>8;
			*T2=Ztemp&0x00000000000000FF;T2++;
			Ztemp=Ztemp>>8;
			*T2=Ztemp&0x00000000000000FF;T2++;
			Ztemp=Ztemp>>8;
			*T2=Ztemp&0x00000000000000FF;T2++;
			Ztemp=Ztemp>>8;
			*T2=Ztemp&0x00000000000000FF;T2++;
		}
	}

	//Verify the validity of the signature
	impl->hash_sign(&h_sign, T, T_length);
	for (uint32_t i = 0; i < m_sign; i++)
		for (uint32_t j = 0; j < n_sign; j++)
			x_norm += double_t(X_sign[i][j] * X_sign[i][j]);
	x_norm = sqrt(x_norm);
	if (x_norm > beta)
	{
		std::cout << "The signature is invalid. The norm of the signature vector exceeds beta."<< std::endl;
		integrity = 0;
	}
	for (uint32_t i = 0; i < m_sign; i++)
	{
		for (uint32_t nnum = 0; nnum < n_sign; nnum++)
			if (X_sign[i][nnum] >= 0)
				temp_sign[nnum] = (uint32_t)X_sign[i][nnum];
			else
				temp_sign[nnum] = q_sign - (uint32_t)(-X_sign[i][nnum]);
		temp_s.set(temp_sign, temp_sign + n_sign, false); 
		
		temp_s.ntt_pow_phi();
		spkx = spkx + (impl->spk[i]) * temp_s;
	}
	
	error_sign = spkx - h_sign;
	error_sign.invntt_pow_invphi();
	beta = (double_t)sqrt(pow(base_sign,2*approx_sign)*(base_sign*base_sign+1)/(base_sign*base_sign-1))* sqrt(n_sign)*r_sign;
	x_norm = 0.0;
	for (auto &e_i : error_sign.poly_obj()){
		if (e_i >= q_sign / 2)
			e_i = (q_sign-e_i);
		x_norm +=  double_t(e_i*e_i);
	}		
	x_norm = sqrt(x_norm);
	
	if (x_norm > beta)
	{
		std::cout << "The signature is invalid. The signature does not satisfy the verification equation." << std::endl;
		integrity = 0;
	}

	//Verify the validity of the proof
	uint64_t B = 2*c * knum * p_POR;
	Poly_t temp, V_temp;
	
	for (uint32_t i = 0; i < c; i++)
	{
		impl->hash(&hh, ucname, I[i]); 
		temp = {hh};
		uint64_t tempV[n_POR] = {0};
		for (uint32_t j = 0; j < 2 * knum; j++)
		{
			if (V[i][j + 1] == -1)
				tempV[V[i][j++]] = q_POR - 1;
			else
				tempV[V[i][j++]] = 1;
		}

		V_temp.set(tempV, tempV + n_POR, false);
		V_temp.ntt_pow_phi();
		temp = V_temp * temp;
		h = h + temp;
	}
	beta = (double_t)c * knum * (impl->stag) * sqrt(n_POR * m_POR);
	x_norm = 0.0;
	for (uint32_t i = 0; i < m_POR; i++)
	{
		temp = {P_x[i]};
		for (auto &temp_i : temp.poly_obj())
		{
			if (temp_i > q_POR / 2)
				temp_i = q_POR - temp_i;

			x_norm += double_t(temp_i * temp_i);
		}
	}
	x_norm = sqrt(x_norm);
	if (x_norm > beta)
	{
		std::cout << "The proof is invalid. The norm of the homomorphic tag vector exceeds beta." << std::endl;
		integrity = 0;
	}
	error = {h};
	for (uint32_t i = 0; i < d_POR; i++)
	{
		temp = {P_m[i]};
		for (auto &temp_i : temp.poly_obj())
		{
			if (temp_i > q_POR / 2)
				temp_i = q_POR - temp_i;
			if (temp_i > B)
			{
				std::cout << "The norm of the homomorphic block vector exceeds the bound." << std::endl;
				integrity = 0;
			}
		}
		temp = {P_m[i]};
		temp.ntt_pow_phi();
		error = error + U[i] * temp;
	}
	for (uint32_t i = 0; i < m_POR; i++)
	{
		temp = {P_x[i]};
		temp.ntt_pow_phi();
		error = error - (impl->A[i]) * temp;
	}
	error.invntt_pow_invphi();
	beta = (double_t)c * knum *sqrt(pow(base_POR,2*approx)*(base_POR*base_POR+1)/(base_POR*base_POR-1)) * sqrt(n_POR)*r_POR;
	x_norm = 0.0;
	for (auto &e_i : error.poly_obj()){
		if (e_i >= q_POR / 2)
			e_i = (q_POR-e_i);
		x_norm +=  double_t(e_i*e_i);
	}		
	x_norm = sqrt(x_norm);
	
	if (x_norm > beta)
	{
		std::cout << "The proof is invalid. The error of the proof does not satisfy the verification equation." << std::endl;
		integrity = 0;
	}
	delete[] T;
	return integrity;
}
