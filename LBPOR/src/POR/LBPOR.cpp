/*
API for the implementation of our lattice-based POR scheme. The code contains all the algorithms of our lattice-based POR scheme. The trapdoor generation algorithm (GenTrap) and the preimage sampling algorithm (SamplePre) used in the scheme are based on [MP12]. Digital signature is implemented by using a ring variant of GPV signature due to [23], which also calls GenTrap and SamplePre.

[MP12] Micciancio, D., Peikert, C.: Trapdoors for lattices: Simpler, tighter, faster,
smaller. In: EUROCRYPT 2012. LNCS, vol. 7237, pp. 700-718. Springer (2012).
[GPRR19] G¨ur, K.D., Polyakov, Y., Rohloff, K., Ryan, G.W., Sajjadpour, H., Savas, E.: Practical applications of improved gaussian sampling for trapdoor lattices. IEEE Trans. Computers 68(4), 570-584 (2019).
*/

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
	// Gauss parameters
	const double_t sigma;
	const double_t stag;
	const double_t sigma_sign;
	const double_t ssign;
	
	const uint32_t security;
	
	//Public and private keys
	Poly_t A[m]; 
	Poly_t R[w][k];
	Poly_s spk[m_sign];
	Poly_s ssk[w][k_sign];
	
	Poly_t gs[k];
	Poly_s gs_sign[k_sign];
	std::string index1, index0;
	uint64_t tempz[n];
	uint64_t lim;
	uint64_t lim_sign;
	shake_context sc;
	
	// Gaussians
	std::unique_ptr<GaussianSampling> gaussianSampler;
	std::unique_ptr<GaussianSampling_sign> gaussianSampler_sign;
	std::unique_ptr<Gauss_t> gaussianNoise;
	std::unique_ptr<fastGauss_t> fastGaussNoise;
	std::unique_ptr<Gauss_s> gaussianNoise_sign;
	std::unique_ptr<fastGauss_s> fastGaussNoise_sign;

	// Main functions
	void hash(Poly_t *output, std::string name, uint32_t j) noexcept;
	void hash_sign(Poly_s *output, std::string T) noexcept;
	void samplePre(Poly_t *x, Poly_t u) noexcept;
	void samplePre_sign(Poly_s *x, Poly_s u) noexcept;
	void RLWE_GenTrap(void) noexcept;
	void Sign_GenTrap(void) noexcept;
};

/*The hash function is implemented using SHKE128 and the mapping function HashToPoint in [Falcon18].
[Falcon18] Fouque, P.A., Hoffstein, J., Kirchner, P., Lyubashevsky, V., Pornin, T., Prest, T., Ricosset, T., Seiler, G., Whyte, W., Zhang, Z.: Falcon: Fast-fourier latticebased compact signatures over ntru (2018), submission to NIST Post-Quantum Competition*/
void POR::LBPOR::encapsulation::hash(Poly_t *output, std::string name, uint32_t j) noexcept
{
	std::string h, hj = std::to_string(j);
	uint16_t modhash = 0;
	h = name + hj + index0;
	
	uint64_t temp[n];
	size_t ni = n;
	shake_init(&sc, 256); //256:shake128, 512:shake256.
	shake_inject(&sc, (unsigned char *)h.data(), h.size());
	shake_flip(&sc); 
	while (ni > 0)
	{
		unsigned char buf[8];
		uint64_t tempZ;
		shake_extract(&sc, buf, 8); //Extract 8-byte hash strings at a time
		tempZ = ((uint64_t)buf[0] << 56) |((uint64_t)buf[1] << 48) | ((uint64_t)buf[2] << 40) | ((uint64_t)buf[3] << 32) | ((uint64_t)buf[4] << 24) | ((uint64_t)buf[5] << 16) | ((uint64_t)buf[6] << 8) | ((uint64_t)buf[7]);
		if (tempZ < lim)
		{
			temp[n - ni] = tempZ % q;
			ni--;
		}
	}
	output->set(temp, temp + n, false);
}
void POR::LBPOR::encapsulation::hash_sign(Poly_s *output, std::string T) noexcept
{
	uint32_t temp[n];
	size_t ni = n;
	shake_init(&sc, 256);
	shake_inject(&sc, (unsigned char *)T.data(), T.size());
	shake_flip(&sc);
	while (ni > 0)
	{
        unsigned char buf[8];
		uint64_t tempZ;
		shake_extract(&sc, buf, 8);
		tempZ = ((uint64_t)buf[0] << 56) |((uint64_t)buf[1] << 48) | ((uint64_t)buf[2] << 40) | ((uint64_t)buf[3] << 32) | ((uint64_t)buf[4] << 24) | ((uint64_t)buf[5] << 16) | ((uint64_t)buf[6] << 8) | ((uint64_t)buf[7]);
		if (tempZ < lim_sign)
		{
			temp[n - ni] = tempZ % q_sign;
			ni--;
		}
	}
	output->set(temp, temp + n, false);
}

//The function is used to sampling preimages as block tags. The basic process of preimage sampling algorithm see [MP12].
void POR::LBPOR::encapsulation::samplePre(Poly_t x[], Poly_t u) noexcept
{
	//Compute perturbation vector
	Poly_t perturbation[m+1];
    gaussianSampler->sampleP(perturbation, A);
	
	// compute v
	Poly_t v;
	v = u - perturbation[m];
	v.invntt_pow_invphi();
	
	// compute z
	uint64_t vCoefs[n];
	uint32_t index = 0;
	for (auto &v_i : v.poly_obj())
		vCoefs[index++] = (uint64_t)v_i;
	Poly_t z[k];
	gaussianSampler->sampleGPoly(z, vCoefs);

	// compute x
	for (uint32_t i = 0; i < w; i++)
	{
		x[i] = {perturbation[i]};
		for (uint32_t l = 0; l < k; l++)
		{
			x[i] = x[i] + (R[i][l]) * z[l];
		}
	}
	for (uint32_t i = w; i < m; i++)
		x[i] = perturbation[i] + z[i - w];
}

//The function is used to sampling preimages as signatures. The basic process of preimage sampling algorithm see [MP12].
void POR::LBPOR::encapsulation::samplePre_sign(Poly_s x[], Poly_s u) noexcept
{
	//Compute perturbation vector
	Poly_s perturbation[m_sign+1];
	gaussianSampler_sign->sampleP(perturbation, spk);
	
	// compute v
	Poly_s v;
	v = u - perturbation[m_sign];
	v.invntt_pow_invphi();
	
	// compute z
	uint64_t vCoefs[n];
	uint32_t index = 0;
	for (auto &v_i : v.poly_obj())
		vCoefs[index++] = (uint64_t)v_i;
	Poly_s z[k_sign];
	gaussianSampler_sign->sampleGPoly(z, vCoefs);

	// compute x
	for (uint32_t i = 0; i < w; i++)
	{
		x[i] = {perturbation[i]};
		for (uint32_t l = 0; l < k_sign; l++)
		{
			x[i] = x[i] + (ssk[i][l]) * z[l];
		}
	}
	for (uint32_t i = w; i < m_sign; i++)
		x[i] = perturbation[i] + z[i - w];
}

//The function is used to generate public-private key pair for our POR scheme. The construction of trapdoor is based on the security of Ring-LWE problem.
void POR::LBPOR::encapsulation::RLWE_GenTrap(void) noexcept
{
    Gauss_t noise = *(gaussianNoise.get());
	for (uint32_t j = 0; j < w; j++)
	{
		for (uint32_t i = 0; i < k; i++)
		{
			R[j][i].set(noise);
			R[j][i].ntt_pow_phi();
		}
	}

	A[0].set(1,false);
	A[0].ntt_pow_phi();
	A[1] = nfl::uniform();
	A[1].ntt_pow_phi();

	for (uint32_t j = w; j < m; j++)
	{
		A[j] = {gs[j - w]};
		A[j] = A[j] - R[0][j - w] - A[1] * R[1][j - w];
	}
}

//The function is used to generate public-private key pair for the signature scheme in [MP12]. The construction of trapdoor is based on the security of Ring-LWE problem.
void POR::LBPOR::encapsulation::Sign_GenTrap(void) noexcept
{
	Gauss_s noise = *(gaussianNoise_sign.get());
	for (uint32_t j = 0; j < w; j++)
	{
		for (uint32_t i = 0; i < k_sign; i++)
		{
			ssk[j][i].set(noise);
			ssk[j][i].ntt_pow_phi();
		}
	}

	spk[0].set(1,false);
	spk[0].ntt_pow_phi();
	spk[1] = nfl::uniform();
	spk[1].ntt_pow_phi();

	for (uint32_t j = w; j < m_sign; j++)
	{
		spk[j] = {gs_sign[j - w]};
		spk[j] = spk[j] - ssk[0][j - w] - spk[1] * ssk[1][j - w];
	}
}

/* Public Interface ---------------------------------------------------------------------------------- */
POR::LBPOR::LBPOR(const double_t sigma, const double_t stag, const double_t sigma_sign, const double_t ssign, const uint32_t lambda) noexcept : impl(new encapsulation{.sigma = sigma, .stag = stag, .sigma_sign=sigma_sign, .ssign = ssign, .security = lambda})
{
	//prebuild the Gaussian noise
	impl->fastGaussNoise.reset(new fastGauss_t(sigma*O, lambda, n));
	impl->gaussianNoise.reset(new Gauss_t(impl->fastGaussNoise.get()));
	impl->fastGaussNoise_sign.reset(new fastGauss_s(sigma_sign*O, lambda, n));
	impl->gaussianNoise_sign.reset(new Gauss_s(impl->fastGaussNoise_sign.get()));

	//Initialize the Gadget matrices
	NFL_POLY_COEF_TYPE gi = 1;
	for (uint32_t i = 0; i < k; ++i)
	{
		impl->gs[i].set(gi,false);
		impl->gs[i].ntt_pow_phi();
		gi <<= logbase;
	}
	uint64_t gi_sign = 1;
	for (uint32_t i = 0; i < k_sign; ++i)
	{
		impl->gs_sign[i].set(gi_sign,false);
		impl->gs_sign[i].ntt_pow_phi();
		gi_sign <<= logbase_sign;
	}

	//Values of padding and related parameters in hash function
	impl->index0 = std::to_string(0L);
	impl->index1 = std::to_string(1L);
	impl->lim = 4 * q;
	impl->lim_sign = 4 * q_sign;
}

POR::LBPOR::~LBPOR(void) noexcept
{
	//Destructor
}

/************************** Setup algorithm *************************/
void POR::LBPOR::setup(void) noexcept
{
	//Generate public-private key pairs.
	impl->RLWE_GenTrap();
	impl->Sign_GenTrap();
}

/************************** File-storing algorithm *************************/
void POR::LBPOR::Store(std::string &String_name, Poly_t U[], uint16_t ***M, int64_t ***X, const uint32_t l, std::string &T, int64_t X_sign[][n]) const noexcept
{
	Poly_t h, u, output[m], temp;
	Poly_s h_sign, output_sign[m_sign];
	
	//Randomly select a auxiliary polynomial vector
	for (uint32_t i = 0; i < d; i++){
		U[i] = nfl::uniform(); 
		U[i].ntt_pow_phi();
	}
	mpz_t name;
	mpz_init(name);
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	mpz_rrandomb(name, state, impl->security);
	
	//Concatenate file parameters as a string name||l||U.
	std::string TL, TU;
	char cname[32];
	mpz_get_str(cname, 16, name);
	String_name = cname;
	TL = std::to_string(l);
	T = impl->index1 + String_name + TL;
	for (int i = 0; i < d; i++)
	{
		for (auto &temp_i : U[i].poly_obj())
		{
			TU = std::to_string(temp_i);
			T += TU;
		}
	}
	
	//perform digital signature
	impl->hash_sign(&h_sign, T);
	impl->samplePre_sign(output_sign, h_sign);
	
	//Transform the data type of signature vector from R_q^m_sign to Z^(n*m_sign)
	for (uint32_t j = 0; j < m_sign; j++)
	{
		output_sign[j].invntt_pow_invphi();
		uint32_t nnum = 0;
		for (auto &temp_i : output_sign[j].poly_obj())
		{
			if (temp_i > q_sign / 2)
				X_sign[j][nnum++] = int64_t(temp_i) - q_sign;
			else
				X_sign[j][nnum++] = (int64_t)temp_i;
		}
	}
	
	//Generate l block tags 
	for (uint32_t i = 0; i < l; i++)
	{
		impl->hash(&h, String_name, i + 1);
		u = {h};
		
		for (uint32_t j = 0; j < d; j++)
		{
			temp.set(M[i][j], M[i][j] + n,false);
			temp.ntt_pow_phi();
			u = u + U[j] * (temp);
		}
		impl->samplePre(output, u);
		for (uint32_t j = 0; j < m; j++)
		{
			output[j].invntt_pow_invphi();
			uint32_t nnum = 0;
			for (auto &temp_i : output[j].poly_obj())
			{
				if (temp_i > q / 2)
					X[i][j][nnum++] = int64_t(temp_i) - q;
				else
					X[i][j][nnum++] = (int64_t)temp_i;
			}
		}
	}
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
	//Generate a challenge
	srand(time(NULL));
	int32_t temp;
	int16_t block[l] = {0};
	for (uint32_t i = 0; i < c; i++)
	{
		do
		{
			I[i] = (rand() % l) + 1;
		} while (block[I[i] - 1] == 1); //Random selection without putting back
		block[I[i] - 1] = 1;


		//Array V[][] is used to store the positions and values of non-zero coefficient in the challenge weight. For example, if V[i][j] stores a non-zero coefficient position, then V[i][j+1] stores the value of the coefficient.

			V[i][0]=0;

			V[i][1] = rand() % 3-1;



	}
}

/****************************** Prove algorithm ********************************/
void POR::LBPOR::Prove(Poly_t P_x[], Poly_t P_m[], uint32_t I[], int16_t **V, uint32_t c, uint16_t ***M, int64_t ***X) noexcept
{
	Poly_32 sP_m[d]{0};
	for (uint32_t j = 0; j < m; j++)
		P_x[j] = {0};
	uint32_t index;
	uint64_t temp_X[n];
    Poly_t temp, V_temp;
	Poly_32 temps, VS_temp;
	
	for (uint32_t i = 0; i < c; i++)
	{
		//
		index = I[i] - 1;
		uint64_t tempV[n] = {0};
		uint32_t tempVS[n] = {0};
		for (uint32_t j = 0; j < 2 * knum; j++)
		{
			//Compute the challenge weight
			if (V[i][j + 1] == -1)
			{
				tempV[V[i][j]] = q - 1;
				tempVS[V[i][j++]] = q_small - 1;
			}
			else
			{
				tempV[V[i][j]] = 1;
				tempVS[V[i][j++]] = 1;
			}
		}
		V_temp.set(tempV, tempV + n, false);
		V_temp.ntt_pow_phi();
		VS_temp.set(tempVS, tempVS + n,false);
		VS_temp.ntt_pow_phi();

		//compute homomorphic data block
		for (uint32_t j = 0; j < m; j++)
		{
			for (uint32_t nnum = 0; nnum < n; nnum++)
				if (X[index][j][nnum] >= 0)
					temp_X[nnum] = (uint64_t)X[index][j][nnum];
				else
					temp_X[nnum] = q - (uint64_t)(-X[index][j][nnum]);
			temp.set(temp_X, temp_X + n,false);
			temp.ntt_pow_phi();
			temp = V_temp * temp;
			P_x[j] = P_x[j] + temp;
		}

		//compute homomorphic tag in Poly_32 with 30-bit modulus
		for (uint32_t j = 0; j < d; j++)
		{
			temps.set(M[index][j], M[index][j] + n,false);
			temps.ntt_pow_phi();
			temps = VS_temp * temps;
			sP_m[j] = sP_m[j] + temps;
		}
	}

	for (uint32_t j = 0; j < m; j++)
		P_x[j].invntt_pow_invphi();

	//convert the data type of P_m from Poly_32 to Poly_t
    uint64_t temp_q= q - (uint64_t)q_small;
	for (uint32_t j = 0; j < d; j++)
	{
		sP_m[j].invntt_pow_invphi();
		index = 0;
		for (auto &h_i : sP_m[j].poly_obj())
			if (h_i > q_small / 2)
				temp_X[index++] = (uint64_t)h_i + temp_q;
			else
				temp_X[index++] = h_i;
		P_m[j].set(temp_X, temp_X + n,false);
	}
}

/***************************** Verify algorithm ******************************/
int32_t POR::LBPOR::Verify(std::string name, std::string T, Poly_t U[], Poly_t P_x[], Poly_t P_m[], uint32_t c, uint32_t I[], int16_t **V, int64_t X_sign[][n]) noexcept
{
	Poly_t h{0}, hh{0}, Ax{0}, Um{0};
	Poly_s h_sign{0}, spkx{0}, temp_s;
	uint64_t temp_sign[n];
    double_t beta = impl->ssign * sqrt(n * m_sign); //compute the norm bound of signature vector 
	double_t x_norm = 0.0;
	int32_t integrity = 1;
	
	//Verify the validity of the signature
	impl->hash_sign(&h_sign, T);
	for (uint32_t i = 0; i < m_sign; i++)
		for (uint32_t j = 0; j < n; j++)
			x_norm += double_t(X_sign[i][j] * X_sign[i][j]);
	x_norm = sqrt(x_norm);
	if (x_norm > beta)
	{
		std::cout << "The signature is invalid. The norm of the signature vector exceeds beta."<< std::endl;
		integrity = 0;
	}
	for (uint32_t i = 0; i < m_sign; i++)
	{
		//Convert vectors in Z^n to polynomials in R_q
		for (uint32_t nnum = 0; nnum < n; nnum++)
			if (X_sign[i][nnum] >= 0)
				temp_sign[nnum] = (uint64_t)X_sign[i][nnum];
			else
				temp_sign[nnum] = q_sign - (uint64_t)(-X_sign[i][nnum]);
		temp_s.set(temp_sign, temp_sign + n, false); 
		
		temp_s.ntt_pow_phi();
		spkx = spkx + (impl->spk[i]) * temp_s;
	}
	if (spkx - h_sign != 0)
	{
		std::cout << "The signature is invalid. The signature does not satisfy the verification equation." << std::endl;
		integrity = 0;
	}

	//Verify the validity of the proof
	uint64_t B = c * knum * p; //compute the norm bound of homomorphic data block
	Poly_t temp, V_temp;
	
	for (uint32_t i = 0; i < c; i++)//Compute the linear combination result of hashes
	{
		impl->hash(&hh, name, I[i]); 
		temp = {hh};
		uint64_t tempV[n] = {0};
		for (uint32_t j = 0; j < 2 * knum; j++)
		{
			if (V[i][j + 1] == -1)
				tempV[V[i][j++]] = q - 1;
			else
				tempV[V[i][j++]] = 1;
		}

		V_temp.set(tempV, tempV + n, false);
		V_temp.ntt_pow_phi();
		temp = V_temp * temp;
		h = h + temp;
	}
	beta = (double_t)c * knum * (impl->stag) * sqrt(n * m); //compute the norm bound of homomorphic tag
	x_norm = 0.0;
	for (uint32_t i = 0; i < m; i++)
	{
		temp = {P_x[i]};
		for (auto &temp_i : temp.poly_obj())
		{
			if (temp_i > q / 2)
				temp_i = q - temp_i;

			x_norm += double_t(temp_i * temp_i);
		}
	}
	x_norm = sqrt(x_norm);
	if (x_norm > beta)
	{
		std::cout << "The proof is invalid. The norm of the homomorphic tag vector exceeds beta." << std::endl;
		integrity = 0;
	}
	for (uint32_t i = 0; i < m; i++)
	{
		temp = {P_x[i]};
		temp.ntt_pow_phi();
		Ax = Ax + (impl->A[i]) * temp;
	}
	Um = {h};
	for (uint32_t i = 0; i < d; i++)
	{
		temp = {P_m[i]};
		for (auto &temp_i : temp.poly_obj())
		{
			if (temp_i > q / 2)
				temp_i = q - temp_i;
			if (temp_i > B)
			{
				std::cout << "The norm of the homomorphic block vector exceeds the bound." << std::endl;
				integrity = 0;
			}
		}
		temp = {P_m[i]};
		temp.ntt_pow_phi();
		Um = Um + (U[i]) * temp;
	}
	if (Ax - Um != 0)
	{
		std::cout << "The proof is invalid. The proof does not satisfy the verification equation." << std::endl;
		integrity = 0;
	}
	return integrity;
}
