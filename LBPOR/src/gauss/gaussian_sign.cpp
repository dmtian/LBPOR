#include<cstdint>
#include<chrono>
#include<cstring>
#include<memory>
#include<cmath>
#include<random>
#include<stack>
#include<algorithm>
#include<complex>
#include<mutex>
#include<atomic>
#include "gaussian_sign.hpp"
#include "FFT/fast_fft.hpp"
#include <string>
using namespace FFT;
/* Static Functions ---------------------------------------------------------------------------------- */
inline static void getBits(uint8_t * bits, uint32_t nl, uint16_t length) {
	uint32_t mask = nl;
	for(uint32_t i = 0; i < length; ++i) {
		bits[i] = mask & ((uint8_t)base_sign-1);
		mask >>= logbase_sign;
	}
}
/* Static Functions End ---------------------------------------------------------------------- */

/*
 * Private fields
 */
struct Gaussian::GaussianSampling_sign::encapsulation {
	const double_t s;
	uint8_t modulus[k_sign];
	double_t sigma; 
	Poly_s R[w_POR][k_approx_sign];  

	std::mt19937 * randomGenerator;
	std::unique_ptr<std::mt19937> randomGenerator_t;

    const FastFFT<n_sign> fastFft{};
	
	double_t G_l[k_sign];
	double_t G_h[k_sign];
	double_t G_d[k_sign];
	
	mutable std::stack<int32_t*> G_precomputedPerturbations;
	mutable std::mutex G_precomputed_mutex;
	void G_preCompute(int32_t *perturbation) const noexcept;
	int32_t * G_Perturb(void) const noexcept;
	inline void G_sampleD(int32_t * output, const double_t sigma, const double_t * c) const noexcept;
	void G_sampleG(int32_t * output, const uint32_t coset_i) const noexcept;
	
	double_t P_z0;
	double_t P_z1;
	double_t P_z2;
	std::complex<double_t> P_a[n_sign];
	std::complex<double_t> P_b[n_sign];
	std::complex<double_t> P_d[n_sign];
	
	mutable std::stack<Poly_s *> P_precomputedPerturbations;
	mutable std::mutex P_precomputed_mutex;
	void P_samplePz(Poly_s * output) const noexcept;
	void P_sampleFz(std::complex<double_t> * output, const std::complex<double_t> * f, const std::complex<double_t> * c, const uint32_t length) const noexcept;
	void P_sample2z(std::complex<double_t> * output_1, std::complex<double_t> * output_2, const std::complex<double_t> * f_a, const std::complex<double_t> * f_b, const std::complex<double_t> * f_d, const std::complex<double_t> * c, const uint32_t length) const noexcept;
    
	int32_t GenerateIntegerKarney(const double_t stddev, const double_t mean) const noexcept;
	bool AlgorithmP(int input) const noexcept;
	int32_t AlgorithmG(void) const noexcept;
	bool AlgorithmH(void) const noexcept;
	bool AlgorithmHDouble(void) const noexcept;
	bool AlgorithmB(int32_t D_k, double x) const noexcept ;
    bool AlgorithmBDouble(int32_t D_k, double x) const noexcept;
	void GenerateVectorKarney(uint32_t *noise) const noexcept;
};

/*Random sampling an integer on a discrete Gaussian distribution with Gaussian parameter stddev and center mean.*/
int32_t Gaussian::GaussianSampling_sign::encapsulation::GenerateIntegerKarney(const double_t stddev, const double_t mean) const noexcept{
		int32_t result;
		double_t sample_sigma=stddev*O;
		std::uniform_int_distribution<int32_t> uniform_sign(0, 1);
		std::uniform_int_distribution<int32_t> uniform_j(0, ceil(sample_sigma)-1);
		bool flagSuccess = false;
		int32_t D_k;
		while (!flagSuccess) {
			
			// STEP D1
			D_k = AlgorithmG();
			// STEP D2
			if (!AlgorithmP(D_k * (D_k - 1))) continue;
			// STEP D3
			int32_t sign = uniform_sign(*randomGenerator);
			if (sign == 0)
				sign = -1;
			// STEP D4
			double di0 = sample_sigma * D_k + sign * mean;
			int32_t i0 = std::ceil(di0);
			double x0 = (i0 - di0) / sample_sigma;
			int32_t j = uniform_j(*randomGenerator);
			double x = x0 + j / sample_sigma;
			// STEPS D5 and D6
			if (!(x < 1) || (x == 0 && sign < 0 && D_k == 0))
				continue;
			// STEP D7
			int32_t h = D_k + 1; while (h-- && AlgorithmB(D_k, x)) {};
			if (!(h < 0)) continue;
			// STEP D8
			result = sign*(i0 + j);
			flagSuccess = true;
		}
		return result;
	}

/*Random sampling vector on a discrete Gaussian distribution with Gaussian parameter P_z1 and center 0.*/
void Gaussian::GaussianSampling_sign::encapsulation::GenerateVectorKarney(uint32_t *noise) const noexcept{
	double_t sample_sigma=P_z1*O;
	std::uniform_int_distribution<int32_t> uniform_sign(0, 1);
	std::uniform_int_distribution<int32_t> uniform_j(0, ceil(sample_sigma)-1);
	bool flagSuccess = false;
	int32_t result;
	int32_t D_k;
	for(int32_t i=0;i<n_sign;i++){
	    flagSuccess = false;
		while (!flagSuccess) {
			// STEP D1
			D_k = AlgorithmG();
			// STEP D2
			if (!AlgorithmP(D_k * (D_k - 1))) continue;
			// STEP D3
			int32_t sign = uniform_sign(*randomGenerator);
			if (sign == 0)
				sign = -1;
			// STEP D4
			double di0 = sample_sigma * D_k;
			int32_t i0 = std::ceil(di0);
			double x0 = (i0 - di0) / sample_sigma;
			int32_t j = uniform_j(*randomGenerator);
			double x = x0 + j / sample_sigma;
			// STEPS D5 and D6
			if (!(x < 1) || (x == 0 && sign < 0 && D_k == 0))
				continue;
			// STEP D7
			int32_t h = D_k + 1; while (h-- && AlgorithmB(D_k, x)) {};
			if (!(h < 0)) continue;
			// STEP D8
			result = sign*(i0 + j);
			flagSuccess = true;
		}
		if(result<0)
		    noise[i]=uint32_t (q_sign+result);
	    else
		    noise[i]=(uint32_t) result;	
	}
}
		
	bool Gaussian::GaussianSampling_sign::encapsulation::AlgorithmP(int input)const noexcept{
		while (input-- && AlgorithmH()){}; return input < 0;
	}
	int32_t Gaussian::GaussianSampling_sign::encapsulation::AlgorithmG(void)const noexcept
	{
		int temp = 0; while (AlgorithmH()) ++temp; return temp;
	}
	bool Gaussian::GaussianSampling_sign::encapsulation::AlgorithmH(void)const noexcept{
		std::uniform_real_distribution<float> dist(0,1);
		float h_a, h_b;
		h_a = dist(*randomGenerator);
		if (h_a > 0.5) 
			return true;
		else if (h_a < 0.5)
		{
			for (;;) {
				h_b = dist(*randomGenerator);
				if (h_b > h_a)
					return false;
				else if (h_b < h_a)
					h_a = dist(*randomGenerator);
				else
					return AlgorithmHDouble();
				if (h_a > h_b)
					return true;
				else if (h_a == h_b)
					return AlgorithmHDouble();
			}
		}
		else
			return AlgorithmHDouble();
	}
	bool Gaussian::GaussianSampling_sign::encapsulation::AlgorithmHDouble(void)const noexcept {
		std::uniform_real_distribution<double> dist(0, 1);
		double h_a, h_b;
		h_a = dist(*randomGenerator);
		if (!(h_a < 0.5)) return true;
		for (;;) {
			h_b = dist(*randomGenerator);
			if (!(h_b<h_a))
				return false;
			else
				h_a = dist(*randomGenerator);
			if (!(h_a<h_b)) return true;
		}
	}
	bool Gaussian::GaussianSampling_sign::encapsulation::AlgorithmB(int32_t D_k, double x)const noexcept {
		std::uniform_real_distribution<float> dist(0.0, 1.0);
		float y = x;
		int32_t D_n = 0, D_m = 2 * D_k + 2;
		float z, D_r;
		float rTemp;
		for (;; ++D_n) {	
			z = dist(*randomGenerator);
			if (z > y)
				break;
			else if (z < y)
			{
				D_r = dist(*randomGenerator);
				rTemp = (2 * D_k + x) / D_m;
				if (D_r > rTemp)
					break;
				else if (D_r < rTemp)
					y = z;
				else
					return AlgorithmBDouble( D_k, x);
			}
			else
				return AlgorithmBDouble( D_k, x);
		}
    	return (D_n % 2) == 0;
	}
bool Gaussian::GaussianSampling_sign::encapsulation::AlgorithmBDouble(int32_t D_k, double x)const noexcept {
		std::uniform_real_distribution<double> dist(0.0, 1.0);
		double y = x;
		int32_t D_n = 0, D_m = 2 * D_k + 2;
		double z, D_r;
		for (;; ++D_n) {
			z = dist(*randomGenerator);
			if (!(z < y))
				break;
			D_r = dist(*randomGenerator);
			if (!(D_r < (2 * D_k + x) / D_m))
				break;
			y = z;
		}
		return (D_n % 2) == 0;
	}

/*SampleG is used to generate the preimage t with gt = u.*/
//A stack for storing the perturbation vectors of the SampleG algorithm
int32_t * Gaussian::GaussianSampling_sign::encapsulation::G_Perturb(void) const noexcept {
	std::lock_guard<std::mutex> lock(G_precomputed_mutex);
	int32_t * result = G_precomputedPerturbations.top();
	G_precomputedPerturbations.pop();
	return result;
}

void Gaussian::GaussianSampling_sign::encapsulation::G_preCompute(int32_t perturbation[]) const noexcept { 
	double_t sigmas[k_sign];
	for(uint32_t i = 0; i < k_sign; ++i) {
		sigmas[i] = sigma/G_l[i];
	}

	double_t beta = 0.0;
	int32_t z[k_sign+1];
	for(uint32_t i = 0; i < k_sign; ++i) {
		double_t c  = beta/G_l[i];
		int32_t  cz = (int32_t) floor(c);
		z[i] = cz + (int32_t)GenerateIntegerKarney(sigmas[i], c - cz);
		beta = -z[i]*G_h[i];
	}
    z[k_sign]=0;
	perturbation[0] = z[0] + ((z[0]<<1)<<logbase_sign) + z[1]<<logbase_sign;
	for(uint32_t i = 1; i < k_sign; ++i) {
		perturbation[i] = (z[i-1] + (z[i]<<1) + z[i+1])<<logbase_sign;
	}
}

inline void Gaussian::GaussianSampling_sign::encapsulation::G_sampleD(int32_t * output, const double_t sigmaa, const double_t * c) const noexcept {
	const double_t c_d = -c[k_sign - 1]/G_d[k_sign - 1];
	const int32_t c_dz = (int32_t) floor(c_d);
	output[k_sign-1] = c_dz + GenerateIntegerKarney(sigmaa/G_d[k_sign - 1], c_d - c_dz);
	double_t cs[k_sign];
	const int32_t zk_1 = output[k_sign - 1];
	for(uint32_t i = 0; i < k_sign; ++i) {
	    cs[i] =c[i]+ zk_1*G_d[i];
	}
    int32_t cs_iz;  
	for(uint32_t i = 0; i < k_sign - 1; ++i) {
		cs_iz = (int32_t) floor(-cs[i]); 
		output[i] = cs_iz + GenerateIntegerKarney(sigmaa, -cs[i] - cs_iz); 
	}
}
void Gaussian::GaussianSampling_sign::encapsulation::G_sampleG(int32_t * output, const uint32_t coset_i) const noexcept {
	uint8_t u[k_sign];
	getBits(u, coset_i, k_sign);

	int32_t pp[k_sign];
	G_preCompute(pp);

	double_t cs[k_sign];
	cs[0] = (u[0] - pp[0])/base_sign;
	for(uint32_t i = 1; i < k_sign; ++i) {
		cs[i] = (cs[i - 1] + u[i] - pp[i])/base_sign;}

	int32_t z[k_sign];
	G_sampleD(z, sigma, cs);

	output[0] = (z[0]<<logbase_sign)+modulus[0]*z[k_sign-1]+u[0];
	for(uint32_t i = 1; i < k_sign -1; ++i) {
		output[i] =  (z[i]<<logbase_sign)  - z[i-1] + modulus[i]*z[k_sign-1] + u[i];
	}
	output[k_sign-1] = modulus[k_sign-1]*z[k_sign-1] - z[k_sign-2] + u[k_sign-1];
}

/* Perturbation sampling algorithm is used to generate the perturbation vector needed in SamplePre.*/
void Gaussian::GaussianSampling_sign::encapsulation::P_samplePz(Poly_s *output) const noexcept {
	uint32_t coefs_noise[n_sign];
	for(uint32_t i = w_POR; i < m_sign; ++i) {
		GenerateVectorKarney(coefs_noise);
		output[i].set(coefs_noise,coefs_noise+n_sign,false);
		output[i].ntt_pow_phi();
	}

    Poly_s poly_c0{0};	
	Poly_s poly_c1{0};

	for(uint32_t i = 0; i < k_approx_sign; ++i) {
		poly_c0 = poly_c0 + R[0][i]*output[w_POR+i];
		poly_c1 = poly_c1 + R[1][i]*output[w_POR+i]; 
	}
	poly_c0.invntt_pow_invphi();
	poly_c1.invntt_pow_invphi();

	uint32_t index = 0;
	std::complex<double_t> * fft_c0 = new std::complex<double_t>[n_sign];
	std::complex<double_t> * fft_c1 = new std::complex<double_t>[n_sign];
	Poly_s poly_c;

	std::complex<double_t> * c = new std::complex<double_t>[n_sign];
	
	for(auto & coef : poly_c0.poly_obj()) {
		if(coef>q_sign/2){
			c[index] = std::complex<double_t>{-P_z2*(q_sign-coef), 0};}
		else{
			c[index] = std::complex<double_t>{P_z2*coef, 0};}
		index++;
	}	
	fastFft.fftForward(fft_c0, c);
	index = 0;
	for(auto & coef : poly_c1.poly_obj()) {
		if(coef>q_sign/2)
			c[index++] = std::complex<double_t>{-P_z2*(q_sign-coef), 0};
		else
			c[index++] = std::complex<double_t>{P_z2*coef, 0};
	}
	
	fastFft.fftForward(fft_c1, c);
	std::complex<double_t> * fft_c = new std::complex<double_t>[n_sign << 1];
	fastFft.ffmerge(fft_c, fft_c0, fft_c1, n_sign);
	delete[] c;
	delete[] fft_c0;
	delete[] fft_c1;
	
    std::complex<double_t> * fft_p0 = new std::complex<double_t>[n_sign];
	std::complex<double_t> * fft_p1 = new std::complex<double_t>[n_sign];
	P_sample2z(fft_p0, fft_p1, P_a, P_b, P_d, fft_c, n_sign);

	std::complex<double_t> ifft_p[n_sign];
	uint32_t coefs_p[n_sign];

	fastFft.fftInverse(ifft_p, fft_p0);
	for(uint32_t i = 0; i < n_sign; ++i) {
		if(ifft_p[i].real()<0)
		    coefs_p[i] = q_sign - (uint32_t) (-ifft_p[i].real());
		else
			coefs_p[i] = (uint32_t) ifft_p[i].real();
//std::cout<<"   "<<ifft_p[i].real()<<std::endl;
	}
	output[0].set(coefs_p, coefs_p + n_sign,false);
	output[0].ntt_pow_phi();

	fastFft.fftInverse(ifft_p, fft_p1);
	for(uint32_t i = 0; i < n_sign; ++i) {
		if(ifft_p[i].real()<0)
		    coefs_p[i] = q_sign - (uint32_t) (-ifft_p[i].real());	
		else
			coefs_p[i] = (uint32_t) ifft_p[i].real();
		
	}
	output[1].set(coefs_p, coefs_p + n_sign,false);
	output[1].ntt_pow_phi();

	delete[] fft_p0;
	delete[] fft_p1;
}

void Gaussian::GaussianSampling_sign::encapsulation::P_sample2z(std::complex<double_t> *output_1, std::complex<double_t> *output_2, const std::complex<double_t> *f_a, const std::complex<double_t> * f_b, const std::complex<double_t> * f_d, const std::complex<double_t> * c, const uint32_t length) const noexcept {

	std::complex<double_t> c0[length];
	std::complex<double_t> c1[length];
	fastFft.ffsplit(c0, c1, c, length << 1);

	P_sampleFz(output_2, f_d, c1, length);

	const std::complex<double_t> identity{1, 0};
	std::complex<double_t> *bd_1 = new std::complex<double_t>[length];
	for(uint32_t i = 0; i < length; ++i){
		bd_1[i] = f_b[i]*(identity/f_d[i]);
	}

	for(uint32_t i = 0; i < length; ++i) {
		c0[i] += bd_1[i]*(output_2[i]-c1[i]);
	}

	std::complex<double_t> tmp[length];
	for(uint32_t i = 0; i < length; ++i) {
		tmp[i] = f_a[i]- bd_1[i]*std::conj(f_b[i]);
	}
	delete[] bd_1;

	P_sampleFz(output_1, tmp, c0, length);
}

void Gaussian::GaussianSampling_sign::encapsulation::P_sampleFz(std::complex<double_t> *output, const std::complex<double_t> *f, const std::complex<double_t> *c, const uint32_t length) const noexcept {
	if (length == 1) {
		
		*output = (double_t) GenerateIntegerKarney(sqrt(fabs(f[0].real())), c[0].real());
	
	}
	else {
		const uint32_t length_2 = length >> 1;

		std::complex<double_t> f0[length_2];
		std::complex<double_t> f1[length_2];
		fastFft.ffsplit(f0, f1, f, length);

		std::complex<double_t> * q0 = new std::complex<double_t>[length_2];
		std::complex<double_t> * q1 = new std::complex<double_t>[length_2];
		P_sample2z(q0, q1, f0, f1, f0, c, length_2);

		fastFft.ffmerge(output, q0, q1, length_2);
		delete[] q0;
		delete[] q1;
	}
}

/* Private Interface End ----------------------------------------------------------------------------- */

/* Public Interface ---------------------------------------------------------------------------------- */
Gaussian::GaussianSampling_sign::GaussianSampling_sign(void) noexcept {}
Gaussian::GaussianSampling_sign::GaussianSampling_sign(const double_t s, Poly_s R[][k_approx_sign]) noexcept : \
impl(new encapsulation {.s=s}) {
    for(uint32_t i = 0; i < w_POR; ++i)
	    for(uint32_t j = 0; j<k_approx_sign ; ++j)
	        impl->R[i][j]={R[i][j]};
				
	impl->sigma = r_sign*sqrt(base_sign*base_sign+1.0)/(base_sign+1.0);
		
	getBits(impl->modulus, q_sign, k_sign);
	
	std::random_device rd;
	impl->randomGenerator_t.reset(new std::mt19937(rd()));
	impl->randomGenerator = impl->randomGenerator_t.get();
	
	impl->G_l[0] = sqrt(base_sign*(1.0 + 1.0/k_sign) + 1.0);
	impl->G_h[0] = 0;
	for(uint32_t i = 1; i < k_sign; ++i) {
		impl->G_l[i] = sqrt(base_sign*(1.0 + 1.0/(k_sign - i)));
		impl->G_h[i] = sqrt(base_sign*(1.0 - 1.0/(k_sign - i + 1)));
	}
	impl->G_d[0] = (double_t)impl->modulus[0]/base_sign;
	for(uint32_t i = 1; i < k_sign; ++i) {
		impl->G_d[i] = (impl->G_d[i - 1] + impl->modulus[i])/(double_t) base_sign;
	}

	const double_t s_2 = s*s;
	const double_t r_2 = r_sign*r_sign;
	impl->P_z0 = (s_2 - r_2)/(r_2*s_2);
	impl->P_z1 = sqrt(s_2 - r_2);
	impl->P_z2 = (-r_2)/(s_2 - r_2);

	Poly_s rr{0};
	Poly_s re{0};
	Poly_s ee{0};
	for(uint32_t i = 0; i < k_approx_sign; ++i) {
		rr = rr + R[0][i]*R[0][i];
		re = re + R[0][i]*R[1][i];
		ee = ee + R[1][i]*R[1][i];
	}
	rr.invntt_pow_invphi();
	re.invntt_pow_invphi();
	ee.invntt_pow_invphi();
	auto rrCoefs = rr.poly2mpz();
	auto reCoefs = re.poly2mpz();
	auto eeCoefs = ee.poly2mpz();
	const uint32_t coefsSize = rrCoefs.size();
	for(uint32_t i = 0; i < coefsSize; ++i){
	    if(mpz_cmp_ui(rrCoefs[i],q_sign/2)>0)
	        mpz_sub_ui(rrCoefs[i],rrCoefs[i],q_sign);
		if(mpz_cmp_ui(reCoefs[i],q_sign/2)>0)
	        mpz_sub_ui(reCoefs[i],reCoefs[i],q_sign);
		if(mpz_cmp_ui(eeCoefs[i],q_sign/2)>0)
	        mpz_sub_ui(eeCoefs[i],eeCoefs[i],q_sign);
	}
	std::complex<double_t> tmp_a[coefsSize];
	std::complex<double_t> tmp_b[coefsSize];
	std::complex<double_t> tmp_d[coefsSize];
	for(uint32_t i = 0; i < coefsSize; ++i) {
		tmp_a[i] = std::complex<double_t>(s_2 - impl->P_z0*mpz_get_ui(rrCoefs[i]), 0);
		tmp_b[i] = std::complex<double_t>(-impl->P_z0*mpz_get_ui(reCoefs[i]), 0);
		tmp_d[i] = std::complex<double_t>(s_2 - impl->P_z0*mpz_get_ui(eeCoefs[i]), 0);
	}

	impl->fastFft.fftForward(impl->P_a, tmp_a);
	impl->fastFft.fftForward(impl->P_b, tmp_b);
	impl->fastFft.fftForward(impl->P_d, tmp_d);
}

Gaussian::GaussianSampling_sign::GaussianSampling_sign(GaussianSampling_sign && other) noexcept {
	impl.reset(other.impl.get());
}
Gaussian::GaussianSampling_sign::~GaussianSampling_sign(void) noexcept {
}
Gaussian::GaussianSampling_sign & Gaussian::GaussianSampling_sign::operator= (GaussianSampling_sign && other) noexcept {
	impl.reset(other.impl.get());
	return *this;
}

/*SampleGPoly is used to generate the preimage z with Gz = u and convert z from vector to polynomial.*/
void Gaussian::GaussianSampling_sign::sampleGPoly(Poly_s * polys, const uint32_t * coset) const noexcept {
	uint32_t * samplingz = new uint32_t[n_sign*(k_approx_sign)];
	uint32_t * samplingzInd  = samplingz;
	int32_t samplingzIndex[k_sign];
	for(uint32_t i = 0; i < n_sign; ++i) {
		impl->G_sampleG(samplingzIndex, coset[i]);
		for(uint32_t j = approx_sign; j < k_sign; ++j){
			if(samplingzIndex[j]<0)
				*samplingzInd = samplingzIndex[j] + q_sign;
			else
			    *samplingzInd = samplingzIndex[j];
		samplingzInd++;}
    }
	
	uint32_t * samplingzT = new uint32_t[n_sign*(k_approx_sign)];
	const uint32_t * samplingzInd2  = samplingz;
	uint32_t * samplingzTInd = samplingzT;
	for(uint32_t i = 0; i < k_approx_sign; ++i) {
		for(uint32_t j = 0; j < n_sign; ++j) {
			*samplingzTInd++ = (uint32_t) *(samplingzInd2 + (k_approx_sign)*j);
		}
		++samplingzInd2;
	}
	delete[] samplingz;
	
	samplingzTInd = samplingzT;
	for(uint32_t i = 0; i < k_approx_sign; ++i) {
	    polys[i].set(samplingzTInd, samplingzTInd + n_sign,false);
		polys[i].ntt_pow_phi();
		samplingzTInd += n_sign;
	}
	delete[] samplingzT;
}

//Generating a perturbation vector for SamplePre
void Gaussian::GaussianSampling_sign::sampleP(Poly_s *pSamples, const Poly_s A[])  noexcept {
	    impl->P_samplePz(pSamples);

		pSamples[m_sign].set(0);
		for(uint32_t i=0;i<m_sign;i++)
		    pSamples[m_sign] = pSamples[m_sign] + A[i]*pSamples[i];
		
}
/* @Override */
Poly_s * Gaussian::GaussianSampling_sign::samplePz(void) const noexcept {
	std::lock_guard<std::mutex> lock(impl->P_precomputed_mutex);
	Poly_s * result = impl->P_precomputedPerturbations.top();
	impl->P_precomputedPerturbations.pop();
	return result;
}
