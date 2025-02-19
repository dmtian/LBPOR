/*
The code mainly implements SampleG and SampleP included in the algorithm SamplePre. The sampling algorithm, SampleG, is used to generate the preimage z with Gz = u and described in Figure 2 of [GM18]. The perturbation sampling algorithm, SampleP, is used to generate the perturbation vectors for SamplePre and described in Figure 4 of [GM18].

[GM18] Genise, N., Micciancio, D.: Faster gaussian sampling for trapdoor lattices with arbitrary modulus. In: EUROCRYPT 2018. LNCS, vol. 10820, pp. 174-203.
*/

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
#include "gaussian.hpp"
#include "FFT/fast_fft.hpp"
#include <string>
using namespace FFT;
/* Static Functions ---------------------------------------------------------------------------------- */
//Decomposing the integer nl into bits[0]*base^0+bits[1]*base^1+...+bits[length-1]*base^(length-1)
inline static void getBits(uint8_t * bits, uint64_t nl, uint32_t length) {
	uint64_t mask = nl;
	for(uint32_t i = 0; i < length; ++i) {
		bits[i] = mask & ((uint8_t)base-1);
		mask >>= logbase;
	}
}
/* Static Functions End --------------------------------------------------------------------- */
/*
 * Private fields
 */
struct Gaussian::GaussianSampling::encapsulation {
	const double_t s;
	uint8_t modulus[k];  // The base-ary decomposition sequence of module q
	double_t sigma;    
	Poly_t R[w][k];  

	//Random number generator
	std::mt19937 * randomGenerator;
	std::unique_ptr<std::mt19937> randomGenerator_t;
	
    //fast Fourier transform
	const FastFFT<n> fastFft{};
	
	// Variables associated to sampleG
	double_t G_l[k];
	double_t G_h[k];
	double_t G_d[k];
	
	// Functions associated to sampleG
	mutable std::stack<int32_t*> G_precomputedPerturbations;
	mutable std::mutex G_precomputed_mutex;
	void G_preCompute(int32_t *perturbation) noexcept;
	int32_t * G_Perturb(void)  noexcept;
	inline void G_sampleD(int64_t * output, const double_t sigma, const double_t * c) noexcept;
	void G_sampleG(int64_t * output, const uint64_t coset_i) noexcept;
	
	// Variables associated to sampleP
	double_t P_z0;                // z
	double_t P_z1;                // sqrt(s^2 - alpha^2)
	double_t P_z2;                // (-alpha^2)/(s^2 - alpha^2)
	std::complex<double_t> P_a[n];     // a
	std::complex<double_t> P_b[n];     // b
	std::complex<double_t> P_d[n];     // d

	// Functions associated to sampleP
	mutable std::stack<Poly_t *> P_precomputedPerturbations;
	mutable std::mutex P_precomputed_mutex;
	void P_samplePz(Poly_t * output)  noexcept;
	void P_sampleFz(std::complex<double_t> * output, const std::complex<double_t> * f, const std::complex<double_t> * c, const uint32_t length) noexcept;
	void P_sample2z(std::complex<double_t> * output_1, std::complex<double_t> * output_2, const std::complex<double_t> * f_a, const std::complex<double_t> * f_b, const std::complex<double_t> * f_d, const std::complex<double_t> * c, const uint32_t length) noexcept;
    
    //Random sampling functions of discrete Gaussian distribution
	int64_t GenerateIntegerKarney(const double_t stddev, const double_t mean) const noexcept;
	void GenerateVectorKarney(uint64_t *noise) const noexcept;
	bool AlgorithmP(int input) const noexcept;
	int32_t AlgorithmG(void) const noexcept;
	bool AlgorithmH(void) const noexcept;
	bool AlgorithmHDouble(void) const noexcept;
	bool AlgorithmB(int32_t D_k, double x) const noexcept ;
    bool AlgorithmBDouble(int32_t D_k, double x) const noexcept;
};

/*Random sampling an integer on a discrete Gaussian distribution with Gaussian parameter stddev and center mean. This function uses the algorithm in [kar16], including subfunctions AlgorithmP(), AlgorithmG(), AlgorithmH(), AlgorithmHDouble() and AlgorithmB().

[Kar16] Karney, C.F.F.: Sampling exactly from the normal distribution. ACM Trans. Math.
Softw. 42(1), 3:1-3:14 (2016).*/
int64_t Gaussian::GaussianSampling::encapsulation::GenerateIntegerKarney(const double_t stddev, const double_t mean) const noexcept{
		int64_t result;
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
	
/*Random sampling vector on a discrete Gaussian distribution with Gaussian parameter P_z1 and center 0. This function uses the algorithm in [kar16].*/
void Gaussian::GaussianSampling::encapsulation::GenerateVectorKarney(uint64_t *noise) const noexcept{
	double_t sample_sigma=P_z1*O;
	std::uniform_int_distribution<int32_t> uniform_sign(0, 1);
	std::uniform_int_distribution<int32_t> uniform_j(0, ceil(sample_sigma)-1);
	bool flagSuccess = false;
	int64_t result;
	int32_t D_k;
	for(int32_t i=0;i<n;i++){
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
			double x0 = (i0 - di0) /sample_sigma;
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
		    noise[i]=uint64_t (q+result);
	    else
		    noise[i]=(uint64_t) result;	
	}
}
		
bool Gaussian::GaussianSampling::encapsulation::AlgorithmP(int input)const noexcept{
	while (input-- && AlgorithmH()){}; return input < 0;
}

int32_t Gaussian::GaussianSampling::encapsulation::AlgorithmG(void)const noexcept
{
	int temp = 0; while (AlgorithmH()) ++temp; return temp;
}

bool Gaussian::GaussianSampling::encapsulation::AlgorithmH(void)const noexcept{
	std::uniform_real_distribution<float> dist(0,1);
	float h_a, h_b;
	h_a = dist(*randomGenerator);

	// less than the half
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
			else //numbers are equal - need higher precision
				return AlgorithmHDouble();
			if (h_a > h_b)
				return true;
			else if (h_a == h_b) //numbers are equal - need higher precision
				return AlgorithmHDouble();
		}
	}
	else //numbers are equal - need higher precision
		return AlgorithmHDouble();
}
	
bool Gaussian::GaussianSampling::encapsulation::AlgorithmHDouble(void)const noexcept {
	std::uniform_real_distribution<double> dist(0, 1);
	double h_a, h_b;
	h_a = dist(*randomGenerator);
	// less than the half
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
	
bool Gaussian::GaussianSampling::encapsulation::AlgorithmB(int32_t D_k, double x)const noexcept {
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
			else // D_r == Temp - need double precision
				return AlgorithmBDouble( D_k, x);
		}
		else // z == x - need double precision
			return AlgorithmBDouble( D_k, x);
	}
    return (D_n % 2) == 0;
}

bool Gaussian::GaussianSampling::encapsulation::AlgorithmBDouble(int32_t D_k, double x)const noexcept {
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

//A stack for storing the perturbation vectors of the SampleG algorithm
int32_t * Gaussian::GaussianSampling::encapsulation::G_Perturb(void)  noexcept {
	std::lock_guard<std::mutex> lock(G_precomputed_mutex);
	int32_t * result = G_precomputedPerturbations.top();
	G_precomputedPerturbations.pop();
	return result;
}

/*SampleG is used to generate the preimage t with gt = u. The algorithm SampleG includes G_sampleG, G_sampleD and G_Perturb(G_preCompute), corresponding to sampleG, sampleD and perturb in [GM18], respectively.
*/
void Gaussian::GaussianSampling::encapsulation::G_preCompute(int32_t *perturbation)  noexcept {
	// perturb(sigma) in SampleG
	double_t sigmas[k];
	for(uint32_t i = 0; i < k; ++i) {
		sigmas[i] = sigma/G_l[i];
	}

	double_t beta = 0.0;
	int32_t z[k+1];
	for(uint32_t i = 0; i < k; ++i) {
		double_t c  = beta/G_l[i];
		int32_t cz = (int32_t) floor(c);
		z[i] = cz + (int32_t)GenerateIntegerKarney(sigmas[i], c - cz);
		beta = -z[i]*G_h[i];
	}
    z[k]=0;
	perturbation[0] = z[0] + ((z[0]<<1)<<logbase) + z[1]<<logbase;
	for(uint32_t i = 1; i < k; ++i) {
		perturbation[i] = (z[i-1] + (z[i]<<1) + z[i+1])<<logbase;
	}
}

inline void Gaussian::GaussianSampling::encapsulation::G_sampleD(int64_t * output, const double_t sigmaa, const double_t * c)noexcept {
	// SampleD(sigma,c) in SampleG
	const double_t c_d = -c[k - 1]/G_d[k - 1];
	const int64_t c_dz = (int64_t) floor(c_d);
	output[k-1] = c_dz + GenerateIntegerKarney(sigmaa/G_d[k - 1], c_d - c_dz);
	double_t cs[k];
	const int64_t zk_1 = output[k - 1];
	for(uint32_t i = 0; i < k; ++i) {
	    cs[i] =c[i]+ zk_1*G_d[i];
	}
    int64_t cs_iz;  
	for(uint32_t i = 0; i < k - 1; ++i) {
		cs_iz = (int64_t) floor(-cs[i]); 
		output[i] = cs_iz + GenerateIntegerKarney(sigmaa, -cs[i] - cs_iz);  
	}
}
void Gaussian::GaussianSampling::encapsulation::G_sampleG(int64_t * output, const uint64_t coset_i)  noexcept {
	// get bits of coset_i
	uint8_t u[k];
	getBits(u, coset_i, k);
	// get perturbations
	int32_t pp[k];
	G_preCompute(pp);

	// compute cs
	double_t cs[k];
	cs[0] = (u[0] - pp[0])/base;
	for(uint32_t i = 1; i < k; ++i) {
		cs[i] = (cs[i - 1] + u[i] - pp[i])/base;}
	// compute z
	int64_t z[k];
	G_sampleD(z, sigma, cs);
	// compute t
	output[0] = (z[0]<<logbase)+modulus[0]*z[k-1]+u[0];
	for(uint32_t i = 1; i < k -1; ++i) {
		output[i] = (z[i]<<logbase) - z[i-1] + modulus[i]*z[k-1] + u[i];
	}
	output[k-1] = modulus[k-1]*z[k-1] - z[k-2] + u[k-1];
}

/* Perturbation sampling algorithm is used to generate the perturbation vector needed in SamplePre. The algorithm SampleP includes P_samplePz, P_sample2z and P_sampleFz, corresponding to samplePz, sample2z and sampleFz in [GM18], respectively. The algorithm needs to call FFT. 
*/
void Gaussian::GaussianSampling::encapsulation::P_samplePz(Poly_t *output)  noexcept {
	uint64_t coefs_noise[n];
	for(uint32_t i = w; i < m; ++i) {
		GenerateVectorKarney(coefs_noise);
		output[i].set(coefs_noise,coefs_noise+n,false);
		output[i].ntt_pow_phi();
	}
	
    Poly_t poly_c0{0};	
	Poly_t poly_c1{0};

	for(uint32_t i = 0; i < k; ++i) {
		poly_c0 = poly_c0 + R[0][i]*output[w+i];
		poly_c1 = poly_c1 + R[1][i]*output[w+i]; 
	}
	poly_c0.invntt_pow_invphi();
	poly_c1.invntt_pow_invphi();

	uint32_t index = 0;
	std::complex<double_t> * fft_c0 = new std::complex<double_t>[n];
	std::complex<double_t> * fft_c1 = new std::complex<double_t>[n];
	Poly_t poly_c;
	std::complex<double_t> * c = new std::complex<double_t>[n];
	
	for(auto & coef : poly_c0.poly_obj()) {
		if(coef>q/2){
			c[index] = std::complex<double_t>{-P_z2*(q-coef), 0};}
		else{
			c[index] = std::complex<double_t>{P_z2*coef, 0};}
		index++;
	}	
	fastFft.fftForward(fft_c0, c);
	index = 0;
	for(auto & coef : poly_c1.poly_obj()) {
		if(coef>q/2)
			c[index++] = std::complex<double_t>{-P_z2*(q-coef), 0};
		else
			c[index++] = std::complex<double_t>{P_z2*coef, 0};
	}
	fastFft.fftForward(fft_c1, c);
	std::complex<double_t> * fft_c = new std::complex<double_t>[n << 1];
	fastFft.ffmerge(fft_c, fft_c0, fft_c1, n);
	delete[] c;
	delete[] fft_c0;
	delete[] fft_c1;
	
	// compute p0 and p1
    std::complex<double_t> * fft_p0 = new std::complex<double_t>[n];
	std::complex<double_t> * fft_p1 = new std::complex<double_t>[n];
	P_sample2z(fft_p0, fft_p1, P_a, P_b, P_d, fft_c, n);
 
	// fft inverse
	std::complex<double_t> ifft_p[n];
	NFL_POLY_COEF_TYPE coefs_p[n];
	
	fastFft.fftInverse(ifft_p, fft_p0);
	for(uint32_t i = 0; i < n; ++i) {
		if(ifft_p[i].real()<0)
		    coefs_p[i] = q - (NFL_POLY_COEF_TYPE) (-ifft_p[i].real());
		else
			coefs_p[i] = (NFL_POLY_COEF_TYPE) ifft_p[i].real();
	}
	output[0].set(coefs_p, coefs_p + n,false);
	output[0].ntt_pow_phi();
		
	fastFft.fftInverse(ifft_p, fft_p1);
	for(uint32_t i = 0; i < n; ++i) {
		if(ifft_p[i].real()<0)
		    coefs_p[i] = q - (NFL_POLY_COEF_TYPE) (-ifft_p[i].real());	
		else
			coefs_p[i] = (NFL_POLY_COEF_TYPE) ifft_p[i].real();
	}
	output[1].set(coefs_p, coefs_p + n,false);
	output[1].ntt_pow_phi();
	delete[] fft_p0;
	delete[] fft_p1;
}

void Gaussian::GaussianSampling::encapsulation::P_sample2z(std::complex<double_t> *output_1, std::complex<double_t> *output_2, const std::complex<double_t> *f_a, const std::complex<double_t> * f_b, const std::complex<double_t> * f_d, const std::complex<double_t> * c, const uint32_t length)  noexcept {

	//let c(x) = c0(x^2) + x.c1(x^2)
	std::complex<double_t> c0[length];
	std::complex<double_t> c1[length];
	fastFft.ffsplit(c0, c1, c, length << 1);

	//q1 <-- sampleFz(d, c1)
	P_sampleFz(output_2, f_d, c1, length);

	//compute bd^{-1}
	const std::complex<double_t> identity{1, 0};
	std::complex<double_t> *bd_1 = new std::complex<double_t>[length];
	for(uint32_t i = 0; i < length; ++i){
		bd_1[i] = f_b[i]*(identity/f_d[i]);
	}

	//c0 := c0 + bd^{-1}(q1 - c1)
	for(uint32_t i = 0; i < length; ++i) {
		c0[i] += bd_1[i]*(output_2[i]-c1[i]);
	}

	//compute a - bd^{-1}b*
	std::complex<double_t> tmp[length];
	for(uint32_t i = 0; i < length; ++i) {
		tmp[i] = f_a[i]- bd_1[i]*std::conj(f_b[i]);
	}
	delete[] bd_1;

	//q0 <-- sampleFz(a - bd^{-1}b*, c0)
	P_sampleFz(output_1, tmp, c0, length);
}

void Gaussian::GaussianSampling::encapsulation::P_sampleFz(std::complex<double_t> *output, const std::complex<double_t> *f, const std::complex<double_t> *c, const uint32_t length) noexcept {
	if (length == 1) {
		*output = (double_t) GenerateIntegerKarney(sqrt(fabs(f[0].real())), c[0].real());
	}
	else {
		const uint32_t length_2 = length >> 1;
		// let f(x) = f0(x^2) + x.f1(x^2)
		std::complex<double_t> f0[length_2];
		std::complex<double_t> f1[length_2];
		fastFft.ffsplit(f0, f1, f, length);

		// (q0, q1) Sample2Z(f0, f1, f0, c)
		std::complex<double_t> * q0 = new std::complex<double_t>[length_2];
		std::complex<double_t> * q1 = new std::complex<double_t>[length_2];
		P_sample2z(q0, q1, f0, f1, f0, c, length_2);

		// let q(x) = q0(x^2) + x.q1(x^2)
		fastFft.ffmerge(output, q0, q1, length_2);
		delete[] q0;
		delete[] q1;
	}
}

/* Private Interface End ----------------------------------------------------------------------------- */

/* Public Interface ---------------------------------------------------------------------------------- */
Gaussian::GaussianSampling::GaussianSampling(void) noexcept {}
Gaussian::GaussianSampling::GaussianSampling(const double_t s, Poly_t R[][k]) noexcept : \
impl(new encapsulation {.s=s}) {
    for(uint32_t i = 0; i < w; ++i)
	    for(uint32_t j = 0; j< k; ++j)
	        impl->R[i][j]={R[i][j]};
				
	impl->sigma = r*sqrt(base*base+1.0)/(base+1.0);
		
	// get the base-ary decomposition sequence of module q
	getBits(impl->modulus, q, k);
	
	// Initialize random generator
	std::random_device rd;
	impl->randomGenerator_t.reset(new std::mt19937(rd()));
	impl->randomGenerator = impl->randomGenerator_t.get();
	
	// compute G_l, G_h and G_d
	impl->G_l[0] = sqrt(base*(1.0 + 1.0/k) + 1.0);
	impl->G_h[0] = 0;
	for(uint32_t i = 1; i < k; ++i) {
		impl->G_l[i] = sqrt(base*(1.0 + 1.0/(k - i)));
		impl->G_h[i] = sqrt(base*(1.0 - 1.0/(k - i + 1)));
	}
	impl->G_d[0] = (double_t)impl->modulus[0]/base;
	for(uint32_t i = 1; i < k; ++i) {
		impl->G_d[i] = (impl->G_d[i - 1] + impl->modulus[i])/(double_t) base;
	}

	// compute P_z0, P_z1 and P_z2
	const double_t s_2 = s*s;
	const double_t r_2 = r*r;
	impl->P_z0 = (s_2 - r_2)/(r_2*s_2);
	impl->P_z1 = sqrt(s_2 - r_2);
	impl->P_z2 = (-r_2)/(s_2 - r_2);

	// compute P_a, P_b and P_d
	Poly_t rr{0};
	Poly_t re{0};
	Poly_t ee{0};
	for(uint32_t i = 0; i < k; ++i) {
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
	    if(mpz_cmp_ui(rrCoefs[i],q/2)>0)
	        mpz_sub_ui(rrCoefs[i],rrCoefs[i],q);
		if(mpz_cmp_ui(reCoefs[i],q/2)>0)
	        mpz_sub_ui(reCoefs[i],reCoefs[i],q);
		if(mpz_cmp_ui(eeCoefs[i],q/2)>0)
	        mpz_sub_ui(eeCoefs[i],eeCoefs[i],q);
	}
	std::complex<double_t> tmp_a[coefsSize];
	std::complex<double_t> tmp_b[coefsSize];
	std::complex<double_t> tmp_d[coefsSize];
	for(uint32_t i = 0; i < coefsSize; ++i) {
		tmp_a[i] = std::complex<double_t>(s_2 - impl->P_z0*mpz_get_ui(rrCoefs[i]), 0);
		tmp_b[i] = std::complex<double_t>(-impl->P_z0*mpz_get_ui(reCoefs[i]), 0);
		tmp_d[i] = std::complex<double_t>(s_2 - impl->P_z0*mpz_get_ui(eeCoefs[i]), 0);
	}

	// compute the fft for P_a, P_b and P_d
	impl->fastFft.fftForward(impl->P_a, tmp_a);
	impl->fastFft.fftForward(impl->P_b, tmp_b);
	impl->fastFft.fftForward(impl->P_d, tmp_d);
}

Gaussian::GaussianSampling::GaussianSampling(GaussianSampling && other) noexcept {
	impl.reset(other.impl.get());
}

Gaussian::GaussianSampling::~GaussianSampling(void) noexcept {
}

Gaussian::GaussianSampling & Gaussian::GaussianSampling::operator= (GaussianSampling && other) noexcept {
	impl.reset(other.impl.get());
	return *this;
}

/*SampleGPoly is used to generate the preimage z with Gz = u and convert z from vector to polynomial.
*/
void Gaussian::GaussianSampling::sampleGPoly(Poly_t * polys, const uint64_t * coset)  noexcept {
	uint64_t * samplingz = new uint64_t[n*k];
	uint64_t * samplingzInd  = samplingz;
	int64_t samplingzIndex[k];
	for(uint32_t i = 0; i < n; ++i) {
		impl->G_sampleG(samplingzIndex, coset[i]);

		for(uint32_t j = 0; j < k; ++j){
			if(samplingzIndex[j]<0)
				*samplingzInd = samplingzIndex[j] + q;
			else
			    *samplingzInd = samplingzIndex[j];
		samplingzInd++;}
    }
	// transpose the samplingz matrix
	uint64_t  * samplingzT = new uint64_t [n*k];
	const uint64_t * samplingzInd2  = samplingz;
	NFL_POLY_COEF_TYPE * samplingzTInd = samplingzT;
	for(uint32_t i = 0; i < k; ++i) {
		for(uint32_t j = 0; j < n; ++j) {
			*samplingzTInd++ = (NFL_POLY_COEF_TYPE) *(samplingzInd2 + k*j);
		}
		++samplingzInd2;
	}
	delete[] samplingz;

	//Transforming the vector samplingzT of n*k long into a polynomial vector of k long.
	samplingzTInd = samplingzT;
	for(uint32_t i = 0; i < k; ++i) {
	    polys[i].set(samplingzTInd, samplingzTInd + n,false);
		polys[i].ntt_pow_phi();
		samplingzTInd += n;
	}
	delete[] samplingzT;
}

//Generating a perturbation vector for SamplePre
void Gaussian::GaussianSampling::sampleP(Poly_t *pSamples, const Poly_t A[]) noexcept {
	    impl->P_samplePz(pSamples);
		
        //The last polynomial of vector pSamples is the product of A and p.
		pSamples[m].set(0);
		for(uint32_t i=0;i<m;i++) {
		    pSamples[m] = pSamples[m] + A[i]*pSamples[i];
        }
}

//A stack for storing the perturbation vectors for SamplePre. The stack was not executed because the generation algorithm of perturbation vector was not run offline.
Poly_t * Gaussian::GaussianSampling::samplePz(void) const noexcept {
	std::lock_guard<std::mutex> lock(impl->P_precomputed_mutex);
	Poly_t * result = impl->P_precomputedPerturbations.top();
	impl->P_precomputedPerturbations.pop();
	return result;
}
