#ifndef PRIVATE_KEY_GENERATOR_HPP
#define PRIVATE_KEY_GENERATOR_HPP
#include<cstdint>
#include<cmath>
#include<memory>
#include "param.hpp"
#include<iostream>
#include <chrono>
#include<algorithm>
#include<complex>

namespace POR{

class LBPOR {

public:
	LBPOR(const double_t sigma, const double_t stag, const double_t sigma_sign, const double_t ssign,const uint32_t labmda) noexcept; //Initialize scheme parameters
	~LBPOR(void) noexcept;
	
	//The algorithms of our lattice-based POR scheme.
    void setup(void) noexcept;
    void Store(std::string &String_namename, Poly_t U[], uint16_t ***M, int64_t ***X, const uint32_t l, std::string& T, int64_t X_sign[][n]) const noexcept;
	void Audit(uint32_t I[],int16_t **V, const uint32_t c, const uint32_t l) noexcept;
	void Prove(Poly_t P_x[], Poly_t P_m[], uint32_t I[], int16_t **V, uint32_t c, uint16_t ***M, int64_t ***X) noexcept;
	int32_t Verify(std::string name, std::string T, Poly_t U[], Poly_t P_x[], Poly_t P_m[],uint32_t c,uint32_t I[], int16_t **V,int64_t X_sign[][n]) noexcept;
	
	void setGaussian(void) noexcept; //Initialize SamplePre algorithm

private:
	struct encapsulation;
	std::unique_ptr<encapsulation> impl;
};
}
#endif
