#ifndef LBPOR_HPP
#define LBPOR_HPP
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
	
    void setup(void) noexcept;
    void Store(char *name, int8_t ***M, int64_t ***X, const uint32_t l, int32_t X_sign[][n_sign]) const noexcept;
	void Audit(uint32_t I[],int16_t **V, const uint32_t c, const uint32_t l) noexcept;
	void Prove(Poly_t P_x[], Poly_t P_m[], uint32_t I[], int16_t **V, uint32_t c, int8_t ***M, int64_t ***X) noexcept;
	int32_t Verify(char *name, uint32_t l, Poly_t P_x[], Poly_t P_m[],uint32_t c,uint32_t I[], int16_t **V,int32_t X_sign[][n_sign]) noexcept;
	
	void setGaussian(void) noexcept;

private:
	struct encapsulation;
	std::unique_ptr<encapsulation> impl;
};
}
#endif
