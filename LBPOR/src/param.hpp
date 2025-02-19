#ifndef PARAM_HPP
#define PARAM_HPP

#include<cstdint>
#include<nfl.hpp>

/*********** Selection of scheme parameters for 128-bit security level ***********/

#define k 16 //k=log_base(q)
#define qlog 62 //log_2(q)
#define base 16 //base of G-trapdoor
#define logbase 4 //log_2(base)

#define n 2048UL //Security parameter
#define w 2 
#define m 18 //m=w+k
#define d 2
#define knum 1 //l_1 norm of a challenge weight
#define p 262147 //Data value range
#define q 4611686018326724609ULL //modulus
#define pi 3.14159265
#define O 0.39894228 //O=1/sqrt(2*pi)
#define NFL_POLY_COEF_TYPE uint64_t
#define r 4.5054//5.5578 //r = w(sqrt(log_2(n)))=sqrt(ln(2*n*(1+1/ceta))/pi) for ceta =2^-80

//Parameter definition of signature scheme
#define k_sign  16//9
#define qlog_sign 62
#define base_sign 16//128
#define logbase_sign 4//7
#define m_sign  18//11
#define q_sign  4611686018326724609ULL 
#define r_sign  4.5054//5.5578

//Definition of polynomial data type in POR scheme
using poly_tt = nfl::poly_from_modulus<NFL_POLY_COEF_TYPE, n, qlog>;
using Poly_t  = nfl::poly_p<typename poly_tt::value_type, poly_tt::degree, poly_tt::nmoduli>;

//Definition of polynomial data type in signature scheme
using poly_ss= nfl::poly_from_modulus<uint64_t, n, qlog_sign>;
using Poly_s  = nfl::poly_p<typename poly_ss::value_type, poly_ss::degree, poly_ss::nmoduli>;

//Definition of polynomial data type with 30-bit modulus
#define q_small 1073479681UL
#define qlog_small 30
using poly_32s= nfl::poly_from_modulus<uint32_t, n, qlog_small>;
using Poly_32  = nfl::poly_p<typename poly_32s::value_type, poly_32s::degree, poly_32s::nmoduli>;

/******** Define the parameters for the Gaussians **************/
//Definition of Gaussians in POR scheme
using Gauss_t = nfl::gaussian<uint16_t, NFL_POLY_COEF_TYPE, 2>;
using fastGauss_t = nfl::FastGaussianNoise<uint16_t, NFL_POLY_COEF_TYPE, 2>;

//Definition of Gaussians in signature scheme
using Gauss_s = nfl::gaussian<uint16_t, uint64_t, 2>;
using fastGauss_s = nfl::FastGaussianNoise<uint16_t, uint64_t, 2>;

#endif
