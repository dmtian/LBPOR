#ifndef PARAM_HPP
#define PARAM_HPP

#include<cstdint>
#include<nfl.hpp>

#define k_POR 8
#define qlog 62
#define base_POR 256
#define logbase 8

#define n_POR 2048UL
#define w_POR 2 
#define m_POR 8
#define d_POR 8
#define approx 2

#define k_approx 6
#define knum 14
#define p_POR 128
#define q_POR 4611686018326724609ULL
#define pi 3.14159265
#define O 0.39894228
#define NFL_POLY_COEF_TYPE uint64_t
#define r_POR 5.5578

#define n_sign 1024UL
#define k_sign  6
#define k_approx_sign 3
#define qlog_sign 30
#define approx_sign 3 
#define base_sign 32
#define logbase_sign 5
#define m_sign  5
#define q_sign  1073479681UL 
#define r_sign  5.5379

using poly_tt = nfl::poly_from_modulus<NFL_POLY_COEF_TYPE, n_POR, qlog>;
using Poly_t  = nfl::poly_p<typename poly_tt::value_type, poly_tt::degree, poly_tt::nmoduli>;

using poly_ss= nfl::poly_from_modulus<uint32_t, n_sign, qlog_sign>;
using Poly_s  = nfl::poly_p<typename poly_ss::value_type, poly_ss::degree, poly_ss::nmoduli>;

#define q_small 1073479681UL
#define qlog_small 30
using poly_32s= nfl::poly_from_modulus<uint32_t, n_POR, qlog_small>;
using Poly_32  = nfl::poly_p<typename poly_32s::value_type, poly_32s::degree, poly_32s::nmoduli>;

using Gauss_t = nfl::gaussian<uint16_t, NFL_POLY_COEF_TYPE, 2>;
using fastGauss_t = nfl::FastGaussianNoise<uint16_t, NFL_POLY_COEF_TYPE, 2>;

using Gauss_s = nfl::gaussian<uint16_t, uint32_t, 2>;
using fastGauss_s = nfl::FastGaussianNoise<uint16_t, uint32_t, 2>;

#endif
