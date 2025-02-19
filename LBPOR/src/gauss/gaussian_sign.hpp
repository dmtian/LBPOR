#ifndef GAUSSIAN_SIGN_HPP
#define GAUSSIAN_SIGN_HPP
#include<cstdint>
#include<cmath>
#include<memory>
#include<algorithm>
#include "param.hpp"
namespace Gaussian {
	class GaussianSampling_sign {

		public:

			/** 
			 * @Default Constructor
			 */
			GaussianSampling_sign(void) noexcept;

			/** 
			 * @Constructor
			 */
			GaussianSampling_sign(const double_t s, Poly_s R[][k_sign]) noexcept;

			/** 
			 * @Move Constructor
			 */
			GaussianSampling_sign(GaussianSampling_sign && other) noexcept;

			/** 
			 * @Destructor
			 */
			~GaussianSampling_sign(void) noexcept;

			/** 
			 * @Move assignment
			 */
			GaussianSampling_sign & operator=(GaussianSampling_sign && other) noexcept;
			
			//The preimage sampling algorithm of g-trapdoor
			void sampleGPoly(Poly_s * polys,  const uint64_t * coset) const noexcept;

            //Perturbation vector generation algorithm			
		    Poly_s * samplePz(void) const noexcept;
            void sampleP(Poly_s *pSamples, const Poly_s A[])  noexcept;
		private:
			struct encapsulation;
			std::unique_ptr<encapsulation> impl;
	};
}

#endif
