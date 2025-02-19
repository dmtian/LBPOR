#ifndef GAUSSIAN_HPP
#define GAUSSIAN_HPP
#include<cstdint>
#include<cmath>
#include<memory>
#include<algorithm>
#include "param.hpp"
namespace Gaussian {
	class GaussianSampling {

		public:

			/** 
			 * @Default Constructor
			 */
			GaussianSampling(void) noexcept;

			/** 
			 * @Constructor
			 */
			GaussianSampling(const double_t s, Poly_t R[][k]) noexcept;

			/** 
			 * @Move Constructor
			 */
			GaussianSampling(GaussianSampling && other) noexcept;

			/** 
			 * @Destructor
			 */
			~GaussianSampling(void) noexcept;

			/** 
			 * @Move assignment
			 */
			GaussianSampling & operator=(GaussianSampling && other) noexcept;
			
			//The preimage sampling algorithm of g-trapdoor
			void sampleGPoly(Poly_t * polys,  const uint64_t * coset)  noexcept;	
		    
			//Perturbation vector generation algorithm
			Poly_t * samplePz(void) const noexcept;
            void sampleP(Poly_t *pSamples, const Poly_t A[])  noexcept;
		private:
			struct encapsulation;
			std::unique_ptr<encapsulation> impl;
	};
}

#endif
