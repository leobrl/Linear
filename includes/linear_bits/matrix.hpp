#pragma once

#include "block.hpp"
#include "typedef_elem.hpp"

#include <memory>

namespace linear{
	
	template<typename T> class Matrix;
	template <typename T> std::ostream& operator<<( std::ostream&, const Matrix<T>& );

	template<typename T>
	class Matrix{
		typedef std::unique_ptr<Block<T>> block;
		#define MAKE_POINTER_TO_BLOCK std::make_unique<Block<T>>
		
		public:
			const natural n_row;
			const natural n_col;
		
		private:
			block mem; 

		public:

			inline Matrix() : n_row(0), n_col(0) {};
			inline ~Matrix()	= default;
			
			explicit Matrix(natural, natural);
			explicit Matrix(const Matrix&);
			explicit Matrix(const std::vector<T>&);
			
			Matrix(Matrix&&)	= default;

			Matrix& operator= (const Matrix&)	= delete;
			Matrix& operator= (Matrix&&)		= delete;

			inline const T&		operator() (natural r, natural c) const {return mem->operator[](r + n_col*c);}
			inline T& 			operator() (natural r, natural c) 		{return mem->operator[](r + n_col*c);}
			
			// Friends
			friend std::ostream& operator<< <T> (std::ostream&, const Matrix<T>&);
	};
}
