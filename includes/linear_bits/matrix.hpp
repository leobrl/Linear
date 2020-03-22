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
			explicit Matrix(const natural, const natural, const std::vector<T>&);
			
			Matrix(Matrix&&)	= default;

			Matrix& operator= (const Matrix&)	= delete;
			Matrix& operator= (Matrix&&)		= delete;

			inline const T&		operator() (natural r, natural c) const {return mem->operator[](c + n_col*r);} // row major
			inline T& 			operator() (natural r, natural c) 		{return mem->operator[](c + n_col*r);} // row major

			class row_iterator{
				
				public:
					typedef std::bidirectional_iterator_tag iterator_category;
					typedef T 								value_type;
					typedef T*                             	pointer;
					typedef T&                             	reference;

				private:
					Matrix& mat;
					natural current_row;
					natural current_col;

				public:
					explicit 	row_iterator(Matrix<T>&);
					explicit 	row_iterator(Matrix<T>&, natural, natural);
					inline 		row_iterator(const row_iterator&) = default;

					reference operator* ();

					row_iterator& operator++ ();
					row_iterator& operator-- ();

					bool operator!= (const row_iterator&) const;
					bool operator== (const row_iterator&) const;
			};

			row_iterator begin();
			row_iterator end();

			// Friends
			friend std::ostream& operator<< <T> (std::ostream&, const Matrix<T>&);
	};
}
