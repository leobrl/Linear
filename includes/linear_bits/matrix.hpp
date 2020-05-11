#pragma once

#include "block.hpp"
#include "typedef_elem.hpp"

#include <memory>
#include <emmintrin.h>

namespace linear{
	
	template<typename T> class Matrix;
	template <typename T> std::ostream& operator<<( std::ostream&, const Matrix<T>& );

	template<typename T>
	class Matrix{
		
		using block = std::unique_ptr<Block<T>>;

		private:
			natural n_row;
			natural n_col;

			block mem; 

		public:

			inline 			Matrix() noexcept : n_row(0), n_col(0) {}
			inline 			~Matrix() = default;
			
							Matrix(const Matrix&);
							Matrix(Matrix&&) noexcept;
			explicit 		Matrix(natural, natural);							
			explicit 		Matrix(const std::vector<T>&);
			explicit 		Matrix(const natural, const natural, const std::vector<T>&);
			
			Matrix& 		operator= (const Matrix&)	= delete;
			Matrix& 		operator= (Matrix&&)		= delete;

			inline const T&	operator() (natural r, natural c) const {return mem->operator[](c + n_col*r);} // row major
			inline T& 		operator() (natural r, natural c) 		{return mem->operator[](c + n_col*r);} // row major

			Matrix& 		transpose();

			Matrix& 		operator*= 	(Matrix&);
			Matrix& 		operator+= 	(const Matrix&);
			Matrix&			operator-= 	(const Matrix&);
			
			Matrix			operator* 	(const Matrix&);

			Matrix&			operator>>	(const std::vector<T>);
		
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
					row_iterator operator++ (int);
					row_iterator& operator-- ();
					row_iterator operator-- (int);

					bool operator!= (const row_iterator&) const;
					bool operator== (const row_iterator&) const;
			};

			class col_iterator{
				
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
					explicit 	col_iterator(Matrix<T>&);
					explicit 	col_iterator(Matrix<T>&, natural, natural);
					inline 		col_iterator(const col_iterator&) = default;

					reference operator* ();

					col_iterator& operator++ ();
					col_iterator operator++ (int);
					col_iterator& operator-- ();
					col_iterator operator-- (int);

					bool operator!= (const col_iterator&) const;
					bool operator== (const col_iterator&) const;
			};

			row_iterator 	begin();
			row_iterator 	end();

			col_iterator 	col_begin();
			col_iterator 	col_end();

			inline natural	nrow() const {return n_row;}
			inline natural 	ncol() const {return n_col;}

			// Friends
			friend std::ostream& operator<< <T> (std::ostream&, const Matrix<T>&);

		private:
			Matrix 		multiply_naive(const Matrix&);
			Matrix 		multiply_tiled(const Matrix&);
	};
}
