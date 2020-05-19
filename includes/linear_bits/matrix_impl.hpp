#pragma once

#include "matrix.hpp"

namespace linear{

	template<typename T>
	Matrix<T>::Matrix(natural n_row_, natural n_col_) : n_row(n_row_), n_col(n_col_){
		mem = std::make_unique<Block<T>>(n_col*n_row);
	}

	template<typename T>
	Matrix<T>::Matrix(const Matrix& rhs) : n_row(rhs.n_row), n_col(rhs.n_col){
		mem = std::make_unique<Block<T>>(*rhs.mem);
	}
	
	template<typename T> 
	Matrix<T>::Matrix(const std::vector<T>& diag) : n_row(diag.size()), n_col(diag.size()) {

		mem = std::make_unique<Block<T>>(n_row*n_col);

		for (natural r = 0; r < n_row; ++r){
			mem->operator[](r + r*n_col) = diag[r];
		}
	}

	template<typename T>
	Matrix<T>::Matrix(const natural n_row_, const natural n_col_, const std::vector<T>& elem) : n_row(n_row_), n_col(n_col_)
	{
		if (n_row*n_col != elem.size()) { 
			throw std::invalid_argument("Number of elements is inconsistent with dimensions provided.");
		}

		mem = std::make_unique<Block<T>>(elem);
	}

	template<typename T>
	Matrix<T>::Matrix(Matrix&& rhs) noexcept : n_row(rhs.n_row), n_col(rhs.n_col), mem(std::move(rhs.mem)){
		rhs.n_row = 0;
		rhs.n_col = 0;
	}

	template<typename T>
	typename Matrix<T>::row_iterator
	Matrix<T>::begin(){
		return typename Matrix<T>::row_iterator(*this);
	}

	template<typename T>
	typename Matrix<T>::row_iterator
	Matrix<T>::end(){
		return typename Matrix<T>::row_iterator(*this, n_row, 0); // one elemen after end
	}

	template<typename T>
	Matrix<T>::row_iterator::row_iterator(Matrix<T>& m) : mat(m){
		current_row = natural(0);
		current_col = natural(0);
	}

	template<typename T>
	Matrix<T>::row_iterator::row_iterator(Matrix<T>& m, natural row, natural col) : mat(m), current_row(row), current_col(col){
	}

	template<typename T>
	T& Matrix<T>::row_iterator::operator* (){
		return mat.operator()(current_row, current_col);
	}

	template<typename T>
	typename Matrix<T>::row_iterator& 
	Matrix<T>::row_iterator::operator++ (){

		current_col++;
		if(current_col >= mat.n_col){
			current_row++;
			current_col = 0;
		}

		return *this;
	}

	template<typename T>
	typename Matrix<T>::row_iterator 
	Matrix<T>::row_iterator::operator++ (int){

		typename Matrix<T>::row_iterator temp {*this}; 

		++(*this);

		return temp;
	}

	template<typename T>
	typename Matrix<T>::row_iterator& 
	Matrix<T>::row_iterator::operator-- (){

		if(current_col > 0){
			current_col--;
		} else{
			current_row--;
			current_col = mat.n_col - 1;
		}
		return *this;
	}

	template<typename T>
	typename Matrix<T>::row_iterator
	Matrix<T>::row_iterator::operator-- (int){

		typename Matrix<T>::row_iterator temp {*this}; 

		--(*this);

		return temp;
	}

	template<typename T>
	bool Matrix<T>::row_iterator::operator!=(const Matrix<T>::row_iterator& it) const{
		return( (current_row != it.current_row) || (current_col != it.current_col) );
	}

	template<typename T>
	bool Matrix<T>::row_iterator::operator==(const Matrix<T>::row_iterator& it) const{
		return((current_row == it.current_row) & (current_col == it.current_col));
	}

	template<typename T>
	typename Matrix<T>::col_iterator
	Matrix<T>::col_begin(){
		return typename Matrix<T>::col_iterator(*this);
	}

	template<typename T>
	typename Matrix<T>::col_iterator
	Matrix<T>::col_end(){
		return typename Matrix<T>::col_iterator(*this, 0, n_col); // one elemen after end
	}

	template<typename T>
	Matrix<T>::col_iterator::col_iterator(Matrix<T>& m) : mat(m){
		current_row = natural(0);
		current_col = natural(0);
	}

	template<typename T>
	Matrix<T>::col_iterator::col_iterator(Matrix<T>& m, natural row, natural col) : mat(m), current_row(row), current_col(col){
	}

	template<typename T>
	T& Matrix<T>::col_iterator::operator* (){
		return mat.operator()(current_row, current_col);
	}

	template<typename T>
	typename Matrix<T>::col_iterator& 
	Matrix<T>::col_iterator::operator++ (){

		current_row++;
		if(current_row >= mat.n_row){
			current_col++;
			current_row = 0;
		}

		return *this;
	}

	template<typename T>
	typename Matrix<T>::col_iterator 
	Matrix<T>::col_iterator::operator++ (int){

		typename Matrix<T>::col_iterator temp {*this}; 

		++(*this);

		return temp;
	}

	template<typename T>
	typename Matrix<T>::col_iterator& 
	Matrix<T>::col_iterator::operator-- (){

		if(current_row > 0){
			current_row--;
		} else{
			current_col--;
			current_row = mat.n_row - 1;
		}

		return *this;
	}

	template<typename T>
	typename Matrix<T>::col_iterator
	Matrix<T>::col_iterator::operator-- (int){

		typename Matrix<T>::col_iterator temp {*this}; 
		--(*this);

		return temp;
	}

	template<typename T>
	bool Matrix<T>::col_iterator::operator!=(const Matrix<T>::col_iterator& it) const{
		return( (current_row != it.current_row) || (current_col != it.current_col) );
	}

	template<typename T>
	bool Matrix<T>::col_iterator::operator==(const Matrix<T>::col_iterator& it) const{
		return((current_row == it.current_row) & (current_col == it.current_col));
	}

	template<typename T>
	Matrix<T>& Matrix<T>::operator*= (Matrix& rhs){
		
		/*
			in place matrix multiplication. Rhs is transposed to reduce 
			cache misses.
		*/

		if((n_row != rhs.n_row) | (n_col != rhs.n_col)){
			throw std::invalid_argument("Invalid matrix multiplication.");
		}
		
		rhs.transpose();

		std::vector<T> row (n_row);
		for(natural r = 0; r < n_row; ++r){
			for(natural c = 0; c < n_col; ++c){
				T prod {0};
				for (natural j = 0; j < n_row; ++j){
					prod += this->operator()(r, j) * rhs(c, j);
				}
				row[c] = prod;
			}

			for(natural c = 0; c < n_col; ++c){
				this->operator()(r, c) = row[c];
			}
		}

		rhs.transpose();

		return *this;
	}

	template<typename T>
	Matrix<T>& Matrix<T>::operator+= (const Matrix& rhs){
		if((n_row != rhs.n_row) | (n_col != rhs.n_col)){
			throw std::invalid_argument("Invalid matrix dimensions.");
		}

		for(natural r = 0; r < n_row; ++r){
			for(natural c = 0; c < n_col; ++c){
				this->operator()(r, c) += rhs.operator()(r, c);
			}
		}

		return *this;
	}

	template<typename T>
	Matrix<T>& Matrix<T>::operator-= (const Matrix& rhs){
		if((n_row != rhs.n_row) | (n_col != rhs.n_col)){
			throw std::invalid_argument("Invalid matrix dimensions.");
		}

		for(natural r = 0; r < n_row; ++r){
			for(natural c = 0; c < n_col; ++c){
				this->operator()(r, c) -= rhs.operator()(r, c);	
			}
		}

		return *this;
	}

	template<typename T>
	Matrix<T>& Matrix<T>::transpose(){
			
		std::vector<T> t;
		for(natural c = 0; c < n_col; ++c){
			for(natural r = 0; r < n_row; ++r){
				t.push_back(mem->operator[](c + n_col*r));
			}
		}

		mem->fill(t);

		auto sw {n_row};
		n_row = n_col;
		n_col = sw;

		return *this;
	}

	template<typename T>
	Matrix<T> Matrix<T>::operator* (const Matrix<T>& lhs){
		return multiply_tiled(lhs);
	}

	template<typename T>
	Matrix<T> Matrix<T>::multiply_naive(const Matrix<T>& rhs){
		natural nr = n_row;
		natural nc = rhs.ncol();
		
		auto res = Matrix<T>(nr, nc);

		natural r, c, i;
		for( r = 0; r < nr; ++r){
			for( c = 0; c < nc; ++c){
				for( i=0; i< n_col; ++i){
					res(r, c) += this->operator()(r, i) * rhs(i, c);
				}
			}
		}
		return res;
	}

	
	template<typename T>
	Matrix<T> Matrix<T>::multiply_tiled(const Matrix<T>& rhs){
		natural nr = n_row;
		natural nc = rhs.ncol();
		
		auto res = Matrix<T>(nr, nc);

		// tile sizes
		natural I_b = 64;
		natural J_b = 32;
		natural K_b = 8;

		// tiling
		natural I, J, K;
		natural i, j, k;
		for( I = 0; I < nr; I += I_b){
			natural ni{std::min(I+I_b, nr)};
			natural NI = ni & 1? ni-1: ni; // one extra loop if ni is odd
			
			for( J = 0; J < nc; J += J_b){
				natural nj{std::min(J+J_b, nc)};
				natural NJ = nj & 1? nj-1 : nj; // one extra loop if ni is odd

				for( K = 0; K < n_col; K += K_b){		
					natural nk{std::min(K+K_b, n_col)};
					
					// tile again to fill memory					 
					for( i = I; i < NI; i += 2){
						for( j = J; j < NJ; j += 2){							
							for( k = K; k < nk; ++k){
								res(i, j) += this->operator()(i, k) * rhs(k, j);
								res(i+1, j) += this->operator()(i+1, k) * rhs(k, j);
								res(i, j+1) += this->operator()(i, k) * rhs(k, j+1);								
								res(i+1, j+1) += this->operator()(i+1, k) * rhs(k, j+1);
							}						
						}
					}
					
					// calculate extra row if rhs has odd number of rows
					for(i = NI; i < ni; i += 1){
						// if NJ != nj it will be updated in next loop
						for( j = J; j < NJ; j += 1){ 
							for( k = K; k < nk; ++k){
								res(i, j) += this->operator()(i, k) * rhs(k, j);
							}
						}
					}

					// calculate extra column if lhs has odd number of columns
					for( j = NJ; j < nj; j += 1){
						for( i = I; i < ni; i += 1){
							for( k = K; k < nk; ++k){
								res(i, j) += this->operator()(i, k) * rhs(k, j);
							}
						}
					}
				}

			}
		}
		

		return res;
	}
	
	/* Template specialization */

	template<>
	inline Matrix<double> Matrix<double>::multiply_tiled(const Matrix<double>& rhs){
		natural nr = n_row;
		natural nc = rhs.ncol();
		
		auto res = Matrix<double>(nr, nc);

		// tile sizes
		natural I_b = 256;
		natural J_b = 256;
		natural K_b = 16;

		natural i_b = 4;
		natural j_b = 4;
		natural i_b_shift = 2;
		natural j_b_shift = 2;


		__m128d a0, a1, a2, a3;
		__m128d prod;
		__m128d b01, b23;
		__m128d ab0, ab1, ab2, ab3;
		__m128d ab4, ab5, ab6, ab7;

		double partial_input[8] __attribute__ ((aligned (16)));
		double partial_output[16] __attribute__ ((aligned (16)));
		
		// tiling
		natural I, J, K;
		natural i, j, k;
		for( I = 0; I < nr; I += I_b){
			natural ni{std::min(I+I_b, nr)};
			natural ni_ = ((ni >> i_b_shift) << i_b_shift);
			natural NI = (ni_ == ni)? ni : ni_; // extra loop needed if ni is odd
			
			for( J = 0; J < nc; J += J_b){
				natural nj{std::min(J+J_b, nc)};
				natural nj_ = ((nj >> j_b_shift) << j_b_shift);
				natural NJ = (nj_ == nj)? nj : nj_; // extra loop needed if ni is odd

				for( K = 0; K < n_col; K += K_b){		
					natural nk = {std::min(K+K_b, n_col)};

					// tile again to fill memory					 
					for( i = I; i < NI; i += i_b){

						for( j = J; j < NJ; j += j_b){

							ab0 = _mm_setzero_pd();
							ab1 = _mm_setzero_pd();
							ab2 = _mm_setzero_pd();
							ab3 = _mm_setzero_pd();
							ab4 = _mm_setzero_pd();
							ab5 = _mm_setzero_pd();
							ab6 = _mm_setzero_pd();
							ab7 = _mm_setzero_pd();

							for( k = K; k < nk; ++k){

								partial_input[0] = this->operator()(i, k);
								partial_input[1] = this->operator()(i+1, k);
								partial_input[2] = rhs(k, j);
								partial_input[3] = rhs(k, j+1);
								partial_input[4] = rhs(k, j+2);
								partial_input[5] = rhs(k, j+3);
								partial_input[6] = this->operator()(i+2, k);
								partial_input[7] = this->operator()(i+3, k);

								a0 = _mm_load_pd1(partial_input);
								a1 = _mm_load_pd1(partial_input+1);
								
								b01 = _mm_load_pd(partial_input+2);		
								b23 = _mm_load_pd(partial_input+4);

								a2 = _mm_load_pd1(partial_input+6);
								a3 = _mm_load_pd1(partial_input+7);

								// a(k, 0) *b(k,:)
								prod = _mm_mul_pd(a0, b01);
								ab0 = _mm_add_pd(ab0, prod);

								prod = _mm_mul_pd(a0, b23);
								ab1 = _mm_add_pd(ab1, prod);
								
								// a(k, 1) *b(k,:)
								prod = _mm_mul_pd(a1, b01);
								ab2 = _mm_add_pd(ab2, prod);

								prod = _mm_mul_pd(a1, b23);
								ab3 = _mm_add_pd(ab3, prod);

								//////////////////////////////

								// a(k, 0) *b(k,:)
								prod = _mm_mul_pd(a2, b01);
								ab4 = _mm_add_pd(ab4, prod);

								prod = _mm_mul_pd(a2, b23);
								ab5 = _mm_add_pd(ab5, prod);
								
								// a(k, 1) *b(k,:)
								prod = _mm_mul_pd(a3, b01);
								ab6 = _mm_add_pd(ab6, prod);

								prod = _mm_mul_pd(a3, b23);
								ab7 = _mm_add_pd(ab7, prod);
							}			
						
							_mm_store_pd(partial_output, ab0);
							_mm_store_pd(partial_output+2, ab1);
							_mm_store_pd(partial_output+4, ab2);
							_mm_store_pd(partial_output+6, ab3);
							_mm_store_pd(partial_output+8, ab4);
							_mm_store_pd(partial_output+10, ab5);
							_mm_store_pd(partial_output+12, ab6);
							_mm_store_pd(partial_output+14, ab7);	

							
							res(i, j) += partial_output[0];
							res(i, j+1) += partial_output[1];
							res(i, j+2) += partial_output[2];
							res(i, j+3) += partial_output[3];
							
							res(i+1, j) += partial_output[4];
							res(i+1, j+1) += partial_output[5];
							res(i+1, j+2) += partial_output[6];
							res(i+1, j+3) += partial_output[7];

							res(i+2, j) += partial_output[8];
							res(i+2, j+1) += partial_output[9];
							res(i+2, j+2) += partial_output[10];
							res(i+2, j+3) += partial_output[11];
							
							res(i+3, j) += partial_output[12];
							res(i+3, j+1) += partial_output[13];
							res(i+3, j+2) += partial_output[14];
							res(i+3, j+3) += partial_output[15];
						} 							
					}	
				
					// calculate extra row if rhs has odd number of rows
					for(i = NI; i < ni; i += 1){
						// if NJ != nj it will be updated in next loop
						for( j = J; j < NJ; j += 1){ 
							for( k = K; k < nk; ++k){
								res(i, j) += this->operator()(i, k) * rhs(k, j);
							}
						}
					}

					// calculate extra column if lhs has odd number of columns
					for( j = NJ; j < nj; j += 1){
						for( i = I; i < ni; i += 1){
							for( k = K; k < nk; ++k){
								res(i, j) += this->operator()(i, k) * rhs(k, j);
							}
						}
					}
				}

			}
		}
		return res;
	}
	
	/**/
	template<typename T>
	Matrix<T>& Matrix<T>::operator>> (const std::vector<T> rhs){

		if((n_row*n_col) != rhs.size()){
			throw std::invalid_argument("Invalid vector dimensions.");
		}
		mem->fill(rhs);

		return *this;
	}
}
