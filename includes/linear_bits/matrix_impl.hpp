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
	Matrix<T>& Matrix<T>::operator*= (const Matrix& rhs){
		
		if(n_row != rhs.n_col){
			throw std::invalid_argument("Invalid matrix multiplication.");
		}

		std::vector<T> row (n_row);
		for(natural r = 0; r < n_row; ++r){
			for(natural c = 0; c < n_col; ++c){
				T prod {0};
				for (natural j = 0; j < n_row; ++j){
					prod += this->operator()(r, j) * rhs(j, c);
				}
				row[c] = prod;
			}

			for(natural c = 0; c < n_col; ++c){
				this->operator()(r, c) = row[c];
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

		for(size_t i=0; i < t.size(); ++i){
			mem->operator[](i) = t[i];
		}

		auto sw {n_row};
		n_row = n_col;
		n_col = sw;

		return *this;
	}
}
