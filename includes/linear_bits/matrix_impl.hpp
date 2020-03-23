#pragma once

#include "matrix.hpp"

namespace linear{

	template<typename T>
	Matrix<T>::Matrix(natural n_row, natural n_col) : n_row(n_row), n_col(n_col){
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
	Matrix<T>::Matrix(const natural n_row, const natural n_col, const std::vector<T>& elem) : n_row(n_row), n_col(n_col)
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

		current_col--;
		if(current_col < 0){
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

}