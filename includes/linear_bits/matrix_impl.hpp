#pragma once

#include "matrix.hpp"

namespace linear{

	template<typename T>
	Matrix<T>::Matrix(natural n_row, natural n_col) : n_row(n_row), n_col(n_col){
		mem = MAKE_POINTER_TO_BLOCK(n_col*n_row);
	}

	template<typename T>
	Matrix<T>::Matrix(const Matrix& rhs) : n_row(rhs.n_row), n_col(rhs.n_col){
		mem = MAKE_POINTER_TO_BLOCK(Block<T>(*rhs.mem));
	}
	
	template<typename T> 
	Matrix<T>::Matrix(const std::vector<T>& diag) : n_row(diag.size()), n_col(diag.size()) {

		mem = MAKE_POINTER_TO_BLOCK(n_row*n_col);

		for (natural r = 0; r < n_row; ++r){
			mem->operator[](r + r*n_col) = diag[r];
		}
	}
}