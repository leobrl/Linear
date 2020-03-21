

namespace linear{

    template<typename T>
    std::ostream& operator<< (std::ostream& os, const Block<T>& b){ 
        for (auto& e: b.buffer){
            os << e;
            if (&e != &b.buffer.back()) os << " ";
        }
        return os;
    };

	template<typename T>
    std::ostream& operator<< (std::ostream& os, const Matrix<T>& mat){

		for (natural row = 0; row < mat.n_row; ++row){		
			for (natural col = 0; col < mat.n_col-1; ++col){
				os << mat(row, col) << " ";
			}
			os << mat(row, mat.n_col-1);

			if(row != (mat.n_col - 1)){
				os << std::endl;
			}
		}
		return os;
	};

}
