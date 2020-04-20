#pragma once

#include "block.hpp"

namespace linear{
    
    template<typename T>
    Block<T>::Block(const size_t n) : buffer(std::vector<T>(n)) {
        buffer.shrink_to_fit();
    };

    template<typename T>
    Block<T>::Block(const std::vector<T>& v) : buffer(v){
        buffer.shrink_to_fit();
    };

    template<typename T>
    void Block<T>::fill(const T& value){
        std::fill(buffer.begin(), buffer.end(), value);
    };

    template<typename T>
    void Block<T>::fill(const std::vector<T>& rhs){

        if(rhs.size() != buffer.size()){
            throw std::invalid_argument("Invalid vector dimensions.");
        }

        std::copy(rhs.begin(), rhs.end(), buffer.begin());
    }
}
