#pragma once

#include "block.hpp"

namespace linear{
    
    template<typename T>
    linear::Block<T>::Block(const size_t n) : buffer(std::vector<T>(n)) {
        buffer.shrink_to_fit();
    };

    template<typename T>
    linear::Block<T>::Block(const std::vector<T>& v) : buffer(v){
        buffer.shrink_to_fit();
    };

    template<typename T>
    void linear::Block<T>::fill(const T& value){
        std::fill(buffer.begin(), buffer.end(), value);
    };
}
