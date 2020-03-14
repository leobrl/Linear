#pragma once

#include<memory>

namespace linear{

    template<typename T>
    class Block{
        
        protected:

            std::vector<T> buffer;

        public:

            inline Block() = default;
            
            explicit Block(const size_t);

            ~Block() = default;

            Block(Block&&) = default;
            
            Block(const Block&) = default;

            Block& operator= (const Block&) = default;

            Block& operator= (Block&&) = default;    

            void fill(const T&);

            inline T operator() (size_t idx) { return buffer[idx]; }

            
    };
}
