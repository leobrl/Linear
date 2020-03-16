#pragma once

#include<memory>

namespace linear{

    template <typename T> class Block;
    template <typename T> std::ostream& operator<<( std::ostream&, const Block<T>& );

    template<typename T>
    class Block{
        
        protected:

            std::vector<T> buffer;

        public:

            inline Block() = default;
            
            explicit Block(const size_t);

            explicit Block(std::vector<T>);

            ~Block() = default;

            Block(const Block&) = default;

            Block(Block&&) = default;

            Block& operator= (const Block&) = default;

            Block& operator= (Block&&) = default;    

            void fill(const T&);

            inline T& operator[] (size_t idx) { return buffer[idx]; }

            inline const T& operator[] (size_t idx) const { return buffer[idx]; } 

            inline size_t size () const { return buffer.size(); }

            friend std::ostream& operator<< <T> (std::ostream&, const Block<T>&);
    };

    template<typename T>
    std::ostream& operator<< (std::ostream& os, const Block<T>& b){ 
        for (auto& e: b.buffer){
            os << e;
            if (&e != &b.buffer.back()) os << " ";
        }
        return os;
    };
}
