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

            inline Block()  = default;
            inline ~Block() = default;

            explicit Block(const size_t);
            explicit Block(std::vector<T>);

            Block(const Block&) = default;
            Block(Block&&)      = default;

            Block& operator= (const Block&) = default;
            Block& operator= (Block&&)      = default;    

            void fill(const T&);
            
            inline const T& operator[] (size_t idx) const   { return buffer[idx]; }
            inline T&       operator[] (size_t idx)         { return buffer[idx]; }

            inline size_t size () const { return buffer.size(); }

            // Friends
            friend std::ostream& operator<< <T> (std::ostream&, const Block<T>&);
    };
}
