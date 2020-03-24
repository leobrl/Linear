#pragma once

#include<memory>

namespace linear{

    template <typename T> class Block;
    template <typename T> std::ostream& operator<<( std::ostream&, const Block<T>& );

    template<typename T>
    class Block{
        
        using raw_memory = std::vector<T>;

        protected:

            raw_memory buffer;

        public:

            inline Block()  = default;
            inline ~Block() = default;

            explicit Block(const size_t);
            explicit Block(const std::vector<T>& );

            Block(const Block&) = default;
            Block(Block&&)      = default;

            Block& operator= (const Block&) = default;
            Block& operator= (Block&&)      = default;    

            void fill(const T&);
            
            inline const T& operator[] (size_t idx) const   { return buffer[idx]; }
            inline T&       operator[] (size_t idx)         { return buffer[idx]; }

            inline size_t size () const { return buffer.size(); }

            inline typename raw_memory::iterator begin()    {return buffer.begin();}
            inline typename raw_memory::iterator end()      {return buffer.end();}

            // Friends
            friend std::ostream& operator<< <T> (std::ostream&, const Block<T>&);
    };
}
